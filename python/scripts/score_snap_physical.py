"""Score paired PHerc1203 snap meshes against the physical recto labels.

This is the frozen scorer for ``PHERC1203_PHYSICAL_SNAP_REPLICATION_PREREG``.
It compares the same vertices and topology in the PR #11 repair-only arm and
the former four-iteration arm.  The reference labels and their provenance are
from 7jycwjmbfn-eng/pherc0139-physical-audit / Villa PR #1382.

Distances are computed independently in every axial L1 label plane, then
sampled trilinearly.  Eligibility and surface-area weights are fixed from the
common pre-snap mesh before either arm is summarized.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import zarr
from scipy.ndimage import distance_transform_edt, map_coordinates

from score_snap_ridge import load_obj, sha256


LABEL_ORIGIN_L1 = np.array([3936.0, 0.0, 0.0], dtype=np.float64)
LABEL_SHAPE = (2016, 3456, 3456)
LABEL_WINDOW = (
    slice(976, 1072),
    slice(1520, 1872),
    slice(1520, 1872),
)
WINDOW_ORIGIN = np.array(
    [part.start for part in LABEL_WINDOW], dtype=np.float64
)
BOX_LOW_L0 = np.array([9856.0, 3072.0, 3072.0], dtype=np.float64)
BOX_HIGH_L0 = np.array([9984.0, 3712.0, 3712.0], dtype=np.float64)
INNER_MARGIN_L0 = 1.0
PRE_DISTANCE_MAX_L1 = 3.0
MIN_VERTICES_PER_CUBE = 100
MIN_CUBES = 20
MIN_CUBES_PER_ROW_COLUMN = 3
EFFECT_MIN_L1 = 0.10
LARGE_WORSENING_L1 = 0.25
BOOTSTRAP_SEED = 20260810
BOOTSTRAP_N = 10_000
EXPECTED_Z = 9856
EXPECTED_YX = tuple(range(3072, 3712, 128))
CUBE_NAME_RE = re.compile(r"^z(\d{5})_y(\d{5})_x(\d{5})$")


def l0_to_window_coordinates(points_l0: np.ndarray) -> np.ndarray:
    """Map z/y/x L0 index coordinates to the frozen label-window indices."""
    points = np.asarray(points_l0, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n, 3) in z/y/x order")
    return (points - 0.5) / 2.0 - LABEL_ORIGIN_L1 - WINDOW_ORIGIN


def all_trilinear_corners_true(mask: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Return whether all eight floor/ceil corners are true for each point."""
    mask = np.asarray(mask, dtype=bool)
    coords = np.asarray(coords, dtype=np.float64)
    if mask.ndim != 3 or coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("mask must be 3-D and coords must have shape (n, 3)")
    finite = np.isfinite(coords).all(axis=1)
    safe = np.where(np.isfinite(coords), coords, -1.0)
    base = np.floor(safe).astype(np.int64)
    inside = finite & np.all(base >= 0, axis=1)
    inside &= np.all(base + 1 < np.asarray(mask.shape), axis=1)
    result = inside.copy()
    rows = np.flatnonzero(inside)
    if not len(rows):
        return result
    b = base[rows]
    corner_ok = np.ones(len(rows), dtype=bool)
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                corner_ok &= mask[b[:, 0] + dz, b[:, 1] + dy, b[:, 2] + dx]
    result[rows] = corner_ok
    return result


def fixed_vertex_area_weights(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Assign one third of every pre-snap triangle area to each incident vertex."""
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices must have shape (n, 3)")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("faces must have shape (m, 3)")
    if len(faces) and (faces.min() < 0 or faces.max() >= len(vertices)):
        raise ValueError("face index outside vertex array")
    tri = vertices[faces]
    area = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1
    )
    weights = np.zeros(len(vertices), dtype=np.float64)
    for corner in range(3):
        np.add.at(weights, faces[:, corner], area / 3.0)
    return weights


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Return a deterministic lower weighted median of finite positive weight."""
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.shape != weights.shape or values.ndim != 1:
        raise ValueError("values and weights must be equal-length 1-D arrays")
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        raise ValueError("weighted median has no finite positive-weight values")
    v = values[valid]
    w = weights[valid]
    order = np.argsort(v, kind="stable")
    v, w = v[order], w[order]
    cutoff = 0.5 * float(w.sum())
    index = int(np.searchsorted(np.cumsum(w), cutoff, side="left"))
    return float(v[min(index, len(v) - 1)])


def bootstrap_median_ci(values: np.ndarray) -> list[float | None]:
    """Frozen paired-cube percentile interval."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("bootstrap values must be one-dimensional")
    if not len(values):
        return [None, None]
    if not np.isfinite(values).all():
        raise ValueError("bootstrap values must all be finite")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(0, len(values), size=(BOOTSTRAP_N, len(values)))
    medians = np.median(values[draws], axis=1)
    return [float(value) for value in np.percentile(medians, [2.5, 97.5])]


def build_axial_recto_distance(labels: np.ndarray) -> np.ndarray:
    """Build a stack of independent 2-D distances to recto-band bit 8."""
    labels = np.asarray(labels, dtype=np.uint8)
    if labels.ndim != 3:
        raise ValueError("labels must be a 3-D uint8 array")
    recto = (labels & np.uint8(8)) != 0
    distances = np.full(labels.shape, np.inf, dtype=np.float32)
    for z_index in range(labels.shape[0]):
        if recto[z_index].any():
            distances[z_index] = distance_transform_edt(
                ~recto[z_index]
            ).astype(np.float32)
    return distances


def sample_distances(distance_stack: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Trilinearly sample the axial distance stack."""
    values = map_coordinates(
        distance_stack,
        np.asarray(coords, dtype=np.float64).T,
        order=1,
        mode="constant",
        cval=np.inf,
        prefilter=False,
    )
    return np.asarray(values, dtype=np.float64)


def read_frozen_label_window(labels_path: Path) -> tuple[np.ndarray, dict]:
    """Read and validate only the preregistered PHerc1203 label window."""
    array = zarr.open_array(str(labels_path), mode="r")
    if tuple(array.shape) != LABEL_SHAPE or np.dtype(array.dtype) != np.uint8:
        raise ValueError(
            f"unexpected label shape/dtype: {array.shape} {array.dtype}"
        )
    attrs = dict(array.attrs)
    if attrs.get("origin_l1") != [3936, 0, 0]:
        raise ValueError("label origin_l1 does not match the frozen reference")
    expected_bits = {
        "valid": 1,
        "material": 2,
        "centerline": 4,
        "recto_band": 8,
        "boundary_poor": 16,
    }
    if attrs.get("bits") != expected_bits:
        raise ValueError("label bit schema does not match the frozen reference")
    labels = np.asarray(array[LABEL_WINDOW], dtype=np.uint8)
    expected_window_shape = tuple(part.stop - part.start for part in LABEL_WINDOW)
    if labels.shape != expected_window_shape:
        raise ValueError(f"short label window: {labels.shape}")
    metadata = {
        "path": labels_path.as_posix(),
        "shape_zyx": list(array.shape),
        "dtype": str(np.dtype(array.dtype)),
        "origin_l1": attrs["origin_l1"],
        "bits": attrs["bits"],
        "registration_heldout_um": attrs.get("registration_heldout_um"),
        "window_local_zyx": [
            [part.start, part.stop] for part in LABEL_WINDOW
        ],
        "zarray_sha256": sha256(labels_path / ".zarray"),
        "zattrs_sha256": sha256(labels_path / ".zattrs"),
    }
    return labels, metadata


def validate_meshes(pre: dict, candidate: dict, former: dict) -> None:
    """Enforce the shared-index contract required for a paired comparison."""
    for name, arm in (("candidate", candidate), ("former", former)):
        if arm["verts"].shape != pre["verts"].shape:
            raise ValueError(f"{name}: vertex count differs from pre")
        if not np.array_equal(arm["faces"], pre["faces"]):
            raise ValueError(f"{name}: topology/face indices differ from pre")
        if not np.array_equal(arm["uv"], pre["uv"]):
            raise ValueError(f"{name}: registered UV differs from pre")
        if arm["cubes"] != pre["cubes"]:
            raise ValueError(f"{name}: cube ranges differ from pre")
    for name, mesh in (("pre", pre), ("candidate", candidate), ("former", former)):
        if not np.isfinite(mesh["verts"]).all():
            raise ValueError(f"{name}: non-finite vertex coordinate")


def cube_coordinates(cube_name: str) -> tuple[int, int, int]:
    match = CUBE_NAME_RE.fullmatch(cube_name)
    if not match:
        raise ValueError(f"malformed cube name: {cube_name}")
    return tuple(int(value) for value in match.groups())


def validate_cube_inventory(cubes: list[tuple[str, int, int]]) -> None:
    names = [item[0] for item in cubes]
    if len(names) != len(set(names)):
        raise ValueError("duplicate cube provenance name")
    expected = {
        f"z{EXPECTED_Z:05d}_y{y:05d}_x{x:05d}"
        for y in EXPECTED_YX
        for x in EXPECTED_YX
    }
    extra = sorted(set(names) - expected)
    if extra:
        raise ValueError(f"cube outside frozen inventory: {extra[0]}")


def _mesh_record(mesh: dict) -> dict:
    return {"path": mesh["path"], "sha256": mesh["sha256"]}


def _row_column_summary(per_cube: list[dict], axis: str) -> dict[str, dict]:
    index = 1 if axis == "y" else 2
    out: dict[str, dict] = {}
    for coordinate in EXPECTED_YX:
        rows = [
            row for row in per_cube
            if row["qualified"] and cube_coordinates(row["cube"])[index] == coordinate
        ]
        effects = np.asarray([row["effect_l1"] for row in rows], dtype=np.float64)
        out[str(coordinate)] = {
            "n_cubes": len(rows),
            "median_effect_l1": float(np.median(effects)) if len(effects) else None,
        }
    return out


def _all_axis_medians(summary: dict[str, dict], positive: bool) -> bool:
    for row in summary.values():
        value = row["median_effect_l1"]
        if row["n_cubes"] < MIN_CUBES_PER_ROW_COLUMN or value is None:
            return False
        if positive and not value > 0:
            return False
        if not positive and not value < 0:
            return False
    return True


def score(args: argparse.Namespace) -> dict:
    labels, label_metadata = read_frozen_label_window(Path(args.labels))
    valid_labels = (labels & np.uint8(1)) != 0
    distance_stack = build_axial_recto_distance(labels)

    pre = load_obj(Path(args.pre))
    candidate = load_obj(Path(args.candidate))
    former = load_obj(Path(args.former))
    validate_meshes(pre, candidate, former)
    validate_cube_inventory(pre["cubes"])

    meshes = {"pre": pre, "iter0": candidate, "iter4": former}
    coords = {
        name: l0_to_window_coordinates(mesh["verts"])
        for name, mesh in meshes.items()
    }
    positions = [mesh["verts"].astype(np.float64) for mesh in meshes.values()]
    weights = fixed_vertex_area_weights(pre["verts"], pre["faces"])
    used = np.zeros(len(pre["verts"]), dtype=bool)
    used[pre["faces"].ravel()] = True
    eligible = used & np.isfinite(weights) & (weights > 0)
    low = BOX_LOW_L0 + INNER_MARGIN_L0
    high = BOX_HIGH_L0 - INNER_MARGIN_L0
    for mesh_positions in positions:
        eligible &= np.all(mesh_positions >= low, axis=1)
        eligible &= np.all(mesh_positions <= high, axis=1)
    for mesh_coords in coords.values():
        eligible &= all_trilinear_corners_true(valid_labels, mesh_coords)

    distances = {
        name: sample_distances(distance_stack, mesh_coords)
        for name, mesh_coords in coords.items()
    }
    for values in distances.values():
        eligible &= np.isfinite(values)
    eligible &= distances["pre"] <= PRE_DISTANCE_MAX_L1

    per_cube: list[dict] = []
    for cube_name, lo, hi in pre["cubes"]:
        mask = eligible[lo:hi]
        local_indices = np.flatnonzero(mask) + lo
        row: dict[str, object] = {
            "cube": cube_name,
            "n_used": int(used[lo:hi].sum()),
            "n_eligible": int(len(local_indices)),
            "used_area_l0_sq": float(weights[lo:hi][used[lo:hi]].sum()),
            "eligible_area_l0_sq": float(weights[local_indices].sum()),
            "qualified": bool(len(local_indices) >= MIN_VERTICES_PER_CUBE),
        }
        if len(local_indices):
            for name in ("pre", "iter0", "iter4"):
                values = distances[name][local_indices]
                row[f"{name}_weighted_median_l1"] = weighted_median(
                    values, weights[local_indices]
                )
                row[f"{name}_unweighted_median_l1"] = float(np.median(values))
            row["effect_l1"] = (
                row["iter4_weighted_median_l1"]
                - row["iter0_weighted_median_l1"]
            )
            row["unweighted_effect_l1"] = (
                row["iter4_unweighted_median_l1"]
                - row["iter0_unweighted_median_l1"]
            )
        else:
            for key in (
                "pre_weighted_median_l1",
                "iter0_weighted_median_l1",
                "iter4_weighted_median_l1",
                "pre_unweighted_median_l1",
                "iter0_unweighted_median_l1",
                "iter4_unweighted_median_l1",
                "effect_l1",
                "unweighted_effect_l1",
            ):
                row[key] = None
        per_cube.append(row)

    qualified = [row for row in per_cube if row["qualified"]]
    effects = np.asarray([row["effect_l1"] for row in qualified], dtype=np.float64)
    unweighted_effects = np.asarray(
        [row["unweighted_effect_l1"] for row in qualified], dtype=np.float64
    )
    row_summary = _row_column_summary(per_cube, "y")
    column_summary = _row_column_summary(per_cube, "x")
    coverage_gate = (
        MIN_CUBES <= len(pre["cubes"]) <= 25
        and len(qualified) >= MIN_CUBES
        and all(
            row["n_cubes"] >= MIN_CUBES_PER_ROW_COLUMN
            for row in [*row_summary.values(), *column_summary.values()]
        )
    )

    median_effect = float(np.median(effects)) if len(effects) else None
    bootstrap95 = bootstrap_median_ci(effects)
    positive_fraction = float((effects > 0).mean()) if len(effects) else None
    negative_fraction = float((effects < 0).mean()) if len(effects) else None
    candidate_large_worse = (
        float((effects < -LARGE_WORSENING_L1).mean()) if len(effects) else None
    )
    former_large_worse = (
        float((effects > LARGE_WORSENING_L1).mean()) if len(effects) else None
    )

    candidate_pass = bool(
        coverage_gate
        and median_effect is not None
        and median_effect >= EFFECT_MIN_L1
        and bootstrap95[0] is not None
        and bootstrap95[0] > 0
        and positive_fraction is not None
        and positive_fraction >= 0.80
        and _all_axis_medians(row_summary, True)
        and _all_axis_medians(column_summary, True)
        and candidate_large_worse is not None
        and candidate_large_worse <= 0.10
    )
    former_pass = bool(
        coverage_gate
        and median_effect is not None
        and median_effect <= -EFFECT_MIN_L1
        and bootstrap95[1] is not None
        and bootstrap95[1] < 0
        and negative_fraction is not None
        and negative_fraction >= 0.80
        and _all_axis_medians(row_summary, False)
        and _all_axis_medians(column_summary, False)
        and former_large_worse is not None
        and former_large_worse <= 0.10
    )
    if candidate_pass and former_pass:
        raise AssertionError("opposed physical-decision gates both passed")
    if not coverage_gate:
        metric_decision = "physical_coverage_gate_fail"
    elif candidate_pass:
        metric_decision = "candidate_iter0_physically_superior"
    elif former_pass:
        metric_decision = "former_iter4_physically_superior"
    else:
        metric_decision = "tie_or_inconclusive_default_decision"

    result = {
        "tool": "score_snap_physical",
        "protocol": {
            "reference_credit": (
                "7jycwjmbfn-eng/pherc0139-physical-audit and Villa PR #1382"
            ),
            "coordinate_mapping": (
                "p_label_local=(p_L0-0.5)/2-[3936,0,0]-[976,1520,1520]"
            ),
            "distance": "per-label-z-plane 2-D EDT to recto_band bit 8",
            "pre_distance_max_l1": PRE_DISTANCE_MAX_L1,
            "weight": "one third pre-snap triangle area per incident vertex",
            "inferential_unit": "provenance cube",
            "effect": "iter4 distance minus iter0 distance; positive favors iter0",
            "bootstrap_seed": BOOTSTRAP_SEED,
            "bootstrap_resamples": BOOTSTRAP_N,
        },
        "label_metadata": label_metadata,
        "inputs": {
            "pre": _mesh_record(pre),
            "candidate_iter0": _mesh_record(candidate),
            "former_iter4": _mesh_record(former),
        },
        "invariants": {
            "vertices": len(pre["verts"]),
            "faces": len(pre["faces"]),
            "cubes_present": len(pre["cubes"]),
            "topology_uv_cube_ranges_identical": True,
            "all_mesh_coordinates_finite": True,
            "common_eligible_vertices": int(eligible.sum()),
            "qualified_cubes": len(qualified),
            "coverage_gate": coverage_gate,
        },
        "gates": {
            "min_vertices_per_cube": MIN_VERTICES_PER_CUBE,
            "min_qualified_cubes": MIN_CUBES,
            "min_cubes_per_row_column": MIN_CUBES_PER_ROW_COLUMN,
            "effect_min_l1": EFFECT_MIN_L1,
            "positive_or_negative_fraction_min": 0.80,
            "large_worsening_l1": LARGE_WORSENING_L1,
            "large_worsening_fraction_max": 0.10,
        },
        "summary": {
            "median_effect_l1": median_effect,
            "bootstrap95_l1": bootstrap95,
            "positive_fraction": positive_fraction,
            "negative_fraction": negative_fraction,
            "candidate_large_worsening_fraction": candidate_large_worse,
            "former_large_worsening_fraction": former_large_worse,
            "median_unweighted_effect_l1": (
                float(np.median(unweighted_effects))
                if len(unweighted_effects) else None
            ),
            "row_y": row_summary,
            "column_x": column_summary,
            "candidate_gate_pass": candidate_pass,
            "former_gate_pass": former_pass,
            "metric_decision": metric_decision,
            "external_pipeline_safety_gates_required": True,
        },
        "per_cube": per_cube,
    }
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--pre", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--former", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = score(args)
    print(
        json.dumps(
            {
                "invariants": result["invariants"],
                "summary": result["summary"],
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
