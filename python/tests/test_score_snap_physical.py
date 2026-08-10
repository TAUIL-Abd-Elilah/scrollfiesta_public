import importlib.util
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import zarr


def _load_script():
    scripts = Path(__file__).resolve().parents[1] / "scripts"
    sys.path.insert(0, str(scripts))
    path = scripts / "score_snap_physical.py"
    spec = importlib.util.spec_from_file_location("score_snap_physical_tested", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_l0_to_label_window_coordinate_mapping_is_center_aligned():
    score = _load_script()
    local = np.array([[0.0, 0.0, 0.0], [47.75, 175.75, 175.75]])
    global_l1 = local + score.LABEL_ORIGIN_L1 + score.WINDOW_ORIGIN
    points_l0 = 2.0 * global_l1 + 0.5
    np.testing.assert_allclose(score.l0_to_window_coordinates(points_l0), local)
    midpoint = np.array([[9920.0, 3392.0, 3392.0]])
    np.testing.assert_allclose(
        score.l0_to_window_coordinates(midpoint),
        np.array([[47.75, 175.75, 175.75]]),
    )


def test_all_trilinear_corners_must_be_valid_and_inside():
    score = _load_script()
    mask = np.ones((4, 4, 4), dtype=bool)
    point = np.array([[1.2, 1.3, 1.4]])
    assert score.all_trilinear_corners_true(mask, point).tolist() == [True]
    mask[2, 2, 2] = False
    assert score.all_trilinear_corners_true(mask, point).tolist() == [False]
    assert score.all_trilinear_corners_true(
        mask, np.array([[-0.1, 1.0, 1.0], [3.0, 1.0, 1.0]])
    ).tolist() == [False, False]


def test_fixed_area_weights_and_weighted_median():
    score = _load_script()
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    faces = np.array([[0, 1, 2]])
    weights = score.fixed_vertex_area_weights(vertices, faces)
    np.testing.assert_allclose(weights, np.full(3, 1.0 / 6.0))
    assert score.weighted_median(
        np.array([0.0, 1.0, 2.0]), np.array([1.0, 5.0, 1.0])
    ) == 1.0
    with pytest.raises(ValueError, match="no finite positive-weight"):
        score.weighted_median(np.array([1.0]), np.array([0.0]))


def test_axial_distance_is_two_dimensional_per_plane():
    score = _load_script()
    labels = np.zeros((2, 5, 5), dtype=np.uint8)
    labels[0, 2, 2] = np.uint8(8)
    labels[1, 0, 0] = np.uint8(8)
    distances = score.build_axial_recto_distance(labels)
    assert distances[0, 2, 3] == pytest.approx(1.0)
    assert distances[1, 2, 3] == pytest.approx(np.sqrt(13.0))
    sampled = score.sample_distances(
        distances, np.array([[0.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
    )
    np.testing.assert_allclose(sampled, [1.0, np.sqrt(13.0)])


def test_bootstrap_interval_is_frozen_and_deterministic():
    score = _load_script()
    effects = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    first = score.bootstrap_median_ci(effects)
    second = score.bootstrap_median_ci(effects)
    assert first == second
    assert first[0] <= 0.3 <= first[1]
    assert score.bootstrap_median_ci(np.array([])) == [None, None]


def test_shared_mesh_contract_and_frozen_cube_inventory():
    score = _load_script()
    mesh = {
        "verts": np.array(
            [[9858.0, 3074.0, 3074.0],
             [9858.0, 3075.0, 3074.0],
             [9858.0, 3074.0, 3075.0]],
            dtype=np.float32,
        ),
        "faces": np.array([[0, 1, 2]], dtype=np.int32),
        "uv": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
        "cubes": [("z09856_y03072_x03072", 0, 3)],
    }
    score.validate_meshes(mesh, deepcopy(mesh), deepcopy(mesh))
    score.validate_cube_inventory(mesh["cubes"])
    changed = deepcopy(mesh)
    changed["uv"][0, 0] = 0.5
    with pytest.raises(ValueError, match="registered UV differs"):
        score.validate_meshes(mesh, changed, deepcopy(mesh))
    with pytest.raises(ValueError, match="outside frozen inventory"):
        score.validate_cube_inventory([("z09856_y02000_x03072", 0, 3)])


def test_label_loader_checks_schema_and_reads_only_window(tmp_path):
    score = _load_script()
    score.LABEL_SHAPE = (4, 5, 6)
    score.LABEL_WINDOW = (slice(1, 3), slice(1, 4), slice(2, 5))
    root = tmp_path / "labels.zarr"
    array = zarr.open_array(
        str(root),
        mode="w",
        shape=score.LABEL_SHAPE,
        chunks=(2, 3, 3),
        dtype="u1",
        zarr_format=2,
    )
    expected = np.arange(18, dtype=np.uint8).reshape(2, 3, 3)
    array[score.LABEL_WINDOW] = expected
    array.attrs.update(
        origin_l1=[3936, 0, 0],
        bits={
            "valid": 1,
            "material": 2,
            "centerline": 4,
            "recto_band": 8,
            "boundary_poor": 16,
        },
        registration_heldout_um=2.38,
    )
    loaded, metadata = score.read_frozen_label_window(root)
    np.testing.assert_array_equal(loaded, expected)
    assert metadata["shape_zyx"] == [4, 5, 6]
    assert metadata["registration_heldout_um"] == 2.38
