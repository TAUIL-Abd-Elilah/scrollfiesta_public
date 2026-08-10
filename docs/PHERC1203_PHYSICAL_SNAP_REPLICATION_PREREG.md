# PHerc1203 physical-truth snap replication preregistration

Status: **frozen before any PHerc1203 voxel, label chunk, mesh, or local
outcome was read** at 2026-08-10T15:38:13+01:00 (Africa/Casablanca).

Before this freeze, only public catalog/Zarr metadata, the physical-audit
README and source code, and its already published whole-scroll aggregate
numbers had been inspected. The `labels1203_L1.tar` asset had not been
downloaded or opened. No result localized to the box below was known.

This is an independent-scroll, independent-physical-reference test of the
default decision in ScrollFiesta PR #11. It explicitly extends, and does not
duplicate or claim authorship of, the PHerc1203 physical labels and binary
volume evaluator published by `7jycwjmbfn-eng/pherc0139-physical-audit` and
proposed in Villa PR #1382. The new question is narrower: which of two
same-topology ScrollFiesta mesh arms is closer to the published physical
recto reference?

## Fixed question and interpretation

On one untouched, metadata-selected PHerc1203 slab, is the PR #11
repair-only default (`snap_recto_iters=0`) closer to an independently derived
physical recto band than the former four-iteration behavior
(`snap_recto_iters=4`, range 3)?

The physical score is the primary decision metric. The already published
CT-ridge score is secondary and mechanistic. If they disagree, the physical
reference controls any practical default claim, and the disagreement must be
reported rather than averaged away.

The physical labels are still a registered, algorithmically derived
reference rather than infallible manual truth. Their published held-out
registration median is 2.38 um. Claims must retain that limitation.

## Frozen code and data

- Candidate implementation: ScrollFiesta PR #11 exact public head
  `4777630e5c81111684d156ef4ffdb5964cedc57d`.
- Reliability patch used for this run: ScrollFiesta PR #12 exact public head
  `f0d9d2e54823e7ba2460725e81290eead8ed6e5e`. The experiment build is the
  PR #12 patch cherry-picked onto PR #11; its resulting commit/tree and
  binary hashes will be recorded before data access. No geometry or scoring
  behavior may otherwise change.
- RAW CT, level 0:
  `s3://vesuvius-challenge-open-data/PHerc1203/volumes/20250820131727-9.362um-1.2m-113keV-masked.zarr`
- m7 surface prediction, level 0:
  `s3://vesuvius-challenge-open-data/PHerc1203/representations/predictions/surfaces/20250820131727-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- Both level-0 arrays have metadata shape `[18977,6844,6844]`; RAW chunks
  are `[128,128,128]` and prediction chunks are `[192,192,192]`.
- Physical reference: release `v1.0`, tag commit
  `6937d846b5cd6cc4fb07dc8eb9770493b91b2128`, asset
  `labels1203_L1.tar`, 515,379,200 bytes, SHA-256
  `32a09f6081342b0f015b258ec577d0296ff23a55892af9785689d8a55bff344c`.
- The contained uint8 label array has L1 origin `[3936,0,0]`, shape
  `[2016,3456,3456]`, and flags `valid=1`, `material=2`, `centerline=4`,
  `recto_band=8`, `boundary_poor=16`.

The fixed L0 box is:

`z [9856,9984), y [3072,3712), x [3072,3712)`

This is `1 x 5 x 5 = 25` cubes of `128^3` voxels. The selection rule uses
metadata only: choose the 128-aligned one-slab z interval nearest the center
of the physical-label z overlap, and the 128-aligned 640-by-640 y/x footprint
nearest the center of the level-0 array. Ties choose the lower coordinate.
The fixed straight axis is the level-0 array center,
`axis_point_zyx=[0,3421.5,3421.5]`, `axis_dir_zyx=[1,0,0]`. The box and axis
may not move after any data value is read.

## Frozen pipeline and arms

1. Build the exact combined source above in Release/x64, run
   `scroll_unroll --selftest` and the relevant common/Python tests, and record
   source, scorer, test, and binary hashes before data access.
2. Carve exactly the fixed RAW and nonzero-to-255 prediction cubes with the
   unmodified public carver. Record source metadata and output manifests.
3. Run `grid_pipeline` with documented defaults, halo 13, one thread per
   cube, and no outcome-dependent retry settings.
4. Run `scroll_whole`, then `--reregister --audit`, with the fixed axis. No
   axis, registration, audit, or placement parameter may be tuned on this
   box.
5. From one common placed input export a pre-snap mesh with `--steps 12` and
   exactly two stage-4 arms, all other settings fixed:
   - candidate: `--steps 124 --snap-recto-iters 0 --snap-recto-range 3`
   - former behavior: `--steps 124 --snap-recto-iters 4 --snap-recto-range 3`
6. Rebuild once and export the candidate without a `--snap-recto-iters`
   flag. It must be byte-identical to explicit iteration 0.

## Frozen primary physical metric

The scorer must be implemented and pass synthetic coordinate, topology,
eligibility, weighted-median, and bootstrap tests before the label asset or
fixed-box chunks are opened.

- OBJ coordinates are z/y/x level-0 index coordinates. Convert them to the
  physical-label local L1 coordinates by
  `p_label = (p_L0 - 0.5) / 2 - [3936,0,0]`. This is the voxel-center mapping
  consistent with the reference evaluator's pairwise L0-to-L1 max pooling.
- Read only the fixed label-local window `z [976,1072), y [1520,1872),
  x [1520,1872)`. It provides at least 16 L1 voxels around the nominal mesh
  footprint.
- For every label z plane, compute the two-dimensional Euclidean distance
  transform to bit 8 (`recto_band`). Sample the resulting distance stack
  trilinearly at each mesh vertex. Using per-plane 2-D distances matches the
  reference's axial side-of-sheet instrument and prevents a nearer surface in
  another z plane from winning through a 3-D shortcut.
- Eligibility is frozen from the pre-snap mesh and common metadata only: the
  vertex is face-used; all coordinates and fixed pre-snap face-area weights
  are finite; every arm stays inside the fixed L0 box with a one-voxel inner
  margin; all eight trilinear label corners are valid for pre and both arms;
  and the pre-snap vertex is within 3.0 L1 voxels (56.2 um) of a physical
  recto band. Arm outcomes may not change the eligible set.
- Give each vertex one third of the area of each incident pre-snap triangle.
  For each provenance cube and arm, the score is the fixed-area-weighted
  median distance to the recto band, in L1 voxels. The paired cube effect is
  `distance(iter4) - distance(iter0)`, so positive favors the candidate.
- The cube is the inferential unit. Use seed `20260810` and 10,000 paired
  bootstrap resamples of the cube effects. Also report every cube, the
  unweighted sensitivity result, vertex/area coverage, and all y-row and
  x-column medians.

## Frozen gates and decision rule

Input and safety gates:

- exactly 25 nonzero uint8 RAW cubes and 25 prediction cube files, with no
  grid holes;
- at least 20 cubes placed, with zero failed, skipped, or low-confidence
  accepted inputs; audit turn-off pairs <=5% and `|du|<2` completeness >=75%;
- at least 20 cubes score, each with at least 100 common eligible face-used
  vertices, and at least three scored cubes in every y row and x column;
- identical vertex count, face indices, UVs, cube ranges and finite
  coordinates across pre and both arms;
- explicit iteration 0 and the no-flag default are byte-identical; and
- two scorer invocations emit byte-identical JSON.

The candidate is **physically superior** only if every input/safety gate
passes and all of these hold:

1. median paired cube effect >=0.10 L1 voxel (1.87 um);
2. paired bootstrap 95% interval lower bound >0;
3. at least 80% of scored cubes have positive effect;
4. every y-row and x-column median is positive; and
5. no more than 10% of cubes favor iter4 by more than 0.25 L1 voxel.

The former four-iteration behavior is **physically superior** only if the
same five gates pass after reversing the sign. Otherwise the physical result
is a tie/inconclusive default decision. No threshold, box, eligibility rule,
or arm may be changed after access.

Run the published sigma-1, fixed-pre-normal CT-ridge scorer as a secondary
metric with the same 25 cubes. It cannot override the physical decision.

## Frozen visual and publication rule

The only qualitative plane is global L0 z=9920, fixed as the midpoint of the
slab. Assemble all 25 RAW tiles and overlay the physical recto band plus the
pre, iter0, and iter4 mesh intersections. Also render the complete fixed 5x5
per-cube effect heatmap. No slice or crop search is allowed.

Publish the complete positive, negative, or inconclusive result with credit
to the physical-audit/#1382 authors, exact commands, hashes, per-cube data,
and limitations. If the physical result favors iter4, amend the PR #11
default claim rather than hiding the contradiction. If an input/safety gate
fails, publish only a clearly labeled feasibility failure; it supports no
physical arm claim.
