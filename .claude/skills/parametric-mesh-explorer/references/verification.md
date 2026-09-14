# Verification

The reason this tool can be trusted is not that the code looks right — it is
that every structural claim has a check that could have failed and did not.
Reproduce this battery for any new geometry.

## The check battery

Run these at build time, not as a separate test suite, so a bad parameter set
surfaces immediately.

| Check | Expected | Catches |
|---|---|---|
| shared-node disagreement | exactly `0` | two blocks using one key for different points |
| 2-D open edges all on a real boundary | `0` cracks | two blocks using different keys for one point |
| meshed volume vs analytic expectation | ~1e-14 relative | node placement, connectivity, extrusion, volume formula — all at once |
| minimum cell area / volume | `> 0` | inverted or degenerate cells |
| faces per cell | exactly 6 | dropped or duplicated faces |
| patch face total | equals boundary face count | a patch silently losing faces |
| `owner < neighbour`, ascending order | all faces | polyMesh ordering rules |
| face normals away from owner | all faces | orientation errors |
| **cell volume rebuilt from the faces** | all `> 0`, total = analytic | inverted faces, degenerate cells, unit slips |
| orphan / duplicate nodes | `0` | welding gone wrong in the other direction |

The cell-volume row is the one that would have caught the Fluent orientation
bug, so make it a gate on export rather than a report line: integrate `x.n dA`
over each cell's faces exactly as a solver does, and refuse to write the file if
any volume is non-positive or the total misses the analytic volume.

### The volume check is the strong one

Compare against the volume of the **faceted** domain, not the ideal one. For a
polygonal inclusion through the arc nodes the removed area is exactly

```
Σ over sectors of ½·R²·n_az·sin(Δθ / n_az)
```

Accumulate it while building blocks. The meshed volume then matches to ~1e-14,
which only happens if node positions, connectivity, extrusion and the volume
integral are all simultaneously right. Comparing against `πR²` instead leaves a
0.04 % residual that hides real errors — report that separately as the
discretisation error it is (and it tells the user what their azimuthal count
costs them in flow area).

## Independent cross-checks

**Two implementations.** When both a Python and a browser build exist, compare
order-sensitive checksums of the node and quad tables:

```
h = Σ (i+1)·(round(x·1e6) + 7919·round(y·1e6))
```

Identical checksums mean identical tables in identical order — far stronger than
each passing its own tests. This confirmed the browser port reproduced the
Python generator exactly.

**External readers.** Round-trip the output through a tool you did not write:
`meshio` for VTU / Fluent / STL, Python's `zipfile` for the archive. Independent
parsers catch format errors that self-consistency never will.

**Self round-trip.** For the zip, inflate every entry back and compare bytes and
CRC against the source text. Cheap and conclusive.

## Test the parameter space, not the default

Two of three real bugs appeared only away from defaults. Always sweep at least:

- coarse resolution (lowest azimuthal / radial counts) — tolerance failures
- `nZ = 1` — quasi-2-D degenerate paths
- single row / single column — loops that assume neighbours exist
- an absent region (`L_in = 0`) — blocks that should not be built
- the toggles that change topology (half inclusions on/off, offset row on/off)
- an invalid configuration — must be refused, not silently exported

## Bugs actually hit, and what they teach

**Corner nodes named by side (2600 cracked faces).** A polygon corner belongs to
several sides at once. Priority-ordered keys, corners first. → any node on more
than one entity needs the higher-priority name.

**Rod-wall faces classified by centroid.** Passed at 36 azimuthal divisions
(0.4 % off), failed at 12 (3.4 % off) with a 2 % tolerance. → classify by
points that lie exactly on the surface, not by derived points that only
approximately do.

**Radial cut tools starting at the inclusion centre.** For a half inclusion the
centre sits on the domain wall, so the tool split the wall edge and produced
phantom segments. → start the tool strictly inside the inclusion (e.g. `0.5R`)
so it crosses the surface transversally.

**Preview and export drifting apart.** They were separate code paths at first.
→ one block list feeds both; the preview is then a genuine prediction, and when
the export says 49.1° max skew the preview said 49.1° too.

**Exactly-30° cells.** 94.8 % vs 94.3 % below 30° between two implementations
looked like a discrepancy; it was 1920 cells sitting exactly on 30.000° (the
hexagon corner angle) falling either side of `<` vs `≤`. → before chasing a
numeric difference, check for a spike at the threshold.

**All cells non-positive, in the solver.** Fluent rejected 610560 of 610560
cells. The share is the diagnosis: a geometry fault hits some cells, a global
orientation flip hits every one. Before hunting for degenerate blocks, check
whether the failure rate is exactly 100 %.

**A flag whose two branches did the same thing.** The Python twin had a
`FLUENT_NORMAL_C0_TO_C1` switch whose "off" branch reversed the node order *and*
swapped `c0`/`c1`. Those cancel on an internal face, so the flag did nothing
where it mattered and corrupted the owner slot on boundary faces. An escape
hatch nobody has exercised is not an escape hatch. → exercise both branches, or
delete the flag and pin the correct behaviour.

**The contour gated on the mesh-preview toggle.** Turning off "grid preview"
blanked the velocity and pressure fields, because the fill lived inside the same
`if` as the mesh lines. → separate "what am I drawing" (the mode) from "do I
also want lines on top" (the toggle); a coloured map is the requested picture,
not an overlay.

**Two switches that need a truth table.** Both bugs above were combinations, not
single states. When a feature is governed by more than one toggle, sweep every
combination and record a number for each — pixel counts work fine — rather than
checking the state you had open.

**Stale browser snapshots.** The preview pane renders local files as snapshots.
Edits made after the last navigation are invisible, so tests "fail" against code
that no longer exists — and worse, "pass" against code you already fixed.
→ re-navigate after every edit, always.

**`localStorage` on `data:` URLs.** Silently unavailable in the preview. Wrap
persistence in try/catch and do not conclude it is broken from the snapshot.

## Reporting

State what was checked with numbers, and state what was not checked and why. A
verified claim and a plausible one should never read the same.

The Fluent `c0`/`c1` convention used to be the standing example of an unverified
claim, and carrying that caveat honestly for months is what made the eventual
failure legible: when Fluent rejected every cell, the caveat said where to look
first. It is now settled (see `mesh-core.md`) — but note *how* it was settled.
It was not settled by reading a specification; it was settled by a solver run.
Keep that distinction in reports. What a local check can prove is that the mesh
is self-consistent under a stated convention; only the target solver can confirm
the convention.

When something is confirmed externally, go back and delete the old caveat from
these notes. A stale "unverified" warning is as misleading as a missing one.
