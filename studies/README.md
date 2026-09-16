# studies — the parametric campaign

This folder is the campaign itself: what is to be run, what came back, and the
report each stage produced. One folder per study, and inside it

| file | what it is | written by |
|---|---|---|
| `study.json` | the definition — the base case and every case in the matrix | committed by hand, or by the 파라메트릭 tab |
| `results.json` | one record per case as it finishes: mesh, Δp, convergence | the app, as the study runs |
| `report.html` | the stage report, standalone | the app, or `python3 -m …` |

A study definition is small, exact and reviewable; the results are the
measurements. Both are worth committing, and both are worth reading before
the next stage is planned. The meshes and the Fluent case/data files are not
here — those live in `runs/` and are ignored by git, because every one of them
is reproducible from `study.json` plus the code revision recorded in
`results.json`.

## The campaign

| stage | folder | what it settles |
|---|---|---|
| 1 | `stage1-correlations/` | which correlations exist, what each is allowed to be asked, and how far apart they already are |
| 2 | `rod-inline-mesh/`, `rod-staggered-mesh/` | how fine the mesh has to be before Δp stops moving — the mesh every later case uses |
| 3 | `rod-inline-sweep/`, `rod-staggered-sweep/` | Δp over the arrangement/condition matrix, against the correlations |
| 4 | *(later)* | the helical coil, and a correlation of our own |

Stage 1 needs no CFD and is done: `stage1-correlations/report.html`.

Stages 2 and 3 need Fluent. **They must be run on a machine with a licence.**
Anything produced by the mock backend is labelled MOCK and is plumbing
evidence, never a result — a MOCK number must not be quoted, plotted against a
correlation, or carried into a fit.
