# studies — the parametric campaign

This folder is the campaign itself: what is to be run, what came back, and the
report each stage produced. One folder per study, and inside it

| file | what it is | written by |
|---|---|---|
| `study.json` | the definition — the base case and every case in the matrix | committed by hand, or by the 파라메트릭 tab |
| `results.json` | one record per case as it finishes: mesh, Δp, convergence | the app, as the study runs |
| `report.html` | the stage report, standalone | the app, or `python3 -m …` |  ← not in git
| `overview.html` | every case on one screen, standalone | the 파라메트릭 tab's **한눈에 보기** button |  ← not in git
| `history/*.csv` | one file per case: the residual and Δp traces as plain text | the 파라메트릭 tab's **CSV로 내보내기** button |

**`study.json` and `results.json` are in git; the HTML is not.** A report is
rebuilt by a button in seconds — the stage-1 one from `correlations.py` alone,
the others from the `results.json` beside them — so committing it buys nothing
and costs something real: the first person to press that button has a dirty
tree and a `git pull` that aborts on a file they never meant to edit. The
definition and the measurement are the things worth keeping.

A study definition is small, exact and reviewable; the results are the
measurements.  `history/` is `results.json` in a second spelling, for plotting
the traces somewhere that is not this app — it holds nothing the record does
not, so it is regenerable at any time and not worth committing. Both are worth committing, and both are worth reading before
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
