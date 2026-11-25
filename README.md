# lrmap

Multi-scale spatial co-location of ligands and receptors between cell types.

For every sender cell type, receiver cell type and ligand–receptor pair, `lrmap()` measures how much ligand from the sender type surrounds receiver cells that express the receptor. It does this at three spatial scales:
- **contact:** the 8 nearest cells;
- **paracrine:** a Gaussian kernel, σ 25 µm;
- **niche:** a Gaussian kernel, σ 150 µm.

All three ranges can be set by the user (`k`, `sigma_p`, `sigma_n`).

## Tests
Significance comes from toroidal shifts of the sender layer against the receiver layer. Each layer keeps its own spatial pattern; only their alignment is broken. Every spatial row carries three tests:

| Columns | Question | Use |
|---|---|---|
| `q_type`, `fold_type` (default) | Do ligand-high senders sit nearer receptor-high receivers than the two cell types do on average? | Pair-specific signal, calibrated by the shifts |
| `q`, `fold` | Do ligand-expressing senders and receptor-expressing receivers co-locate? | Includes the co-location of the two cell types, so any pair expressed in two neighbouring types is called |
| `q_shuffle` | As `q_type`, against ligand and receptor shuffled within their types (layout fixed) | Ignores smooth spatial trends; anti-conservative when they are present |

`lrmap_results()` returns the calls (`test = "type"` by default). The type-level co-location of every sender–receiver pair is in `fit$coloc`.

On a whole breast section (169,528 cells, 15 types, 1,387 CellPhoneDB pairs), the default test called 1.3–2.2% of ligand–receptor tests per scale, against 0.1–0.3% of random gene pairs. The co-location-only test called random pairs at 2.6–4.7%.

## Leak
Transcripts spill between touching cells. Uncorrected, a receptor leaked from a sender into its neighbours looks like an interface signal. `lrmap_leak()` estimates the leak per receiver–sender type pair across many genes:
- λ is the fraction of the sender's expression picked up by the receiver;
- μ is the change in the receiver's own baseline.

Pairs with little contact are pooled toward the median of their receiver type. `lrmap(leak = )` then removes the leak from every gene before forming ligand and receptor values. The model is the one used by [edgeicg](https://github.com/drieslab/edgeicg); lrmap carries its own copy, and an edgeicg fit can be passed as `leak`.

## Expression only
`scales = "expression"` adds the single-cell style test: mean ligand in the sender type × mean receptor in the receiver type, against shuffled cell labels, on the same tests and values. `lrmap_compare()` then classifies each pair as expression only, spatial only, both or neither. On the breast section, 96% of expression-only calls had no spatial support beyond the cell types, and 72% were between cell types that do not co-locate at any scale.

## Use
```r
library(lrmap)
s <- lrmap_simulate("touch", induce = 1, leak = 0.05)
lk <- lrmap_leak(s$counts, s$coords, s$types, lib = s$lib, min_cells = 30)
fit <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, leak = lk, min_shift = 300,
             scales = c("contact", "paracrine", "niche", "expression"))
fit
lrmap_results(fit)
```

- **Input:** `counts` can be a genes × cells matrix or a reader `function(genes)`, for example over a GiottoDisk store. Only the genes of the ligand–receptor table are read; `lib` (library size over all genes) is then required.
- **Ligand–receptor table:** columns `ligand` and `receptor`, and optionally `ligand_subunits` and `receptor_subunits` with genes joined by `_`. A complex ligand is the minimum over its subunits, and a complex receptor the mean.

## Speed
With [lrmaprs](https://github.com/drieslab/lrmaprs) installed, `backend = "auto"` runs the per-shift loops in Rust.

On the breast section:
- **CellPhoneDB pairs:** 4.3 min on two threads (200 shifts, three scales, about 13,000 tests per scale).
- **Leak estimate on all 17,423 genes:** 7.8 min.
- **Peak memory:** 3.1 GB.

## Giotto
lrmap succeeds Giotto's original cell–cell communication test (`spatCellCellcom()`), which was contact-only.
