# poweRanger <img src="media/hexlogo.png" align="right" height="200"/>

> *In the deep forest live the echoes from the past*

[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)](https://www.r-project.org/)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange)](https://lifecycle.r-lib.org/articles/stages.html)
[![platform](https://img.shields.io/badge/platform-linux-lightgrey?logo=linux)](https://www.linux.org/)

---

## What is Gradient Forest?

Genotype-environment association (GEA) methods ask: *which parts of the genome vary with the environment, and why?* Most approaches — LFMM, BayEnv, linear RDA — answer this with a test per locus: is this SNP significantly associated with temperature? precipitation? They are powerful when the signal is strong, the associations are linear, and population structure is well-characterised.

Gradient Forest asks a different question: **where along an environmental gradient does the genome change most?**

Instead of a p-value per locus, it produces a cumulative turnover function F(x) per predictor per site across individuals — a curve that shows at which values of, say, mean annual temperature, allele frequencies shift the fastest across the landscape. The area under that curve, weighted by the out-of-bag R² of each locus model, becomes a biologically meaningful unit of compositional change that is comparable across predictors and across loci.

This makes Gradient Forest particularly strong when:

- **Population structure is weak, absent or simply unknown** — no latent factor correction needed, the signal lives in the split density of the trees
- **Relationships are non-linear or threshold-driven** — random forests capture any shape of response without transformation or model selection
- **Many loci contribute small effects** — the ensemble aggregation rewards polygenic signal that univariate tests miss
- **You need to know *where*, not just *whether*** — the turnover curve localises the genomic response along the gradient, which is essential for genetic offset and assisted gene flow predictions

Where linear RDA wins on interpretability and speed, and LFMM wins on controlling for confounded structure, Gradient Forest wins on flexibility, scalability, and the richness of its output.

---

## What is poweRanger?

**poweRanger** implements the Gradient Forest variance partitioning algorithm as a post-hoc analysis suite for [`ranger`](https://github.com/imbs-hl/ranger) random forest models. It operates entirely on models you have already trained — one per locus, one per genetic PC, or any other univariate response — and reconstructs the turnover curves, importance rankings, and diagnostic summaries that constitute a full Gradient Forest analysis.

It does not wrap the original `gradientForest`/`extendedForest` packages. It reimplements the core algorithm from Ellis et al. (2012) on top of `ranger`, which is faster, actively maintained, and scales to whole-genome datasets.

Because poweRanger receives a plain named list of `ranger` objects, it is also compatible with **RDA-Forest** (Matz & Black 2025) out of the box: train one ranger per genetic PC instead of per SNP, pass the list in, and the turnover curves describe polygenic compositional turnover rather than single-locus turnover.

```r
models |>
  filter_models(min_r2 = 0.01) |>
  turnover_curves(env, vars = c("bio1", "bio12")) -> gf

gf |> gf_summarise() -> sm

plot_curves(gf)
plot_diagnostic(sm)
plot_importance(sm)
```

---

## Installation

```r
remotes::install_github("jcarocont/poweRanger")
```

**Dependencies:** `ranger`, `dplyr`, `tibble`, `ggplot2`, `purrr`

---

## Citation

If you use poweRanger, please cite the original Gradient Forest paper:

> Ellis N, Smith SJ, Pitcher CR (2012). Gradient forests: calculating importance gradients on physical predictors. *Ecology* 93(1): 156–168. https://doi.org/10.1890/11-0252.1

For the RDA-Forest workflow:

> Matz MV, Black KL (2025). RDAforest: Identifying environmental drivers of polygenic adaptation. *Molecular Ecology Resources*. https://doi.org/10.1111/1755-0998.70002
