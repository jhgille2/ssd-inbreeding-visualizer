# Single-Seed Descent: From Cross to Inbred Line

An interactive visualization of how single-seed descent (SSD) turns a segregating
biparental cross into inbred lines. One chromosome of one lineage is tracked from
the F1 to the F8: each generation shows the plant's two homologs as a mosaic of the
two parental genomes, the crossover positions from the meioses that formed them, and
the resulting rise in homozygosity.

**Live site:** https://jhgille2.github.io/ssd-inbreeding-visualizer/

## Background

In a typical soybean breeding program, two inbred parents are crossed and the
segregating progeny are advanced by single-seed descent: each generation, one seed
per plant is grown and selfed, with no selection, until the lines are nearly
homozygous. Two processes run in parallel:

1. **Meiosis with recombination** shuffles the parental genomes, so each gamete
   carries a mosaic of Parent A and Parent B segments.
2. **Repeated selfing** halves the remaining heterozygosity each generation
   (F2 ~50% homozygous, F3 ~75%, F4 ~87.5%, ...), fixing those mosaics into
   inbred lines.

Recombination creates the mosaic; inbreeding fixes it. This repository simulates
both and visualizes them generation by generation.

## Repository layout

| Path | Contents |
|---|---|
| `docs/` | The interactive site (GitHub Pages). `index.html`, `app.js`, `style.css`, and the dataset `ssd_data.json`. No build step; no dependencies. |
| `R/simulate_ssd.R` | Reference simulation in R using [AlphaSimR](https://github.com/gaynorr/AlphaSimR). Writes `docs/ssd_data.json`. |
| `scripts/generate_data.py` | Numpy port of the same model, used to generate the committed dataset. |
| `.github/workflows/pages.yml` | Deploys `docs/` to GitHub Pages. |
| `.github/workflows/regenerate-data.yml` | Re-runs the R/AlphaSimR simulation (manual or monthly) and commits a fresh dataset. |

## The simulation model

- One representative 100 cM chromosome with 1,000 evenly spaced markers. (Soybean
  has 20 chromosome pairs, 2n = 40; every pair follows the same process.)
- Two fully homozygous parents, fixed for contrasting alleles at every marker, so
  allele identity equals parental origin at every position.
- Parent A x Parent B -> F1; the F1 is selfed; one random progeny is selfed each
  generation through F8. Each selfing involves two independent meioses in the
  parent plant; one random chromatid is transmitted per gamete (the single seed).
- Meiosis follows the Haldane model: crossover counts are Poisson with mean equal
  to the chromosome length in Morgans, crossover positions are uniform, and there
  is no crossover interference. No selection, no mutation, no segregation distortion.

The committed dataset is one realized lineage (seed 174), chosen because its
homozygosity trajectory tracks the expectation closely. Any single line deviates
by chance; the "Simulate a new line" button on the site re-runs the same model in
the browser with a fresh seed so that sampling variance can be explored.

## Regenerating the dataset

R (reference implementation; requires AlphaSimR and jsonlite):

```r
# from the repository root
Rscript R/simulate_ssd.R
```

Python (numpy port; used for the committed dataset):

```bash
python3 scripts/generate_data.py
```

Both write `docs/ssd_data.json` in the same schema. To view the site locally,
serve `docs/` over HTTP, e.g. `python3 -m http.server --directory docs`, and open
http://localhost:8000.

## Citation

Simulation engine: Gaynor et al. (2021) AlphaSimR: an R package for breeding
program simulations. *G3* 11(2): jkaa017. doi:10.1093/g3journal/jkaa017.
