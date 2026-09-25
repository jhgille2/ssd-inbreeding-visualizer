#!/usr/bin/env Rscript
# Simulate single-seed descent (SSD) with AlphaSimR and export chromosome mosaics.
#
# Model:
#   * Two fully homozygous inbred parents, fixed for contrasting alleles at every
#     segregating site on a single 1 Morgan chromosome (allele 0 = Parent A,
#     allele 1 = Parent B, so allele identity is parental origin at every site).
#   * Parent A x Parent B -> F1 (100% heterozygous).
#   * SSD: the F1 is selfed; one random progeny (one seed) is selfed; repeated to
#     F8. AlphaSimR simulates meiosis with its default crossover model.
#   * Crossover positions are inferred as the midpoints between flanking markers
#     at each switch in parental origin along a transmitted gamete.
#
# Output: docs/ssd_data.json — the same schema produced by
# scripts/generate_data.py (a numpy port of this model), consumed by the
# interactive visualization in docs/.
#
# Requirements: R (>= 4.0), AlphaSimR (>= 1.0), jsonlite.
# Run from the repository root:  Rscript R/simulate_ssd.R
#
# Reference: Gaynor et al. (2021) AlphaSimR: an R package for breeding program
# simulations. G3 11(2): jkaa017. doi:10.1093/g3journal/jkaa017

suppressPackageStartupMessages({
  library(AlphaSimR)
  library(jsonlite)
})

SEED <- 174L
N_SEG <- 1000L
CHR_LEN <- 1.0   # Morgans
SELFING_GENS <- 7L  # F2..F8 after the F1

set.seed(SEED)

# --- Founders: two inbred parents, contrasting alleles at every site ----------
genMap <- list(seq(0, CHR_LEN, length.out = N_SEG))
hap <- matrix(c(rep(0L, N_SEG), rep(1L, N_SEG)), nrow = 2, byrow = TRUE)
founderPop <- newMapPop(genMap = genMap, haplotypes = list(hap), inbred = TRUE)

SP <- SimParam$new(founderPop)
founders <- newPop(founderPop, simParam = SP)  # ind 1: all-0, ind 2: all-1

# --- Helpers -----------------------------------------------------------------
pos_cm <- round(seq(0, 100, length.out = N_SEG), 4)

# Parental-origin vectors of one individual's two homologs (allele = origin here).
get_haplo_pair <- function(pop) {
  h1 <- pullSegSiteHaplo(pop, haplo = 1, chr = 1, simParam = SP)
  h2 <- pullSegSiteHaplo(pop, haplo = 2, chr = 1, simParam = SP)
  list(h1 = as.integer(h1[1L, ]), h2 = as.integer(h2[1L, ]))
}

# Crossover positions (cM): midpoints between flanking markers at origin switches.
crossovers_cm <- function(h) {
  sw <- which(diff(h) != 0L)
  if (!length(sw)) return(numeric(0))
  round((pos_cm[sw] + pos_cm[sw + 1L]) / 2, 3)
}

generation_record <- function(label, selfing_gens, hp) {
  list(
    label = label,
    selfing_generations = selfing_gens,
    expected_heterozygosity = 0.5^selfing_gens,  # F1: 1.0, F2: 0.5, ...
    homologs = list(
      list(origin = hp$h1, crossovers_cm = crossovers_cm(hp$h1)),
      list(origin = hp$h2, crossovers_cm = crossovers_cm(hp$h2))
    )
  )
}

# --- Pedigree ----------------------------------------------------------------
F1 <- makeCross(founders,
  crossPlan = matrix(c(1L, 2L), nrow = 1),  # female = Parent A, male = Parent B
  simParam = SP
)
generations <- list(generation_record("F1", 0L, get_haplo_pair(F1)))

pop <- F1
for (i in seq_len(SELFING_GENS)) {
  pop <- self(pop, nProgeny = 1, keepParents = FALSE, simParam = SP)
  generations[[i + 1L]] <- generation_record(
    paste0("F", i + 1L), i, get_haplo_pair(pop)
  )
}

# --- Export ------------------------------------------------------------------
data <- list(
  meta = list(
    title = "Single-seed descent: from cross to inbred line",
    organism = paste(
      "Glycine max (soybean); one representative chromosome is shown",
      "(soybean has 20 chromosome pairs, 2n = 40)"
    ),
    chromosome_length_cm = 100.0,
    n_markers = N_SEG,
    parents = list("0" = "Parent A", "1" = "Parent B"),
    model = paste(
      "Single-seed descent simulated with AlphaSimR: F1 selfed; one random",
      "progeny selfed each generation to F8. Meiosis uses AlphaSimR's default",
      "crossover model; one random gamete transmitted per meiosis (single seed)."
    ),
    generator = "R/simulate_ssd.R (AlphaSimR)",
    seed = SEED
  ),
  marker_pos_cm = pos_cm,
  generations = generations
)

# Resolve docs/ssd_data.json relative to this script's location (repo root/R).
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
out_path <- normalizePath(file.path(script_dir, "..", "docs", "ssd_data.json"),
  mustWork = FALSE
)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write_json(data, out_path, auto_unbox = TRUE, digits = NA)

# Sanity report: observed vs expected homozygosity.
cat("Wrote", out_path, "\n")
cat(sprintf("%4s %7s %7s %8s %8s\n", "gen", "n_xo_h1", "n_xo_h2", "obs_hom", "exp_hom"))
for (g in generations) {
  h1 <- g$homologs[[1]]$origin
  h2 <- g$homologs[[2]]$origin
  cat(sprintf(
    "%4s %7d %7d %8.3f %8.3f\n", g$label,
    length(g$homologs[[1]]$crossovers_cm),
    length(g$homologs[[2]]$crossovers_cm),
    mean(h1 == h2), 1 - g$expected_heterozygosity
  ))
}
