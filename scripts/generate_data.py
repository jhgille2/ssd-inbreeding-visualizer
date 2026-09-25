#!/usr/bin/env python3
"""Simulate single-seed descent (SSD) for one chromosome and export the result.

Model (mirrors R/simulate_ssd.R, the AlphaSimR reference implementation):
  * Two fully homozygous inbred parents, fixed for contrasting alleles at every
    marker on a single 100 cM chromosome (allele 0 = Parent A, allele 1 = Parent B,
    so allele identity is parental origin at every marker).
  * Parent A x Parent B -> F1 (100% heterozygous).
  * SSD: the F1 is selfed; one random progeny (one seed) is selfed; repeated to F8.
    Each selfing involves two independent meioses in the parent plant; one random
    chromatid is transmitted per gamete (single seed).
  * Meiosis follows the Haldane model: crossover count ~ Poisson(chromosome length
    in Morgans), crossover positions uniform along the chromosome, no interference.

Output: docs/ssd_data.json, consumed by the interactive visualization in docs/.
The same JSON schema is produced by R/simulate_ssd.R.
"""

import json
from pathlib import Path

import numpy as np

SEED = 174
N_MARKERS = 1000
CHR_LEN_CM = 100.0
SELFING_GENS = 7  # F2..F8 after the F1
OUT_PATH = Path(__file__).resolve().parent.parent / "docs" / "ssd_data.json"


def meiosis(hom1, hom2, marker_pos, rng):
    """Simulate one meiosis in a plant with homologs hom1/hom2.

    Returns (gamete_origin, crossover_positions_cm): the parental-origin array of
    one randomly chosen chromatid and the crossover positions (cM) of that meiosis.
    """
    mean_xo = (marker_pos[-1] - marker_pos[0]) / 100.0
    n_xo = int(rng.poisson(mean_xo))
    if n_xo:
        xo = np.sort(rng.uniform(marker_pos[0], marker_pos[-1], n_xo))
    else:
        xo = np.zeros(0)
    start = int(rng.integers(0, 2))  # which homolog the chromatid starts on
    seg = np.searchsorted(xo, marker_pos, side="right")
    use_h1 = (seg + start) % 2 == 0
    gamete = np.where(use_h1, hom1, hom2).astype(int)
    return gamete, np.round(xo, 3)


def generation_record(label, selfing_gens, g1, xo1, g2, xo2):
    return {
        "label": label,
        "selfing_generations": selfing_gens,
        "expected_heterozygosity": 0.5**selfing_gens,  # F1: 1.0, F2: 0.5, ...
        "homologs": [
            {"origin": [int(x) for x in g1], "crossovers_cm": [float(x) for x in xo1]},
            {"origin": [int(x) for x in g2], "crossovers_cm": [float(x) for x in xo2]},
        ],
    }


def main():
    rng = np.random.default_rng(SEED)
    marker_pos = np.linspace(0.0, CHR_LEN_CM, N_MARKERS)

    parent_a = np.zeros(N_MARKERS, dtype=int)
    parent_b = np.ones(N_MARKERS, dtype=int)

    # F1: intact parental chromosomes; no recombination has occurred yet.
    g1, xo1 = parent_a.copy(), np.zeros(0)
    g2, xo2 = parent_b.copy(), np.zeros(0)
    generations = [generation_record("F1", 0, g1, xo1, g2, xo2)]

    plant = (g1, g2)
    for i in range(1, SELFING_GENS + 1):
        g1, xo1 = meiosis(plant[0], plant[1], marker_pos, rng)
        g2, xo2 = meiosis(plant[0], plant[1], marker_pos, rng)
        plant = (g1, g2)
        generations.append(
            generation_record(f"F{i + 1}", i, g1, xo1, g2, xo2)
        )

    data = {
        "meta": {
            "title": "Single-seed descent: from cross to inbred line",
            "organism": "Glycine max (soybean); one representative chromosome is "
            "shown (soybean has 20 chromosome pairs, 2n = 40)",
            "chromosome_length_cm": CHR_LEN_CM,
            "n_markers": N_MARKERS,
            "parents": {"0": "Parent A", "1": "Parent B"},
            "model": "Single-seed descent: F1 selfed; one random progeny selfed "
            "each generation to F8. Meiosis: Haldane model (Poisson crossover "
            "count with mean = chromosome length in Morgans, uniform crossover "
            "positions, no interference); one random chromatid transmitted per "
            "gamete.",
            "generator": "scripts/generate_data.py (numpy port of the AlphaSimR "
            "model in R/simulate_ssd.R)",
            "seed": SEED,
        },
        "marker_pos_cm": [round(float(x), 4) for x in marker_pos],
        "generations": generations,
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(data))

    # Sanity report: observed vs expected homozygosity.
    print(f"Wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.1f} KiB)")
    print(f"{'gen':>4} {'n_xo_h1':>7} {'n_xo_h2':>7} {'obs_hom':>8} {'exp_hom':>8}")
    for g in generations:
        h1 = np.array(g["homologs"][0]["origin"])
        h2 = np.array(g["homologs"][1]["origin"])
        obs_hom = float(np.mean(h1 == h2))
        exp_hom = 1.0 - g["expected_heterozygosity"]
        print(
            f"{g['label']:>4} {len(g['homologs'][0]['crossovers_cm']):>7} "
            f"{len(g['homologs'][1]['crossovers_cm']):>7} {obs_hom:>8.3f} {exp_hom:>8.3f}"
        )


if __name__ == "__main__":
    main()
