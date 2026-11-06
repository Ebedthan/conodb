# Classifier fit: UniProt 2026_04

- Training: 5604 sequences (case weights 1/group size), 18 families + NONE; features per family: SIG:hmm, SIG:pssm, SIG:pwm, PRO:hmm, PRO:pssm, PRO:pwm, MAT:hmm, MAT:pssm, MAT:pwm, FW:framework.
- baseline: per-profile HMM cutoffs from training negatives, for a 1% false-positive rate per sequence (Bonferroni over profiles); highest summed excess wins.
- map: penalised conditional logit (shrinkage: posterior means of tau); converged: True (154 iterations).
- HMM floor: 0 bits; regions: SIG, PRO, MAT; feature types: hmm, pssm, pwm, framework.
- shipped model (`fit/model_ship.json`): the penalised fit refitted on all 6897 labelled sequences (train + test), same shrinkage; converged: True. Performance reported here is from the train-only fit on held-out data.
- bayes: hierarchical model in Stan, 4000 draws in 4523 s (4 chains at a time x 2 threads); max R-hat 1.007, min bulk ESS 907, divergences 0.

Shared weights (posterior mean, standardised scores) and between-family spread:

| feature | beta | tau |
|---|---|---|
| SIG:hmm | 2.75 | 0.31 |
| SIG:pssm | -0.44 | 0.72 |
| SIG:pwm | -0.51 | 0.44 |
| PRO:hmm | 1.18 | 0.35 |
| PRO:pssm | 0.66 | 0.66 |
| PRO:pwm | -1.16 | 0.76 |
| MAT:hmm | 0.83 | 0.37 |
| MAT:pssm | 0.2 | 0.91 |
| MAT:pwm | -0.18 | 0.27 |
| FW:framework | 0.64 | 0.29 |


## Test set

| method | n | accuracy | balanced accuracy | macro F1 | positives missed (NONE) | positives wrong family | negatives called a family | by negative source |
|---|---|---|---|---|---|---|---|---|
| baseline | 1293 | 0.966 | 0.895 | 0.911 | 13.6% | 0.0% | 1.7% | shuffled 6.1%, uniprot 1.1% |
| map | 1293 | 0.997 | 0.972 | 0.982 | 0.5% | 1.6% | 0.0% | shuffled 0.0%, uniprot 0.0% |
| bayes | 1293 | 0.997 | 0.972 | 0.982 | 0.5% | 1.6% | 0.0% | shuffled 0.0%, uniprot 0.0% |

Recall per class:

| class | baseline | map | bayes |
|---|---|---|---|
| A | 0.94 | 1.00 | 1.00 |
| I1 | 1.00 | 1.00 | 1.00 |
| I2 | 1.00 | 1.00 | 1.00 |
| M | 0.77 | 1.00 | 1.00 |
| NONE | 0.98 | 1.00 | 1.00 |
| O1 | 0.96 | 1.00 | 1.00 |
| O2 | 0.73 | 0.91 | 0.91 |
| S | 1.00 | 1.00 | 1.00 |
| T | 0.68 | 0.84 | 0.84 |


Detailed evaluation (by sequence type, identity to training, calibration, ablations): `evaluation.md`.
