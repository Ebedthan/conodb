# Evaluation: UniProt 2026_04

Settings: feature types hmm, pssm, pwm, framework; regions SIG, PRO, MAT; HMM floor 0 bits. Region prediction: `anatomy.md`.

All numbers are on held-out sequences (test set, and temporal test set if enabled); no model saw them. `positives correct` = right family; `family if known` = right family when NONE is excluded (a sequence known to be a conopeptide).

## Test set

### Overall

| method | subset | n | accuracy | positives correct | family if known | called NONE | negatives called a family | Brier | log loss |
|---|---|---|---|---|---|---|---|---|---|
| baseline | all | 1293 | 0.966 | 86.4% | - | 13.6% | 1.7% | - | - |
| map | all | 1293 | 0.997 | 97.8% | 97.8% | 0.5% | 0.0% | 0.005 | 0.013 |
| bayes | all | 1293 | 0.997 | 97.8% | 97.8% | 0.5% | 0.0% | 0.005 | 0.015 |

### By sequence type

| method | subset | n | accuracy | positives correct | family if known | called NONE | negatives called a family | Brier | log loss |
|---|---|---|---|---|---|---|---|---|---|
| baseline | complete precursor | 130 | 1.000 | 100.0% | - | 0.0% | - | - | - |
| baseline | negative: shuffled | 131 | 0.939 | - | - | - | 6.1% | - | - |
| baseline | negative: uniprot | 978 | 0.989 | - | - | - | 1.1% | - | - |
| baseline | no signal peptide | 50 | 0.500 | 50.0% | - | 50.0% | - | - | - |
| baseline | signal, incomplete | 4 | 1.000 | 100.0% | - | 0.0% | - | - | - |
| map | complete precursor | 130 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.000 | 0.000 |
| map | negative: shuffled | 131 | 1.000 | - | - | - | 0.0% | 0.001 | 0.009 |
| map | negative: uniprot | 978 | 1.000 | - | - | - | 0.0% | 0.000 | 0.004 |
| map | no signal peptide | 50 | 0.920 | 92.0% | 92.0% | 2.0% | - | 0.120 | 0.250 |
| map | signal, incomplete | 4 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.000 | 0.000 |
| bayes | complete precursor | 130 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.000 | 0.001 |
| bayes | negative: shuffled | 131 | 1.000 | - | - | - | 0.0% | 0.002 | 0.012 |
| bayes | negative: uniprot | 978 | 1.000 | - | - | - | 0.0% | 0.001 | 0.004 |
| bayes | no signal peptide | 50 | 0.920 | 92.0% | 92.0% | 2.0% | - | 0.119 | 0.267 |
| bayes | signal, incomplete | 4 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.000 | 0.002 |

### Positives by highest identity to a training conopeptide

| method | subset | n | accuracy | positives correct | family if known | called NONE | negatives called a family | Brier | log loss |
|---|---|---|---|---|---|---|---|---|---|
| baseline | [0, 0.5) | 1 | 1.000 | 100.0% | - | 0.0% | - | - | - |
| baseline | [0.5, 0.6) | 18 | 0.889 | 88.9% | - | 11.1% | - | - | - |
| baseline | [0.6, 0.7) | 47 | 0.894 | 89.4% | - | 10.6% | - | - | - |
| baseline | [0.7, 0.8) | 79 | 0.823 | 82.3% | - | 17.7% | - | - | - |
| baseline | [0.8, 0.9) | 36 | 0.917 | 91.7% | - | 8.3% | - | - | - |
| baseline | [0.9, 1] | 3 | 0.667 | 66.7% | - | 33.3% | - | - | - |
| map | [0, 0.5) | 1 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.000 | 0.003 |
| map | [0.5, 0.6) | 18 | 0.944 | 94.4% | 94.4% | 5.6% | - | 0.073 | 0.136 |
| map | [0.6, 0.7) | 47 | 0.979 | 97.9% | 97.9% | 0.0% | - | 0.038 | 0.064 |
| map | [0.7, 0.8) | 79 | 0.975 | 97.5% | 97.5% | 0.0% | - | 0.036 | 0.082 |
| map | [0.8, 0.9) | 36 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.002 | 0.013 |
| map | [0.9, 1] | 3 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.002 | 0.029 |
| bayes | [0, 0.5) | 1 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.000 | 0.015 |
| bayes | [0.5, 0.6) | 18 | 0.944 | 94.4% | 94.4% | 5.6% | - | 0.071 | 0.140 |
| bayes | [0.6, 0.7) | 47 | 0.979 | 97.9% | 97.9% | 0.0% | - | 0.029 | 0.053 |
| bayes | [0.7, 0.8) | 79 | 0.975 | 97.5% | 97.5% | 0.0% | - | 0.039 | 0.095 |
| bayes | [0.8, 0.9) | 36 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.006 | 0.022 |
| bayes | [0.9, 1] | 3 | 1.000 | 100.0% | 100.0% | 0.0% | - | 0.003 | 0.037 |

### Calibration

Within each confidence bin, accuracy should match the mean confidence.

| method | bin | n | mean_confidence | accuracy |
|---|---|---|---|---|
| map | [0.5, 0.7) | 6 | 0.612 | 0.5 |
| map | [0.7, 0.9) | 12 | 0.839 | 0.917 |
| map | [0.9, 0.99) | 90 | 0.974 | 1.0 |
| map | [0.99, 1) | 1185 | 0.998 | 1.0 |
| bayes | [0, 0.5) | 1 | 0.494 | 0.0 |
| bayes | [0.5, 0.7) | 7 | 0.617 | 0.571 |
| bayes | [0.7, 0.9) | 16 | 0.816 | 1.0 |
| bayes | [0.9, 0.99) | 101 | 0.972 | 1.0 |
| bayes | [0.99, 1) | 1168 | 0.998 | 1.0 |

### Variant grid on the test set (penalised fits)

Reported after the defaults were chosen by cross-validation (`selection.md`); **bold** = the recommended variant. Intervals: bootstrap over entry groups. `features` per variant in `evaluation_ablations.tsv`.

#### all

| variant | note | n | positives correct [95% CI] | family if known | negatives called a family | log loss [95% CI] |
|---|---|---|---|---|---|---|
| censored · HMM |  | 1293 | 87.0% [78.5%, 94.4%] | 95.1% | 0.2% | 0.088 [0.046, 0.143] |
| censored · HMM+FW |  | 1293 | 87.5% [80.4%, 93.9%] | 97.3% | 0.4% | 0.072 [0.038, 0.114] |
| censored · HMM+PSSM |  | 1293 | 93.5% [88.3%, 97.7%] | 97.3% | 0.0% | 0.040 [0.021, 0.063] |
| censored · HMM+PSSM+FW |  | 1293 | 94.6% [89.0%, 98.2%] | 98.9% | 0.2% | 0.030 [0.016, 0.050] |
| censored · HMM+PWM |  | 1293 | 95.7% [92.1%, 98.6%] | 95.7% | 0.1% | 0.028 [0.014, 0.046] |
| censored · HMM+PWM+FW |  | 1293 | 97.3% [93.7%, 99.5%] | 97.3% | 0.0% | 0.019 [0.010, 0.032] |
| censored · HMM+PSSM+PWM |  | 1293 | 97.3% [94.2%, 99.5%] | 97.3% | 0.1% | 0.018 [0.009, 0.030] |
| **censored · HMM+PSSM+PWM+FW** | current settings | 1293 | 97.8% [94.7%, 100.0%] | 97.8% | 0.0% | 0.013 [0.006, 0.025] |
| floor -20 · HMM |  | 1293 | 81.5% [71.5%, 89.7%] | 96.2% | 0.2% | 0.121 [0.071, 0.185] |
| floor -20 · HMM+FW |  | 1293 | 81.0% [71.0%, 89.3%] | 96.2% | 0.4% | 0.108 [0.064, 0.163] |
| floor -20 · HMM+PSSM |  | 1293 | 84.8% [76.6%, 91.0%] | 97.3% | 0.3% | 0.084 [0.051, 0.122] |
| floor -20 · HMM+PSSM+FW |  | 1293 | 84.8% [76.4%, 91.2%] | 97.8% | 0.4% | 0.081 [0.048, 0.120] |
| floor -20 · HMM+PWM |  | 1293 | 84.8% [75.2%, 92.5%] | 95.1% | 0.1% | 0.071 [0.042, 0.110] |
| floor -20 · HMM+PWM+FW |  | 1293 | 90.2% [82.5%, 95.2%] | 95.1% | 0.0% | 0.062 [0.036, 0.096] |
| floor -20 · HMM+PSSM+PWM |  | 1293 | 94.6% [90.0%, 97.8%] | 96.2% | 0.0% | 0.035 [0.020, 0.056] |
| floor -20 · HMM+PSSM+PWM+FW |  | 1293 | 95.1% [90.0%, 98.2%] | 97.3% | 0.0% | 0.032 [0.018, 0.055] |
| censored · HMM+PSSM · SIG/PRO/MAT/POST | previous ConoDictor database | 1293 | 93.5% [88.3%, 97.7%] | 96.7% | 0.1% | 0.040 [0.021, 0.063] |
| censored · HMM+PSSM+PWM+FW · SIG | signal only | 1293 | 83.2% [73.0%, 89.9%] | 93.5% | 0.1% | 0.071 [0.045, 0.104] |
| censored · HMM+PSSM+PWM+FW · PRO/MAT | without signal | 1293 | 88.6% [78.3%, 95.7%] | 94.6% | 0.0% | 0.082 [0.039, 0.147] |

#### complete precursor

| variant | note | n | positives correct [95% CI] | family if known | negatives called a family | log loss [95% CI] |
|---|---|---|---|---|---|---|
| censored · HMM |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| censored · HMM+FW |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| censored · HMM+PSSM |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| censored · HMM+PSSM+FW |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| censored · HMM+PWM |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.001] |
| censored · HMM+PWM+FW |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.001] |
| censored · HMM+PSSM+PWM |  | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| **censored · HMM+PSSM+PWM+FW** | current settings | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| floor -20 · HMM |  | 130 | 95.4% [88.2%, 100.0%] | 100.0% | - | 0.148 [0.025, 0.332] |
| floor -20 · HMM+FW |  | 130 | 95.4% [88.2%, 100.0%] | 100.0% | - | 0.149 [0.024, 0.330] |
| floor -20 · HMM+PSSM |  | 130 | 96.2% [90.6%, 100.0%] | 100.0% | - | 0.099 [0.014, 0.229] |
| floor -20 · HMM+PSSM+FW |  | 130 | 96.2% [90.6%, 100.0%] | 100.0% | - | 0.102 [0.012, 0.237] |
| floor -20 · HMM+PWM |  | 130 | 94.6% [87.1%, 99.3%] | 100.0% | - | 0.126 [0.034, 0.265] |
| floor -20 · HMM+PWM+FW |  | 130 | 95.4% [89.5%, 99.3%] | 100.0% | - | 0.127 [0.035, 0.262] |
| floor -20 · HMM+PSSM+PWM |  | 130 | 98.5% [95.6%, 100.0%] | 100.0% | - | 0.058 [0.009, 0.140] |
| floor -20 · HMM+PSSM+PWM+FW |  | 130 | 98.5% [95.6%, 100.0%] | 100.0% | - | 0.068 [0.011, 0.162] |
| censored · HMM+PSSM · SIG/PRO/MAT/POST | previous ConoDictor database | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.000 [0.000, 0.000] |
| censored · HMM+PSSM+PWM+FW · SIG | signal only | 130 | 100.0% [100.0%, 100.0%] | 100.0% | - | 0.001 [0.000, 0.002] |
| censored · HMM+PSSM+PWM+FW · PRO/MAT | without signal | 130 | 87.7% [75.0%, 96.9%] | 95.4% | - | 0.569 [0.136, 1.149] |

#### no signal peptide

| variant | note | n | positives correct [95% CI] | family if known | negatives called a family | log loss [95% CI] |
|---|---|---|---|---|---|---|
| censored · HMM |  | 50 | 52.0% [31.5%, 74.5%] | 82.0% | - | 1.864 [1.034, 2.711] |
| censored · HMM+FW |  | 50 | 54.0% [36.0%, 73.3%] | 90.0% | - | 1.499 [0.840, 2.231] |
| censored · HMM+PSSM |  | 50 | 76.0% [61.0%, 88.9%] | 90.0% | - | 0.762 [0.355, 1.306] |
| censored · HMM+PSSM+FW |  | 50 | 80.0% [62.5%, 93.1%] | 96.0% | - | 0.549 [0.227, 1.011] |
| censored · HMM+PWM |  | 50 | 84.0% [71.9%, 94.0%] | 84.0% | - | 0.572 [0.257, 0.970] |
| censored · HMM+PWM+FW |  | 50 | 90.0% [78.0%, 98.2%] | 90.0% | - | 0.355 [0.128, 0.670] |
| censored · HMM+PSSM+PWM |  | 50 | 90.0% [78.7%, 98.3%] | 90.0% | - | 0.339 [0.128, 0.616] |
| **censored · HMM+PSSM+PWM+FW** | current settings | 50 | 92.0% [81.1%, 100.0%] | 92.0% | - | 0.250 [0.061, 0.517] |
| floor -20 · HMM |  | 50 | 44.0% [23.6%, 66.1%] | 86.0% | - | 2.231 [1.321, 3.168] |
| floor -20 · HMM+FW |  | 50 | 42.0% [22.0%, 64.0%] | 86.0% | - | 1.901 [1.124, 2.730] |
| floor -20 · HMM+PSSM |  | 50 | 54.0% [36.7%, 72.1%] | 90.0% | - | 1.481 [0.846, 2.347] |
| floor -20 · HMM+PSSM+FW |  | 50 | 54.0% [35.9%, 71.4%] | 92.0% | - | 1.402 [0.759, 2.268] |
| floor -20 · HMM+PWM |  | 50 | 58.0% [37.0%, 79.5%] | 82.0% | - | 1.197 [0.608, 1.887] |
| floor -20 · HMM+PWM+FW |  | 50 | 76.0% [59.1%, 88.9%] | 82.0% | - | 0.957 [0.460, 1.636] |
| floor -20 · HMM+PSSM+PWM |  | 50 | 84.0% [69.4%, 94.0%] | 86.0% | - | 0.561 [0.242, 1.038] |
| floor -20 · HMM+PSSM+PWM+FW |  | 50 | 86.0% [71.7%, 96.4%] | 90.0% | - | 0.485 [0.162, 0.970] |
| censored · HMM+PSSM · SIG/PRO/MAT/POST | previous ConoDictor database | 50 | 76.0% [61.0%, 88.9%] | 88.0% | - | 0.770 [0.353, 1.317] |
| censored · HMM+PSSM+PWM+FW · SIG | signal only | 50 | 40.0% [21.4%, 58.3%] | 76.0% | - | 1.389 [0.993, 1.825] |
| censored · HMM+PSSM+PWM+FW · PRO/MAT | without signal | 50 | 90.0% [78.0%, 98.0%] | 92.0% | - | 0.342 [0.117, 0.708] |

## Misclassified positives (8 rows; all in `evaluation_misses.tsv`)

| method | id | label | predicted | p_true | best_family | best_family_share | type | length | max_identity | own_SIG_hmm | own_MAT_hmm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| map | D6C4J9 | T | NONE | 0.0926 | M | 0.687 | no signal peptide | 28 | 0.571 | -5.11 | 8.12 |
| map | P0DOW6 | O2 | O1 | 0.1162 | O1 | 0.866 | no signal peptide | 25 | 0.64 | -6.91 | 14.09 |
| map | P86256 | T | A | 0.154 | A | 0.533 | no signal peptide | 12 | 0.727 | -1.87 | 2.18 |
| map | P86257 | T | A | 0.0498 | A | 0.692 | no signal peptide | 12 | 0.75 | -1.87 | 1.93 |
| bayes | D6C4J9 | T | NONE | 0.1017 | M | 0.687 | no signal peptide | 28 | 0.571 | -5.11 | 8.12 |
| bayes | P0DOW6 | O2 | O1 | 0.2422 | O1 | 0.703 | no signal peptide | 25 | 0.64 | -6.91 | 14.09 |
| bayes | P86256 | T | A | 0.1322 | A | 0.523 | no signal peptide | 12 | 0.727 | -1.87 | 2.18 |
| bayes | P86257 | T | A | 0.0424 | A | 0.68 | no signal peptide | 12 | 0.75 | -1.87 | 1.93 |

## Notes

- Posterior between-family spread (tau) per feature: SIG:hmm 0.31, SIG:pssm 0.72, SIG:pwm 0.44, PRO:hmm 0.35, PRO:pssm 0.66, PRO:pwm 0.76, MAT:hmm 0.37, MAT:pssm 0.91, MAT:pwm 0.27, FW:framework 0.29; the penalised fit uses these posterior means.
