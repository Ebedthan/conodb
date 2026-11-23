# Benchmark: test set of 2026_04

ConoDictor 3 = the release's evaluation model (training libraries only). Intervals: bootstrap over entry groups; differences to ConoDictor 3 are paired (same sequences, same resamples). Other tools' labels mapped with `benchmark/labels.tsv`; CONFLICT (two superfamilies) counts as a conopeptide call without the right family; a label that is not a release superfamily (OTHER) counts as a wrong family.

## all positives

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 3 | 184 | 97.8% [94.4, 100.0] | 0.5% | 1.6% | - |
| ConoDictor 3 (Bayes) | 184 | 97.8% [94.5, 100.0] | 0.5% | 1.6% | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 184 | 95.7% [91.1, 98.8] | 2.7% | 1.6% | -2.2 [-4.7, -0.4] |
| ConoDictor 2.3.1 | 184 | 92.9% [87.2, 96.6] | 2.7% | 4.3% | -4.9 [-9.0, -1.7] |
| ConoDictor 2.3.1 (HMM only) | 184 | 94.0% [88.5, 97.5] | 6.0% | 0.0% | -3.8 [-7.2, -0.7] |
| ConoSorter | 184 | 90.8% [85.1, 94.7] | 3.8% | 5.4% | -7.1 [-12.3, -1.8] |
| ConoSorter (regex only) | 184 | 88.6% [82.2, 93.0] | 6.5% | 4.9% | -9.2 [-15.6, -4.0] |

## complete precursor

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 3 | 130 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | - |
| ConoDictor 3 (Bayes) | 130 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 130 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 | 130 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 (HMM only) | 130 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoSorter | 130 | 93.1% [87.9, 97.1] | 2.3% | 4.6% | -6.9 [-12.0, -3.0] |
| ConoSorter (regex only) | 130 | 93.1% [87.7, 96.7] | 2.3% | 4.6% | -6.9 [-12.0, -3.0] |

## signal, incomplete

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 3 | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | - |
| ConoDictor 3 (Bayes) | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 (HMM only) | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoSorter | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoSorter (regex only) | 4 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |

## no signal peptide

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 3 | 50 | 92.0% [80.9, 100.0] | 2.0% | 6.0% | - |
| ConoDictor 3 (Bayes) | 50 | 92.0% [81.0, 100.0] | 2.0% | 6.0% | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 50 | 84.0% [68.9, 95.7] | 10.0% | 6.0% | -8.0 [-17.0, -1.8] |
| ConoDictor 2.3.1 | 50 | 74.0% [54.3, 87.5] | 10.0% | 16.0% | -18.0 [-31.8, -6.3] |
| ConoDictor 2.3.1 (HMM only) | 50 | 78.0% [61.4, 90.3] | 22.0% | 0.0% | -14.0 [-24.5, -4.1] |
| ConoSorter | 50 | 84.0% [71.2, 92.7] | 8.0% | 8.0% | -8.0 [-21.4, +7.1] |
| ConoSorter (regex only) | 50 | 76.0% [60.9, 87.2] | 18.0% | 6.0% | -16.0 [-33.3, +2.1] |

## positives, families the tool covers

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 3 | 184 | 97.8% [94.3, 100.0] | 0.5% | 1.6% | - |
| ConoDictor 3 (Bayes) | 184 | 97.8% [94.5, 100.0] | 0.5% | 1.6% | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 184 | 95.7% [91.0, 99.0] | 2.7% | 1.6% | -2.2 [-4.8, -0.5] |
| ConoDictor 2.3.1 | 184 | 92.9% [87.2, 96.8] | 2.7% | 4.3% | -4.9 [-8.6, -1.8] |
| ConoDictor 2.3.1 (HMM only) | 184 | 94.0% [89.0, 97.5] | 6.0% | 0.0% | -3.8 [-7.0, -1.0] |
| ConoSorter | 184 | 90.8% [85.4, 94.7] | 3.8% | 5.4% | -7.1 [-12.6, -2.3] |
| ConoSorter (regex only) | 183 | 89.1% [83.0, 93.6] | 6.6% | 4.4% | -8.7 [-15.0, -3.4] |

## positives first public after 2021-07-07 (unseen by every tool)

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 3 | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | - |
| ConoDictor 3 (Bayes) | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 (HMM only) | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoSorter | 11 | 81.8% [55.6, 100.0] | 0.0% | 18.2% | -18.2 [-44.4, +0.0] |
| ConoSorter (regex only) | 11 | 81.8% [62.5, 100.0] | 0.0% | 18.2% | -18.2 [-44.4, +0.0] |

## all negatives

| tool | n | called a family [95% CI] | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|
| ConoDictor 3 | 1109 | 0.0% [0.0, 0.0] | - |
| ConoDictor 3 (Bayes) | 1109 | 0.0% [0.0, 0.0] | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 1109 | 4.1% [2.1, 6.5] | +4.1 [+2.2, +6.5] |
| ConoDictor 2.3.1 | 1109 | 28.0% [22.3, 33.9] | +28.0 [+22.4, +33.7] |
| ConoDictor 2.3.1 (HMM only) | 1109 | 2.3% [1.3, 3.5] | +2.3 [+1.3, +3.5] |
| ConoSorter | 1109 | 48.6% [42.7, 54.7] | +48.6 [+42.7, +54.4] |
| ConoSorter (regex only) | 1109 | 46.6% [40.6, 53.3] | +46.6 [+40.6, +52.1] |

## negative: shuffled

| tool | n | called a family [95% CI] | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|
| ConoDictor 3 | 131 | 0.0% [0.0, 0.0] | - |
| ConoDictor 3 (Bayes) | 131 | 0.0% [0.0, 0.0] | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 131 | 33.6% [24.1, 41.3] | +33.6 [+24.4, +41.6] |
| ConoDictor 2.3.1 | 131 | 80.9% [71.0, 89.1] | +80.9 [+70.7, +89.5] |
| ConoDictor 2.3.1 (HMM only) | 131 | 5.3% [2.2, 9.0] | +5.3 [+2.0, +9.2] |
| ConoSorter | 131 | 57.3% [43.8, 69.6] | +57.3 [+43.4, +69.9] |
| ConoSorter (regex only) | 131 | 51.1% [36.1, 67.7] | +51.1 [+35.3, +66.2] |

## negative: uniprot

| tool | n | called a family [95% CI] | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|
| ConoDictor 3 | 978 | 0.0% [0.0, 0.0] | - |
| ConoDictor 3 (Bayes) | 978 | 0.0% [0.0, 0.0] | +0.0 [+0.0, +0.0] |
| Nearest neighbour | 978 | 0.1% [0.0, 0.4] | +0.1 [+0.0, +0.3] |
| ConoDictor 2.3.1 | 978 | 20.9% [16.0, 26.6] | +20.9 [+15.7, +26.6] |
| ConoDictor 2.3.1 (HMM only) | 978 | 1.8% [0.8, 3.1] | +1.8 [+0.8, +3.2] |
| ConoSorter | 978 | 47.4% [41.0, 53.9] | +47.4 [+40.6, +53.4] |
| ConoSorter (regex only) | 978 | 46.0% [39.2, 52.6] | +46.0 [+39.1, +52.1] |

## positives first public after 2021-07-07 (unseen by this tool)

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoDictor 2.3.1 | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |
| ConoDictor 2.3.1 (HMM only) | 11 | 100.0% [100.0, 100.0] | 0.0% | 0.0% | +0.0 [+0.0, +0.0] |

## positives first public after 2013-10-01 (unseen by this tool)

| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |
|---|---|---|---|---|---|
| ConoSorter | 49 | 79.6% [65.9, 89.1] | 4.1% | 16.3% | -18.4 [-31.4, -9.6] |
| ConoSorter (regex only) | 49 | 73.5% [57.5, 85.7] | 12.2% | 14.3% | -24.5 [-40.0, -13.8] |

## Unmapped labels (counted as OTHER: a family the release does not have)

- ConoDictor 2.3.1: Conkunitzin (7)
- ConoDictor 2.3.1 (HMM only): Conkunitzin (12), Conopressin-Conophysin (3)
- ConoSorter: MTFLLLLVSV (2), Y2 (1), O4 (5), MKFPLLFISL (1), MKVAVVLLVS (1), N (1)
- ConoSorter (regex only): O4 (5)
