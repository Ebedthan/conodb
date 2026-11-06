# Partition: UniProt 2026_04

- Units: entry groups at 0.7 identity; a group is never split.
- Tier 1: 20% of groups held out as test, 5 grouped folds on the rest. Tier 2: leave-one-group-out. Tier 3: reference only.
- Temporal test: off.
- Seed: `conodictor-v3`. Previous release: none (fresh partition).

Entries by role: train 940, unlabelled 220, test 185, unmapped 81, reference 52, excluded 15

## Per superfamily

`fold_groups`: training groups in each fold (tier 1).

| superfamily | tier | groups | entries | train_groups | test_groups | temporal_test_groups | train_entries | test_entries | fold_groups |
|---|---|---|---|---|---|---|---|---|---|
| O1 | 1 | 91 | 284 | 73 | 18 | 0 | 239 | 45 | 15/15/15/14/14 |
| A | 1 | 62 | 195 | 50 | 12 | 0 | 144 | 51 | 10/10/10/10/10 |
| M | 1 | 54 | 178 | 43 | 11 | 0 | 127 | 51 | 9/9/9/8/8 |
| T | 1 | 38 | 141 | 30 | 8 | 0 | 123 | 18 | 6/6/6/6/6 |
| O2 | 1 | 23 | 97 | 18 | 5 | 0 | 87 | 10 | 4/4/4/3/3 |
| I1 | 1 | 16 | 38 | 13 | 3 | 0 | 31 | 7 | 3/3/3/2/2 |
| I2 | 1 | 14 | 34 | 11 | 3 | 0 | 31 | 3 | 3/2/2/2/2 |
| Con-insulin | 2 | 8 | 25 | 8 | 0 | 0 | 25 | 0 | logo |
| P | 2 | 7 | 8 | 7 | 0 | 0 | 8 | 0 | logo |
| C | 2 | 6 | 19 | 6 | 0 | 0 | 19 | 0 | logo |
| D | 2 | 6 | 32 | 6 | 0 | 0 | 32 | 0 | logo |
| B | 2 | 5 | 15 | 5 | 0 | 0 | 15 | 0 | logo |
| I3 | 2 | 5 | 8 | 5 | 0 | 0 | 8 | 0 | logo |
| L | 2 | 4 | 8 | 4 | 0 | 0 | 8 | 0 | logo |
| O3 | 2 | 4 | 22 | 4 | 0 | 0 | 22 | 0 | logo |
| S | 2 | 4 | 9 | 4 | 0 | 0 | 9 | 0 | logo |
| H | 2 | 3 | 4 | 3 | 0 | 0 | 4 | 0 | logo |
| J | 2 | 3 | 8 | 3 | 0 | 0 | 8 | 0 | logo |
| E | 3 | 2 | 3 | 0 | 0 | 0 | 0 | 0 |  |
| F | 3 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |  |
| Q | 3 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |  |
| U | 3 | 2 | 5 | 0 | 0 | 0 | 0 | 0 |  |
| B2 | 3 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |  |
| G | 3 | 1 | 3 | 0 | 0 | 0 | 0 | 0 |  |
| G2 | 3 | 1 | 2 | 0 | 0 | 0 | 0 | 0 |  |
| K | 3 | 1 | 6 | 0 | 0 | 0 | 0 | 0 |  |
| N | 3 | 1 | 3 | 0 | 0 | 0 | 0 | 0 |  |
| NSf-1 | 3 | 1 | 2 | 0 | 0 | 0 | 0 | 0 |  |
| R | 3 | 1 | 8 | 0 | 0 | 0 | 0 | 0 |  |
| V | 3 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |  |
| Y | 3 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |  |

## Excluded entries (15)

| accession | reason |
|---|---|
| A0A125S9F0 | label conflict (PRO: N x1, O2 x1) |
| A0A3G3C7S2 | label conflict (PRO: T x1, A x1) |
| B2KPN7 | label conflict (FULL: M x1, A x1) |
| C8CK74 | label conflict (PRO: N x1, O2 x1) |
| P0C6S2 | label conflict (FULL: L x1, R x1) |
| P0C6S3 | label conflict (MAT: R x2, L x1) |
| P0C8V8 | label conflict (MAT: O1 x1, M x1) |
| P0CH39 | label conflict (FULL: M x1, A x1) |
| P0DKQ3 | label conflict (FULL: M x1, O1 x1) |
| P0DKQ4 | label conflict (FULL: O1 x4, M x2) |
| P0DKQ5 | label conflict (FULL: O1 x4, M x2) |
| P0DMB3 | label conflict (MAT: T x2, M x1) |
| P0DPL8 | label conflict (SIG: T x2, A x1) |
| P84704 | label conflict (FULL: L x1, R x1) |
| Q9XZK2 | label conflict (FULL: M x1, O1 x1) |

## Tier changes after exclusions

_none_

## Groups mixing superfamily labels

- A0A2I6EDL3 (M, T)
- P0DQQ5 (A, O2)
- A6YR20 (M, O1)
- B0KZ79 (B, O2)
- B0L0Y5 (M, S)
- B2CI27 (A, S)
- B3FIA6 (M, V)
- B6V6L0 (A, M)
- A0A7M4DUE8 (A, M)
- C7DQB8 (I2, M)
- A0A125S9E4 (K, O2)
- D2Y495 (O1, T)
- H2BKJ0 (M, T)
- P0C904 (O1, O2)
- C0HLK1 (A, M)
- C0HKY6 (M, O2)
- C0HKY4 (M, T)
- P0C351 (A, O2)
- P0DM34 (G, T)
- P0DPL4 (A, T)
- P0DPM2 (A, O1)
- P0DPQ0 (C, O2)
- P58848 (O1, T)
- P0DP11 (M, S)

## Groups moved from train to test

_none_

