# conodictor-db-2026_04

ConoDictor reference database built by conodb 0.1.0 from UniProt release 2026_04.

| File | Content |
|---|---|
| `profiles.hmm` | Region HMMs (SIG/PRO/MAT/POST) per superfamily; GA cutoffs calibrated on negatives |
| `profiles.prf` | Region PSSMs (pftools generalized profiles); /CUT_OFF calibrated on negatives |
| `profiles.pwm.json` | Ungapped lightmotif windows per region |
| `precursors.hmm` + `precursors.labels.json` | Full-precursor HMMs and the region of each match state |
| `sites.json` | Cleavage-site windows (signal peptide end, mature start, mature end) |
| `framework.json` | Cysteine-framework counts per family |
| `model.json` | Classifier: conditional logit over A, B, C, Con-insulin, D, H, I1, I2, I3, J, L, M, O1, O2, O3, P, S, T + NONE; features SIG:hmm, SIG:pssm, SIG:pwm, PRO:hmm, PRO:pssm, PRO:pwm, MAT:hmm, MAT:pssm, MAT:pwm, FW:framework |
| `settings.json` | Settings needed to compute features exactly as in training |
| `cutoffs.tsv` | The calibrated cutoffs |
| `manifest.json` | Versions and sha256 of every file |
| `supplementary/` | Tables behind the evaluation (index in its README) |

Prediction must use conodb.core (shared with conodictor) so that features are
computed exactly as during training.
