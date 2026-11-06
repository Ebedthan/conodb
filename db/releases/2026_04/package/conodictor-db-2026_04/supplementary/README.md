# Supplementary material: ConoDictor database, UniProt 2026_04

| File | Content |
|---|---|
| `S01_report.md` | Database audit: superfamilies, tiers, entry groups at each identity, label curation, conflicts |
| `S02_families.tsv` | Entries, entry groups, species, tiers per superfamily and family |
| `S03_region_support.tsv` | Sequences and identity clusters per family and region |
| `S04_label_conflicts.tsv` | Near-identical sequences with different superfamily labels |
| `S05_unmapped_labels.tsv` | UniProt family labels without a curated mapping |
| `S06_partition.md` | Test set, cross-validation folds, exclusions |
| `S07_partition_summary.tsv` | Groups and entries per superfamily in train / test |
| `S08_negatives.tsv` | Every negative sequence: source, role, QC removals and reasons |
| `S09_profiles.tsv` | Every profile built: set, unit, sequences, match states |
| `S10_cutoffs.tsv` | Baseline per-profile HMM cutoffs |
| `S11_selection.md` | Choice of the defaults by grouped cross-validation on training rows (rule, intervals) |
| `S12_selection_classification.tsv` | Cross-validated classifier variants at every tau |
| `S13_selection_regions.tsv` | Cross-validated region-prediction variants |
| `S14_selection_region_errors.tsv` | Boundary errors of region prediction in training cross-validation |
| `S15_selection_recommended.json` | Settings recommended by cross-validation |
| `S16_fit.md` | Classifier settings, posterior summary, held-out results |
| `S17_bayes_summary.csv` | Full posterior summary of the Bayesian model (R-hat, ESS) |
| `S18_bayes_status.json` | Whether the Bayesian posterior was sampled, reused or adopted |
| `S19_predictions.tsv` | Held-out predictions of every classifier, with P(NONE) and best family |
| `S20_evaluation.md` | Held-out performance by sequence type, identity, calibration, misses, ablations |
| `S21_evaluation_metrics.tsv` | All evaluation metrics (long format) |
| `S22_evaluation_calibration.tsv` | Calibration bins |
| `S23_evaluation_misses.tsv` | Every misclassified held-out positive |
| `S24_evaluation_ablations.tsv` | Test-set results of every classifier variant, with intervals |
| `S25_anatomy.md` | Region prediction benchmark on the test set, every region variant |
| `S26_anatomy_metrics.tsv` | Region prediction metrics per region, variant, family source, subset |
| `S27_anatomy_sequences.tsv` | Per-sequence region prediction metrics (mature peptide, residue accuracy) |
| `S28_anatomy_errors.tsv` | Residues between predicted and true boundaries, per variant |
| `S29_anatomy_predictions.tsv` | Every predicted and true region boundary |
| `S30_manifest.json` | Build manifest: versions, settings, UniProt release, file checksums |
| `S31_labels.tsv` | Curated label mapping |
| `S32_conodb.toml` | Configuration used for this build |
