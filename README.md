# conodb

Reproducible builder for the ConoDictor reference database. One command
downloads the current UniProt data, extracts precursor regions, clusters
them, audits how much independent evidence each superfamily has, builds the
profiles, fits and evaluates the classifier, benchmarks region prediction,
and writes a versioned release: the package conodictor loads, a report of
what changed since the previous release, and a supplementary bundle with
every table behind the evaluation.

```bash
pip install ./conodb          # Python >= 3.10; pyopal, pyhmmer, pyfamsa, lightmotif, numpy, scipy
conodb init                   # writes conodb.toml + labels.tsv
conodb run                    # build / update the database
```

Rerunning is cheap: each stage is skipped when its inputs and the config
sections it reads are unchanged, and nothing is downloaded when UniProt has
no new release.

## Commands

| Command | What it does |
|---|---|
| `conodb init [dir]` | Example `conodb.toml` and `labels.tsv` |
| `conodb run [-c conodb.toml]` | Build or update. `--force` reruns all stages, `--from STAGE` reruns from a stage, `--until STAGE` stops early |
| `conodb run --input dump.jsonl.gz --release NAME` | Build from a local UniProt JSON-lines file (offline, reproducing a past build) |
| `conodb check` | Exit 0 if the latest UniProt release is built, **10** if an update is available |
| `conodb status` | Latest release, stage statistics, report path |
| `conodb stages` | Stages, the config sections they depend on, their outputs |

Exit codes: 0 success, 1 error, **2 labels need curation** (strict mode), 10 update available.

## Stages

| Stage | Reads config | Produces |
|---|---|---|
| `fetch` | `[source]` | `raw/uniprot.jsonl.gz`: every matching entry, verified against `X-Total-Results` |
| `parse` | `[labels]`, `[regions]`, labels file | `entries.tsv`, `regions.tsv`, `regions/{SIG,PRO,MAT,POST}.fasta`, `unmapped_labels.tsv`, `parse_warnings.tsv` |
| `cluster` | `[cluster]` | `clusters.tsv` (nested clusters per region and full sequence), `entry_groups.tsv` |
| `audit` | `[tiers]`, `[cluster]`, `[audit]` | `families.tsv`, `region_support.tsv`, `label_conflicts.tsv`, `diff_*.tsv`, **`report.md`** |
| `partition` | `[partition]`, `[tiers]`, `[cluster]`, previous release | `partition.tsv`, `partition_summary.tsv`, **`partition.md`** |
| `negatives` | `[negatives]`, `[partition]`, `[cluster]` | `negatives.tsv`, `negatives.fasta`, `raw/negatives.jsonl.gz` |
| `build` | `[build]`, `[partition]`, `[anatomy]`, pftools on PATH | `profiles.tsv`; libraries `profiles/{final,ship}` and holdouts: region `.hmm` (+ `.prf`), lightmotif `.pwm.json`, full-precursor `.full.hmm` + `.full.labels.json`, cleavage-site windows `.sites.json` |
| `score` | `[score]`, `[anatomy]`, pfsearch on PATH | `scores.tsv` (long, with hit coordinates), **`features.tsv`** (one row per sequence, one column per family x region x method, plus the framework feature) |
| `bayes` | `[fit]` (posterior settings only) | `fit/model_bayes.json`, `fit/bayes_draws.npz`, `fit/bayes_summary.csv`, `fit/bayes_status.json` |
| `select` | `[evaluate]`, `[fit]`, `[anatomy]` | **`selection.md`**, `selection_classification.tsv`, `selection_regions.tsv`, `selection_recommended.json`: defaults chosen by cross-validation on training rows |
| `fit` | `[fit]`, `[score]` | `cutoffs.tsv`, `fit/model_map.json` (evaluation model), `fit/model_ship.json` (refitted on all labelled data), `predictions.tsv`, **`fit.md`** |
| `anatomy` | `[anatomy]`, `[score]` | **`anatomy.md`**, `anatomy_metrics.tsv`, `anatomy_sequences.tsv`, `anatomy_errors.tsv`, `anatomy_predictions.tsv` |
| `evaluate` | `[evaluate]`, `[fit]` | **`evaluation.md`**, `evaluation_metrics.tsv`, `evaluation_misses.tsv`, `evaluation_calibration.tsv`, `evaluation_ablations.tsv` |
| `package` | `[package]` and everything above | `package/conodictor-db-<release>/` + `.tar.gz`: calibrated profiles, model, region models, **`supplementary/`** |

## Output layout

```
db/
  LATEST                 last successfully built release
  latest -> releases/<release>
  releases/<release>/    one folder per UniProt release (e.g. 2026_05)
    manifest.json        versions, config, UniProt release, sha256 of every output
    report.md            human summary: tiers, curation needs, QC, changes
    .stamps/             per-stage input hashes (incremental reruns)
    logs/
```

A release is only published (`LATEST` updated) after every stage succeeds.
Outputs are written atomically, so an interrupted run never leaves a file
that a later run would trust.

## How regions are assigned

From UniProt features (1-based, inclusive): `Signal` -> **SIG**; `Peptide`
-> **MAT** (`Chain` only when an entry has no `Peptide` and the chain does
not overlap the signal or a propeptide); the `Propeptide` immediately before
the first mature peptide -> **PRO**; propeptides after the last mature
peptide -> **POST**; anything else -> **INTER** (kept, not profiled).
Regions with non-`EXACT` boundaries (UniProt `UNSURE`/`OUTSIDE`) are kept
and flagged `exact=false`; later stages decide how to use them. Set
`regions.allow_uncertain = false` to drop them at parse time. Each region records its evidence class
(experimental > curated > similarity > predicted) and, for MAT, its cysteine
framework (e.g. `CC-C-C`).

## Labels

Labels come from UniProt's SIMILARITY comment, e.g. *"Belongs to the
conotoxin C superfamily. Consomatin family."* -> superfamily `C`, family
`Consomatin`. Built-in rules handle `conotoxin|conopeptide X superfamily` and
`conotoxin X family` and bare codes such as `O2 superfamily`. Everything else is reported in `unmapped_labels.tsv`
and the report until you curate it in `labels.tsv`:

```
raw_superfamily   raw_family  superfamily  family  scope        note
insulin family    *           Con-insulin          conopeptide  venom insulins
conoGAY family    *           ConoGAY              conopeptide
actinoporin family *          Actinoporin          other
```

`scope`: `conopeptide` (used to build the classifier), `other` (kept in the
audit; later used as hard negatives), `exclude` (dropped).

`labels.taxon_scope` (e.g. `["Conidae"]`) restricts classifier classes to
those taxonomic families: conopeptide-labelled entries from other taxa
(turrids, terebrids...) get scope `reference`. They stay in the database as
reference sequences but are never classes or negatives. The report's
"Entries by taxonomic family" table shows the effect.

Single mislabelled entries are corrected in `entry_labels.tsv`
(`labels.entry_overrides`), one row per accession with the reason; these win
over both UniProt and `labels.tsv`. Overrides always win over rules. With
`labels.strict = true`, a new uncurated family stops the run (exit 2) after
writing the curation list, so nothing enters the database unreviewed.

## Clustering, groups and tiers

Identity = identical aligned residues / length of the shorter sequence
(local alignment, pyopal), so a mature-peptide fragment has identity 1.0 to
its precursor. Clusters are computed at each `cluster.identities` threshold,
highest first, and are nested.

**Entry groups** are full-sequence clusters at `cluster.group_identity`.
They are the independent units for train/test splitting and for support:

- tier 1: >= `tiers.tier1_min_groups` groups: full profiles + k-fold calibration
- tier 2: >= `tiers.tier2_min_groups` groups: signal-first profiles, strong shrinkage
- tier 3: fewer: reference sequences only, not a classifier class

`region_support.tsv` gives clusters per region at every threshold, which is
the evidence for choosing `group_identity` (signal peptides are conserved
within a superfamily; full precursors much less so).

**Label conflicts**: sequences identical at the highest threshold but
labelled with different superfamilies, listed in `label_conflicts.tsv` most
severe first (FULL, then SIG, MAT, PRO), with the minority labels to review.
Regions shorter than `audit.conflict_min_length` (default 10 aa) are
ignored, since short C-terminal propeptides such as `GKR` are identical
across families by chance.

The report also shows entry groups and tiers at every identity threshold, so
the grouping threshold can be chosen from the data.

## Partition (test set and folds)

Units are entry groups; a group is never split between training and test.
Entries under review in `label_conflicts.tsv` are left out
(`partition.exclude_conflicts`; in a tie every tied entry is left out), and
tiers are recomputed without them.

| Tier | Role | Cross-validation |
|---|---|---|
| 1 | `test_fraction` of groups -> `test`, rest -> `train` | `folds` grouped folds, balanced by group count |
| 2 | `train` | leave-one-group-out (`fold = logo`) |
| 3 | `reference` | not evaluated |

Other roles: `excluded`, `negative` (scope `other`), `reference` (outside the
taxon scope), `unlabelled`, `unmapped`. With `partition.temporal_cutoff`,
groups first made public on or after that date form a `temporal_test` set.

**Stable across releases.** Each group inherits the role and fold its
members had in the previous release (a group containing any former test
entry stays test). Only new groups are assigned, ordered by a hash of their
smallest accession and `partition.seed`, so updating the database adds data
without reshuffling the benchmark. A training group is moved to test only
when a superfamily's quota cannot be met with new groups, and
`partition.md` lists every such move. Changing `seed` draws a new benchmark.
`partition.tsv` must therefore be kept with each release (it is not in the
generated `.gitignore`).

## Build (region profiles)

Profiles are built from `train` entries only, per superfamily and region
(`CONOPEP_<superfamily>_<REGION>`, as in ConoDictor):

- **final**: all training entries -> `profiles/final.hmm` (and `final.prf`)
- **holdout**: the same family rebuilt without one cross-validation unit,
  each fold (tier 1) or each entry group (tier 2) ->
  `profiles/holdout/<superfamily>/<unit>.hmm|.prf`

An entry's out-of-fold features (next stage) come from its own family's
holdout profiles for its unit plus every other family's final profiles,
which never contain it. The stage refuses to build if a profile would
contain a non-training entry or its own held-out unit.

Pipeline: identical sequences collapsed -> FAMSA alignment (pyfamsa) ->
HMM (pyhmmer, no build date in the file, so rebuilds are byte-identical)
and, when pftools is installed, `pfw -m` (sequence weights) -> `pfmake -m`
with `build.matrix` (blosum45.cmp, found in the pftools data folders or
given as a path). Profiles are cached by their member sequences and
settings: an update only rebuilds families whose data changed.
Cutoffs (`GA` for HMMs, `/CUT_OFF` for PSSMs) are written at the `package` stage.

Libraries:

- **final** (training entries): used for every evaluation number.
- **ship** (train + test + temporal test, labelled entries): what the package
  ships, once evaluation is done. Never used to score anything evaluated.
- **holdout**: as above, for out-of-fold training features.

Besides the region profiles, each library has:

- **PWMs** (`build.pwm`): lightmotif position-weight matrices over the
  region alignment's columns occupied in >= `anatomy.min_occupancy` of the
  sequences (gaps replaced by the column consensus). Ungapped, so they are an
  ablation against the pftools PSSMs, not a replacement.
- **Full-precursor HMMs** (`build.full_models`): one per family, from its
  complete precursors aligned whole, with match states = columns occupied
  >= `anatomy.min_occupancy` and each match state labelled with the region
  most of its residues belong to (`.full.labels.json`).
- **Cleavage-site windows**: `anatomy.site_left` residues before and
  `site_right` after each exact boundary (signal cleavage, mature start,
  mature end), pooled over families; lightmotif matrices are built from them
  when there are at least 5 windows per site type.

## Negatives

Sequences that are not conopeptides, for calibration and for the
"none of the known families" class:

| Source | What |
|---|---|
| `uniprot` | second query (`negatives.query`): non-venom *Conus* + small (30–300 aa) secreted metazoan proteins, toxins and venom proteins excluded; capped at `negatives.max_entries`, chosen by accession hash |
| `venom_other` | main-set venom proteins with scope `other` (CRISP, PLA2, actinoporin...) |
| `shuffled` | each complete precursor shuffled within its regions, cysteines kept in place |

Real negatives at least `negatives.max_identity` identical to a conopeptide
are removed (likely unannotated toxins; listed with the nearest conopeptide).
Real negatives are clustered and split into train/test like the positives,
stable across releases; shuffled ones take their source precursor's side, so
nothing derived from a test precursor is trained on. The negatives download
is cached per (query, UniProt release); offline builds use only the other two
sources.

## Score (out-of-fold features)

Every positive (train/test/temporal_test) and negative (train/test) is scored
against every family's final profiles. Training positives are also scored
against their own family's holdout profiles for their unit, and that score is
the one used in `features.tsv`, so no feature comes from a profile that saw
the sequence.

- HMM: pyhmmer, fixed Z (E-values unused), acceleration filters off
  (`score.max_mode`), reporting threshold `score.hmm_min_bits` (default
  -1000, i.e. negative scores are kept); feature = best domain bit score.
- PSSM: `pfsearch -f -k -a -r` raw score; the stage stops if pfsearch's output
  cannot be matched to every sequence.
- PWM (`score.pwm`): best lightmotif log-odds over all windows.
- Framework (`score.framework`): the cysteine framework of the **predicted
  mature region** (region prediction below, method `anatomy.method`) for
  each candidate family, scored as log((count + a) / (n + a V)) in that
  family's training sequences (a = `score.framework_pseudocount`). Counts and
  region models are out-of-fold exactly like the profiles: a training
  sequence's own family uses its holdout unit's full-precursor model and
  counts.
- `features.tsv` values: a score; empty = below the floor; `NA` = no profile
  available for this sequence (e.g. a tier-2 holdout with one sequence left).

## Fit (classifiers)

Trained on `train` rows of `features.tsv`, evaluated on `test` and
`temporal_test`. Classes: every family with final profiles, plus `NONE`.

| Method | What | Needs |
|---|---|---|
| `baseline` | per-profile HMM cutoffs from training negatives, set for a `fit.baseline_fpr` false-positive rate per sequence (Bonferroni over profiles); highest summed excess over the cutoffs wins, `NONE` if none passes | nothing |
| `map` | conditional-logit model below, penalised fit with fixed shrinkage `fit.tau_fixed` (numpy/scipy, seconds) | nothing |
| `bayes` | the same model, hierarchical, full posterior in Stan; predictions average the posterior | `pip install "./conodb[bayes]"` + CmdStan |

Model (`src/conodb/stan/conodictor_mnl.stan`), for sequence *i* and family *k*:

    eta[i,k] = alpha[k] + sum_p X[i,k,p] * B[k,p]      eta[i,NONE] = 0
    B[k,p]   = beta[p] + tau[p] * delta[k,p],  delta ~ N(0, 1)

`X[i,k,p]` is the out-of-fold score of sequence *i* on family *k*'s profile
*p* (region x method, `fit.regions` x `fit.methods`; POST is scored but not
used by default), standardised per *p*. HMM scores are clipped below at
`fit.hmm_floor` bits (default 0: censoring, "no hit is no hit"; -20 keeps
weak negative evidence and was clearly worse on the 2026_04 data), and a
missing non-HMM value takes the minimum seen in the fitting rows. The same
transformation is applied by `core.apply_model`, which conodictor uses. `beta` is shared across
families (a high score on a family's own profile is evidence for it), `tau`
is how much families depart from it, so rare families shrink to the shared
weights. Case weights are 1 / entry-group size. Cells whose holdout profile
could not be built get their own term, unused at prediction.

**Parallelism.** The Stan fit uses all cores (`fit.cores = 0`): chains run
side by side, and leftover cores become threads inside each chain, where
the log-likelihood is split over slices of sequences with `reduce_sum`
(the model is compiled with `STAN_THREADS`). 16 cores and 4 chains run as
4 chains x 4 threads; 2 cores and 4 chains as 2 chains at a time x 1 thread.
`fit.grainsize` sets the slice size (0 = automatic). Thread scheduling
changes the order of floating-point sums, so draws are not bit-identical
across core counts; posterior summaries agree.

**Shrinkage.** With `fit.tau_from_posterior` (default), the penalised fit
uses the posterior mean of each feature's `tau` instead of the single
`fit.tau_fixed`. The Bayesian and penalised fits agree closely, so the
shipped model is the penalised one: `fit/model_ship.json` is refitted on
every labelled sequence (train + test + temporal test, `fit.ship_all_data`)
with the posterior `tau`; `fit/model_map.json` (training rows only) is the
one every evaluation number comes from.

**Output.** For each sequence: probabilities over the families and NONE,
the prediction, `p_none`, and `best_family_conditional`: the best family
when NONE is excluded, with its probability renormalised (for sequences
known to be conopeptides, e.g. from peptidomics).

`fit.bayes = "auto"` samples when CmdStan is installed (install it with
`python -m cmdstanpy.install_cmdstan`), `"false"` keeps rebuilds fast.
`fit.md` reports R-hat, bulk ESS and divergences, and flags the run if
R-hat >= 1.01 or any divergence occurs. Editing the Stan file triggers a refit.

## When the Bayesian fit reruns

Sampling takes up to an hour, so it has its own stage with its own key: the
exact training design the model sees (a hash of the standardised feature
matrix, labels and weights after `methods`, `regions`, `hmm_floor`, group
weighting and negative balancing are applied), the Stan file, and the
sampler settings (the three prior scales, `chains`, `iter_warmup`,
`iter_sampling`, `adapt_delta`, `seed`). Columns the model does not use
(e.g. POST, PWM when not in `fit.methods`) can change without a resample. Changing `cores`, `grainsize`, `cmdstan`, `tau_fixed` or
`baseline_fpr`, moving to another machine, upgrading conodb, or rerunning
with `--from fit` / `--from bayes` never resamples: a stored posterior whose
key matches is reused (`fit/bayes_status.json`: `reused`). A posterior
sampled by the earlier single `fit` stage is adopted when that stage's stamp
matches the current data and settings (`adopted`). Only `--force`, or a
change to the data or the settings above, samples again.

## Evaluate

On the test (and temporal test) sets, without refitting the Bayesian model:

- by sequence type: complete precursors, signal present but incomplete,
  no signal peptide (mature peptides, pro + mature fragments), negatives by source
- positives by highest identity to a training conopeptide (`evaluate.identity_bands`)
- calibration (confidence bins vs accuracy), Brier score, log loss
- `family if known`: the best family when NONE is excluded, for sequences
  known to be conopeptides (e.g. peptidomics)
- every misclassified positive, with its own-family scores and nearest training sequence
- the classifier variant grid (below) on the test set, penalised fits, with
  bootstrap intervals over entry groups; the variant recommended by the
  select stage is marked. Reported after the choice, never used for it.

## Select (choosing the defaults)

The defaults are chosen by grouped cross-validation on **training rows
only**, by a rule fixed in advance, so the test set is evaluated once with
settings it did not influence.

Grids (`src/conodb/variants.py`, shared with evaluate and anatomy):

- classifier: HMM floor {censored at 0, -20} x {HMM, HMM+PSSM, HMM+PWM,
  HMM+PSSM+PWM} x {without, with framework}, regions SIG/PRO/MAT; plus the
  previous ConoDictor database (censored HMM+PSSM with POST), signal only,
  without signal, and the current settings;
- regions: envelope (no sites / all sites), and transfer with extension
  {within, cross} x overhang {0, 5} x sites {none, S, MS, S+MS, S+MS+ME}.

Classifier: 5 grouped folds (tier-1 partition folds; tier-2 groups and
negatives by a stable hash of their group; shuffled negatives follow their
source precursor); each variant is refitted per fold, at every
`evaluate.select_taus`. Metric: log loss weighted as in the fit. Regions:
every training conopeptide is cut with its own family's holdout models and
site matrices rebuilt without its unit's windows. Metric: `boundary_exact`
(signal end, mature start, mature end exact, pooled).

Rule: a variant is eligible when its paired bootstrap difference to the best
(`evaluate.bootstrap` resamples of entry groups) has a 95% interval
including zero, overall **and** on sequences without a signal peptide (so
mature peptides, a minority, cannot be traded away); the simplest eligible
variant wins. The classifier choice is made at `fit.tau_fixed`; the other
taus show whether it depends on shrinkage, which is why the Bayesian fit is
not needed to choose: run with `fit.bayes = "false"`, set the recommended
settings (printed as TOML in `selection.md`), then run once with Bayes for
the final numbers and the shipped model's posterior shrinkage.

## Anatomy (region prediction benchmark)

Every held-out conopeptide (test, temporal test) is cut into SIG/PRO/MAT/POST
with profiles built from training entries only, by every region variant of
the grid:

| Part | How |
|---|---|
| `envelope` | each region HMM's best envelope (score > 0), put in order, overlaps split |
| `transfer` | alignment to the family's full-precursor HMM; residues take the label of their match state |
| `/within`, `/cross` | unaligned ends extended ungapped over the remaining states, within the region of the last aligned state (a mature peptide alone does not grow a propeptide) or across regions |
| `/ohN` | N residues past where the extension stops keep its label (`anatomy.terminal_overhang`) |
| `+ sites S/MS/ME` | signal end / mature start / mature end refined with the lightmotif cleavage-site matrices (`anatomy.sites`); only where residues exist on both sides |

Refinement searches +/- `anatomy.search_radius` residues around each
boundary for the best site score minus `anatomy.distance_penalty` per
residue shifted. Each variant is run with the true family (region prediction
alone) and with the family predicted by the classifier (end to end). POST
is predicted and scored even though it is not a classification feature.

Metrics, against UniProt (coordinates only where the boundary is EXACT):
per region detection, false-region rate, start/end exact, within 1/2/5,
mean absolute error; per sequence exact mature peptide and per-residue
accuracy; `boundary_exact`; overall, by sequence type and by family; and
`anatomy_errors.tsv`, the residues between predicted and true boundaries
(e.g. a C-terminal G left in the mature peptide). `[anatomy]` chooses what
the framework feature and conodictor use: set it from `selection.md`, not
from this test-set table.

## Package

`package/conodictor-db-<release>/` (and a `.tar.gz`) is everything
conodictor needs, built from the **ship** libraries:

| File | What |
|---|---|
| `profiles.hmm`, `profiles.prf`, `profiles.pwm.json` | region profiles; HMM `GA` and PSSM `/CUT_OFF` set from the training negatives at `package.cutoff_fpr` per sequence (Bonferroni over profiles) |
| `precursors.hmm`, `precursors.labels.json` | full-precursor models and their state labels |
| `sites.json`, `framework.json` | cleavage-site windows, framework counts per family |
| `model.json`, `settings.json` | `fit/model_ship.json`; every setting `core.apply_model` and region prediction need |
| `cutoffs.tsv`, `manifest.json`, `README.md` | cutoffs, versions and sha256 of every file |
| `supplementary/` | `S01_...`, `S02_...`: every table behind the reports (labels, tiers, partition, negatives, profiles, evaluation, ablations, calibration, misses, anatomy, Bayesian diagnostics) with the config used, indexed in its `README.md` |

## Shared code

`src/conodb/core.py` holds everything prediction needs (HMM search with
alignments, region prediction, cleavage sites, the framework feature,
`apply_model`). conodictor v3 imports it, so the tool applies exactly what
was evaluated; a test checks that `apply_model` on raw features reproduces
the evaluation predictions.

## Automating updates

`.github/workflows/update-db.yml` runs `conodb check` weekly and, when
UniProt has a new release, runs the build and opens a pull request with the
report. UniProt releases roughly every 8 weeks.

## Tests

```bash
pip install -e "./conodb[test]" && pytest conodb
```

The tests run the full pipeline against a local mock of the UniProt REST API
(pagination, release headers, transient errors, truncated downloads).
