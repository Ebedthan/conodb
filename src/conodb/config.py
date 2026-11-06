"""Configuration: one TOML file describes a whole database build.

Every value that can change a result lives here, so a build is fully
determined by (config, UniProt release, conodb version). Each stage hashes
only the config sections it reads, so editing e.g. tier thresholds reruns
the audit without re-downloading or re-clustering.
"""

from __future__ import annotations

import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

DEFAULT_QUERY = (
    "taxonomy_id:37797 AND (cc_tissue_specificity:venom OR cc_scl_term:SL-0177) AND reviewed:true"
)


@dataclass
class SourceConfig:
    query: str = DEFAULT_QUERY
    base_url: str = "https://rest.uniprot.org/uniprotkb/search"
    page_size: int = 500
    timeout: float = 60.0
    retries: int = 6


@dataclass
class LabelsConfig:
    # Curated TSV mapping raw UniProt family text -> labels and scope.
    overrides: Path | None = None
    # Per-entry corrections (accession -> superfamily/family/scope), for
    # mislabelled entries found e.g. through label conflicts.
    entry_overrides: Path | None = None
    # Taxonomic families (as in the UniProt lineage, e.g. "Conidae") whose
    # entries may be classifier classes. Conopeptide-labelled entries from
    # other taxa get scope 'reference': kept as reference sequences, never
    # used as classes or as negatives. Empty list = no restriction.
    taxon_scope: list[str] = field(default_factory=list)
    # Fail the run when a raw family string matches no rule/override, so a
    # human curates new families before they silently enter the database.
    strict: bool = False


@dataclass
class RegionsConfig:
    # Keep regions whose UniProt boundaries are not EXACT (UNSURE/OUTSIDE).
    # They are recorded with exact=false; later stages decide how to use
    # them (fine for profile training, excluded from boundary evaluation).
    allow_uncertain: bool = True
    # Minimum length (aa) for a region to be recorded. Keep this low: short
    # C-terminal propeptides (e.g. 'GKR') are real; length filters for
    # profile building belong to the build stage.
    min_length: int = 1


@dataclass
class ClusterConfig:
    # Identity thresholds (identical residues / length of the shorter
    # sequence, CD-HIT style) reported in the audit.
    identities: list[float] = field(default_factory=lambda: [0.9, 0.7, 0.5])
    # Threshold used to define entry groups (units for train/test splits and
    # for support tiers). Must be one of `identities`.
    group_identity: float = 0.7
    threads: int = 0


@dataclass
class TiersConfig:
    # Tiers count independent entry groups (see ClusterConfig.group_identity)
    tier1_min_groups: int = 10
    tier2_min_groups: int = 3


@dataclass
class AuditConfig:
    # Label-conflict check: which sequence kinds to compare, and the minimum
    # length for a region to count (short regions such as C-terminal 'GKR'
    # propeptides are identical across families by chance).
    conflict_kinds: list[str] = field(default_factory=lambda: ["FULL", "SIG", "MAT", "PRO"])
    conflict_min_length: int = 10


@dataclass
class PartitionConfig:
    test_fraction: float = 0.2
    folds: int = 5
    # Salt for the hash that orders groups; change it only to draw a new
    # benchmark (it breaks comparability with earlier releases).
    seed: str = "conodictor-v3"
    # Leave out entries under review in label_conflicts.tsv.
    exclude_conflicts: bool = True
    # "YYYY-MM-DD": groups first public on/after this date form a temporal
    # test set (tier 1-2 superfamilies). "" = off.
    temporal_cutoff: str = ""


@dataclass
class NegativesConfig:
    # Second UniProt query: real non-conopeptide proteins. Default: non-venom
    # Conus proteins plus secreted gastropod proteins, toxins excluded.
    # Small secreted metazoan proteins (precursor-sized, signal peptides) plus
    # non-venom Conus proteins; toxins and venom proteins excluded.
    query: str = ("(taxonomy_id:6490 OR (taxonomy_id:33208 AND cc_scl_term:SL-0243 AND length:[30 TO 300])) "
                  "AND reviewed:true NOT cc_tissue_specificity:venom NOT keyword:KW-0800")
    enabled: bool = True
    include_venom_other: bool = True     # scope 'other' venom proteins of the main set
    shuffles: int = 1                    # shuffled copies per complete precursor (0 = none)
    min_length: int = 30
    max_length: int = 1000
    # Real negatives this identical to a conopeptide are removed (likely toxins).
    max_identity: float = 0.9
    # Cap on real UniProt negatives kept (deterministic choice by accession
    # hash, stable as UniProt grows). 0 = no cap.
    max_entries: int = 5000


@dataclass
class AnatomyConfig:
    """Precursor region prediction (build, score, the anatomy benchmark, conodictor)."""
    # Full-precursor HMMs: an alignment column becomes a match state (and
    # carries a region label) when at least this fraction of sequences has a
    # residue there.
    min_occupancy: float = 0.5
    # Predicted regions shorter than this are dropped.
    min_run: int = 3
    # Label transfer: unaligned ends are extended over the remaining model
    # states within the region of the last aligned state ("within") or
    # across regions ("cross"); where the extension stops, up to
    # terminal_overhang more residues keep the label reached.
    extension: str = "within"
    terminal_overhang: int = 5
    # Cleavage-site matrices (lightmotif): residues before / after a boundary.
    site_left: int = 5
    site_right: int = 3
    # Boundary refinement: search +/- radius residues around the predicted
    # boundary; score = site log-odds - distance_penalty x |shift|.
    search_radius: int = 6
    distance_penalty: float = 0.5
    pseudocount: float = 0.5
    # Without a mature-region prediction, sequences up to this length are
    # treated as mature peptides (framework feature).
    short_peptide_length: int = 60
    # Region-prediction method used for the framework feature and by
    # conodictor: "transfer" or "envelope", each with or without sites.
    method: str = "transfer"
    # Boundaries refined with the cleavage-site matrices: any of "signal",
    # "mature_start", "mature_end" (true = all three, false = none).
    sites: list[str] | bool = field(default_factory=lambda: ["signal", "mature_start"])
    # C-terminal processing: a mature peptide that ends the sequence loses a
    # trailing G[KR]{0,3} (amidation glycine and basic residues).
    cterm_trim: bool = False

    @property
    def site_types(self) -> tuple[str, ...]:
        if self.sites is True:
            return ("signal", "mature_start", "mature_end")
        if self.sites is False:
            return ()
        return tuple(t for t in ("signal", "mature_start", "mature_end") if t in self.sites)

    def region_kwargs(self) -> dict:
        """Keyword arguments of core.predict_segments for these settings (sites passed separately)."""
        return {"method": self.method, "radius": self.search_radius, "penalty": self.distance_penalty,
                "min_run": self.min_run, "overhang": self.terminal_overhang, "extension": self.extension,
                "site_types": self.site_types, "cterm": self.cterm_trim}


@dataclass
class BuildConfig:
    regions: list[str] = field(default_factory=lambda: ["SIG", "PRO", "MAT", "POST"])
    # Regions shorter than this are not used for profiles (e.g. 'GKR').
    min_length: int = 5
    # Minimum number of distinct sequences to build a profile.
    min_unique: int = 2
    # Use regions whose UniProt boundaries are not EXACT.
    use_uncertain: bool = True
    # PSSMs with pftools: "auto" (if pfw/pfmake are on PATH), "true" (required), "false".
    pssm: str = "auto"
    pfw: str = "pfw"
    pfmake: str = "pfmake"
    # Score matrix shipped with pftools: a path, or a file name searched in the
    # usual pftools data folders (PFTOOLS_DATA, <prefix>/share, /usr/share...).
    matrix: str = "blosum45.cmp"
    pfw_options: list[str] = field(default_factory=list)
    pfmake_options: list[str] = field(default_factory=list)
    # Full-precursor HMMs with per-state region labels (region prediction).
    full_models: bool = True
    # Ungapped position weight matrices (lightmotif) from each region
    # alignment, compared with the gapped pftools PSSMs.
    pwm: bool = True
    threads: int = 0


@dataclass
class ScoreConfig:
    # Reporting floor in bits. Very low: scores are stored uncensored and the
    # fit decides the floor (fit.hmm_floor); HMMER still reports no domain
    # for some unrelated sequences, which count as below the floor.
    hmm_min_bits: float = -1000.0
    # Turn HMMER's acceleration filters off so short profiles are always scored.
    max_mode: bool = True
    # PSSM scoring with pftools pfsearch: "auto" | "true" | "false".
    pssm: str = "auto"
    pfsearch: str = "pfsearch"
    pfsearch_options: list[str] = field(default_factory=list)
    # lightmotif PWM scores (needs build.pwm) and the cysteine-framework feature.
    pwm: bool = True
    framework: bool = True
    framework_pseudocount: float = 0.5
    threads: int = 0


@dataclass
class FitConfig:
    # Feature methods used by the classifiers.
    methods: list[str] = field(default_factory=lambda: ["hmm", "pssm", "framework"])
    # Regions whose profile scores are features (POST adds nothing: dropped).
    regions: list[str] = field(default_factory=lambda: ["SIG", "PRO", "MAT"])
    # HMM scores below this (or not reported) are set to it. 0 = the
    # censoring of earlier versions; a negative floor keeps weak evidence.
    hmm_floor: float = 0.0
    weight_by_group: bool = True        # case weight 1 / entry-group size
    # Give shuffled and real negatives the same total weight, so a large
    # synthetic set does not decide what "NONE" looks like.
    balance_negative_sources: bool = True
    baseline_fpr: float = 0.01          # baseline cutoffs: false-positive rate per sequence
    # Penalised fit (fast, always run): shrinkage per feature from the
    # Bayesian posterior when one exists, else tau_fixed.
    tau_fixed: float = 0.5
    tau_from_posterior: bool = True
    # Shipped model: the penalised fit refitted on every labelled sequence
    # (train + test), after evaluation on the held-out ones.
    ship_all_data: bool = True
    # Priors (standardised scores), shared by both fits
    beta_prior_sd: float = 2.5
    tau_prior_sd: float = 1.0
    alpha_prior_sd: float = 2.5
    # Bayesian fit with Stan: "auto" (if CmdStan is installed) | "true" | "false"
    bayes: str = "auto"
    cmdstan: str = ""                   # CmdStan path; "" = cmdstanpy default / CMDSTAN
    chains: int = 4
    iter_warmup: int = 1000
    iter_sampling: int = 1000
    adapt_delta: float = 0.9
    seed: int = 20261105
    # Cores for the Stan fit: chains run in parallel, leftover cores become
    # threads inside each chain (reduce_sum). 0 = all cores.
    cores: int = 0
    # Sequences per reduce_sum slice; 0 = automatic (N / (4 x threads)).
    grainsize: int = 0


@dataclass
class EvaluateConfig:
    # Upper edges of the identity-to-training bands (last band ends at 1).
    identity_bands: list[float] = field(default_factory=lambda: [0.5, 0.6, 0.7, 0.8, 0.9])
    # Upper edges of the confidence bins for calibration.
    calibration_bins: list[float] = field(default_factory=lambda: [0.5, 0.7, 0.9, 0.99])
    # Penalised refits of the variant grid on the test set (seconds each).
    ablations: bool = True
    # Selection stage: grouped cross-validation on training rows (classifier
    # variants, region-prediction variants) and a recommended configuration.
    select: bool = True
    bootstrap: int = 1000                 # resamples of entry groups for intervals
    select_taus: list[float] = field(default_factory=lambda: [0.25, 0.5, 1.0])  # tau sensitivity
    seed: int = 20261106


@dataclass
class PackageConfig:
    # Calibrated cutoffs written into the shipped profiles: false-positive
    # rate per sequence on all negatives (Bonferroni over profiles).
    cutoff_fpr: float = 0.01
    archive: bool = True                # also write a .tar.gz


@dataclass
class OutputConfig:
    root: Path = Path("db")


@dataclass
class Config:
    source: SourceConfig = field(default_factory=SourceConfig)
    labels: LabelsConfig = field(default_factory=LabelsConfig)
    regions: RegionsConfig = field(default_factory=RegionsConfig)
    cluster: ClusterConfig = field(default_factory=ClusterConfig)
    tiers: TiersConfig = field(default_factory=TiersConfig)
    audit: AuditConfig = field(default_factory=AuditConfig)
    partition: PartitionConfig = field(default_factory=PartitionConfig)
    negatives: NegativesConfig = field(default_factory=NegativesConfig)
    build: BuildConfig = field(default_factory=BuildConfig)
    score: ScoreConfig = field(default_factory=ScoreConfig)
    fit: FitConfig = field(default_factory=FitConfig)
    evaluate: EvaluateConfig = field(default_factory=EvaluateConfig)
    anatomy: AnatomyConfig = field(default_factory=AnatomyConfig)
    package: PackageConfig = field(default_factory=PackageConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    path: Path | None = None

    def section(self, name: str) -> dict:
        """Plain dict of one section (used for stage input hashes)."""
        d = asdict(getattr(self, name))
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in d.items()}

    def validate(self) -> None:
        c = self.cluster
        if not c.identities or any(not 0 < t <= 1 for t in c.identities):
            raise ValueError("cluster.identities must be values in (0, 1]")
        if c.group_identity not in c.identities:
            raise ValueError("cluster.group_identity must be one of cluster.identities")
        if self.tiers.tier2_min_groups > self.tiers.tier1_min_groups:
            raise ValueError("tiers.tier2_min_groups must be <= tiers.tier1_min_groups")
        for name in ("overrides", "entry_overrides"):
            path = getattr(self.labels, name)
            if path is not None and not path.exists():
                raise FileNotFoundError(f"labels.{name} not found: {path}")
        p = self.partition
        if not 0 < p.test_fraction < 1 or p.folds < 2:
            raise ValueError("partition.test_fraction must be in (0, 1) and partition.folds >= 2")
        if p.temporal_cutoff:
            import datetime as _dt
            _dt.date.fromisoformat(p.temporal_cutoff)
        if self.fit.bayes not in ("auto", "true", "false"):
            raise ValueError('fit.bayes must be "auto", "true" or "false"')
        if not set(self.fit.methods) <= {"hmm", "pssm", "pwm", "framework"} or not self.fit.methods:
            raise ValueError("fit.methods: hmm, pssm, pwm, framework")
        if set(self.fit.regions) - {"SIG", "PRO", "MAT", "POST"}:
            raise ValueError("fit.regions: SIG, PRO, MAT, POST")
        if self.anatomy.method not in ("transfer", "envelope"):
            raise ValueError('anatomy.method must be "transfer" or "envelope"')
        if self.anatomy.extension not in ("within", "cross"):
            raise ValueError('anatomy.extension must be "within" or "cross"')
        if not isinstance(self.anatomy.sites, bool) and (
                set(self.anatomy.sites) - {"signal", "mature_start", "mature_end"}):
            raise ValueError('anatomy.sites: true, false, or a list of "signal", "mature_start", "mature_end"')
        if self.evaluate.bootstrap < 0:
            raise ValueError("evaluate.bootstrap must be >= 0")
        if self.score.pssm not in ("auto", "true", "false"):
            raise ValueError('score.pssm must be "auto", "true" or "false"')
        if self.build.pssm not in ("auto", "true", "false"):
            raise ValueError('build.pssm must be "auto", "true" or "false"')
        if set(self.build.regions) - {"SIG", "PRO", "MAT", "POST"}:
            raise ValueError("build.regions: only SIG, PRO, MAT, POST")
        bad = set(self.audit.conflict_kinds) - {"FULL", "SIG", "PRO", "MAT", "POST"}
        if bad:
            raise ValueError(f"audit.conflict_kinds: unknown kind(s) {sorted(bad)}")


_SECTIONS = {
    "source": SourceConfig,
    "labels": LabelsConfig,
    "regions": RegionsConfig,
    "cluster": ClusterConfig,
    "tiers": TiersConfig,
    "audit": AuditConfig,
    "partition": PartitionConfig,
    "negatives": NegativesConfig,
    "build": BuildConfig,
    "score": ScoreConfig,
    "fit": FitConfig,
    "evaluate": EvaluateConfig,
    "anatomy": AnatomyConfig,
    "package": PackageConfig,
    "output": OutputConfig,
}
_PATH_FIELDS = {("labels", "overrides"), ("labels", "entry_overrides"), ("output", "root")}


def load_config(path: Path | None) -> Config:
    """Load a TOML config; relative paths resolve against the config's folder."""
    if path is None:
        cfg = Config()
        cfg.validate()
        return cfg
    path = Path(path).resolve()
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    unknown = set(raw) - set(_SECTIONS)
    if unknown:
        raise ValueError(f"Unknown config section(s): {', '.join(sorted(unknown))}")
    kwargs = {}
    for name, cls in _SECTIONS.items():
        values = dict(raw.get(name, {}))
        allowed = set(cls.__dataclass_fields__)
        bad = set(values) - allowed
        if bad:
            raise ValueError(f"Unknown key(s) in [{name}]: {', '.join(sorted(bad))}")
        for key in list(values):
            if (name, key) in _PATH_FIELDS and values[key] is not None:
                p = Path(values[key])
                values[key] = p if p.is_absolute() else (path.parent / p).resolve()
        kwargs[name] = cls(**values)
    cfg = Config(**kwargs, path=path)
    cfg.validate()
    return cfg


EXAMPLE_CONFIG = f'''# conodb configuration: one file fully describes a database build.
# Run with:  conodb run -c conodb.toml

[source]
query = "{DEFAULT_QUERY}"
page_size = 500

[labels]
# Curated mapping of raw UniProt family text -> superfamily/family/scope.
overrides = "labels.tsv"
# Per-accession corrections (accession, superfamily, family, scope, note).
# entry_overrides = "entry_labels.tsv"
# Only these taxonomic families provide classifier classes (others: reference).
taxon_scope = ["Conidae"]
# true = stop when UniProt introduces a family no rule or override covers.
strict = false

[regions]
allow_uncertain = true    # keep non-EXACT boundaries, flagged exact=false
min_length = 1            # record everything; filter when building profiles

[cluster]
identities = [0.9, 0.7, 0.5]
group_identity = 0.7      # defines independent entry groups (splits, tiers)
threads = 0               # 0 = all cores

[tiers]
tier1_min_groups = 10     # full profiles + k-fold calibration
tier2_min_groups = 3      # signal-first profiles, strong shrinkage

[audit]
conflict_kinds = ["FULL", "SIG", "MAT", "PRO"]
conflict_min_length = 10

[partition]
test_fraction = 0.2       # share of tier-1 entry groups held out as test
folds = 5                 # grouped cross-validation folds (tier 1)
seed = "conodictor-v3"    # changing it draws a new benchmark
exclude_conflicts = true  # leave out entries under review in label_conflicts.tsv
temporal_cutoff = ""      # e.g. "2024-01-01" for a time-split test set

[negatives]
# Real non-conopeptides: non-venom Conus + small secreted metazoan proteins, toxins excluded.
query = "(taxonomy_id:6490 OR (taxonomy_id:33208 AND cc_scl_term:SL-0243 AND length:[30 TO 300])) AND reviewed:true NOT cc_tissue_specificity:venom NOT keyword:KW-0800"
enabled = true
include_venom_other = true  # CRISP, PLA2, actinoporin... from the main set (scope 'other')
shuffles = 1                # shuffled copies per precursor, cysteines kept in place
min_length = 30
max_length = 1000
max_identity = 0.9          # real negatives this identical to a conopeptide are removed
max_entries = 5000          # cap on real UniProt negatives (stable choice by accession hash)

[build]
regions = ["SIG", "PRO", "MAT", "POST"]
min_length = 5            # shorter regions (e.g. 'GKR') are not profiled
min_unique = 2            # distinct sequences needed for a profile
use_uncertain = true      # use regions with non-EXACT UniProt boundaries
pssm = "auto"             # "auto": PSSMs if pftools is installed; "true": required; "false": HMMs only
matrix = "blosum45.cmp"   # pftools score matrix (file name searched in pftools data dirs, or a path)
pfw_options = []
pfmake_options = []
full_models = true        # full-precursor HMMs with region labels (region prediction)
pwm = true                # lightmotif matrices from each region alignment
threads = 0

[score]
hmm_min_bits = -1000.0    # store HMM scores uncensored; the fit sets the floor
max_mode = true           # no HMMER acceleration filters (short profiles)
pssm = "auto"             # PSSM scores with pftools pfsearch when PSSMs were built
pfsearch = "pfsearch"
pwm = true                # lightmotif ungapped matrices (comparison with pftools)
framework = true          # cysteine framework of the predicted mature region
threads = 0

[fit]
methods = ["hmm", "pssm", "framework"]   # feature types used by the classifiers
regions = ["SIG", "PRO", "MAT"]           # POST profiles add nothing to classification
hmm_floor = 0.0             # HMM scores below this are set to it (0 = censored; -20 keeps weak evidence)
baseline_fpr = 0.01         # baseline HMM cutoffs: false-positive rate per sequence
balance_negative_sources = true  # shuffled and real negatives weigh the same in total
tau_fixed = 0.5             # shrinkage of the fast penalised fit when no posterior exists
tau_from_posterior = true   # per-feature shrinkage from the Bayesian posterior
ship_all_data = true        # shipped model refitted on train + test after evaluation
beta_prior_sd = 2.5
tau_prior_sd = 1.0
alpha_prior_sd = 2.5
bayes = "auto"              # hierarchical Stan fit if CmdStan is installed ("true" = required)
cmdstan = ""                # CmdStan path if not the cmdstanpy default
chains = 4
iter_warmup = 1000
iter_sampling = 1000
adapt_delta = 0.9
seed = 20261105
cores = 0                   # 0 = all cores: parallel chains, then threads within each chain
grainsize = 0               # 0 = automatic

[evaluate]
identity_bands = [0.5, 0.6, 0.7, 0.8, 0.9]  # positives by highest identity to a training conopeptide
calibration_bins = [0.5, 0.7, 0.9, 0.99]
ablations = true            # test-set results of every variant of the grid (penalised fits)
select = true               # grouped CV on training rows: variant comparison + recommended settings
bootstrap = 1000            # resamples of entry groups for intervals
select_taus = [0.25, 0.5, 1.0]  # shrinkage values for the tau sensitivity check
seed = 20261106

[anatomy]
method = "transfer"         # region prediction: "transfer" (full-precursor HMM) or "envelope"
sites = ["signal", "mature_start"]  # boundaries refined with lightmotif site matrices (true = all, false = none)
cterm_trim = false          # drop a trailing G[KR]{0,3} from a C-terminal mature peptide (amidation)
extension = "within"        # label transfer: extend unaligned ends within the region ("within") or across ("cross")
min_occupancy = 0.5
min_run = 3
terminal_overhang = 5       # residues past where the extension stops that keep its label
site_left = 5
site_right = 3
search_radius = 6
distance_penalty = 0.5
pseudocount = 0.5
short_peptide_length = 60

[package]
cutoff_fpr = 0.01           # calibrated cutoffs written into the shipped profiles
archive = true

[output]
root = "db"
'''

EXAMPLE_LABELS = """raw_superfamily\traw_family\tsuperfamily\tfamily\tscope\tnote
# Rows override the built-in rules. raw_family '*' matches any family text.
# scope: conopeptide (used to build the classifier) | other (kept in audit) | exclude (dropped)
insulin family\t*\tCon-insulin\t\tconopeptide\tcone snail venom insulins
"""
