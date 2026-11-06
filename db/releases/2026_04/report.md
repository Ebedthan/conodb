# ConoDictor database audit: UniProt 2026_04

- UniProt release: **2026_04** (04-November-2026)
- Query: `taxonomy_id:37797 AND (cc_tissue_specificity:venom OR cc_scl_term:SL-0177) AND reviewed:true`
- Entries downloaded: **1493**
- Conopeptide superfamilies: **31** (tier 1: 7, tier 2: 11, tier 3: 13)
- Tiers count entry groups at 0.7 identity: tier 1 >= 10, tier 2 >= 3, tier 3 below.

## Entries by taxonomic family

| taxon_family | entries | by_scope |
|---|---|---|
| Conidae | 1394 | conopeptide 1179, unlabelled 164, unmapped 51 |
| Turridae | 53 | unlabelled 30, reference 13, unmapped 10 |
| Terebridae | 46 | unlabelled 26, unmapped 20 |

## Conopeptide superfamilies

`groups_x` = independent entry groups if full sequences are grouped at identity x.

| superfamily | tier | entries | entry_groups | species | complete_precursors | fragments | groups_0.9 | groups_0.7 | groups_0.5 |
|---|---|---|---|---|---|---|---|---|---|
| O1 | 1 | 285 | 93 | 50 | 208 | 51 | 197 | 93 | 21 |
| A | 1 | 191 | 65 | 53 | 102 | 43 | 138 | 65 | 21 |
| M | 1 | 183 | 63 | 41 | 95 | 21 | 141 | 63 | 27 |
| T | 1 | 147 | 44 | 29 | 115 | 7 | 100 | 44 | 18 |
| O2 | 1 | 105 | 30 | 36 | 80 | 5 | 80 | 30 | 12 |
| I1 | 1 | 38 | 16 | 14 | 11 | 10 | 24 | 16 | 10 |
| I2 | 1 | 32 | 14 | 18 | 29 | 2 | 22 | 14 | 6 |
| Con-insulin | 2 | 25 | 8 | 10 | 25 | 0 | 13 | 8 | 3 |
| P | 2 | 8 | 7 | 8 | 6 | 0 | 8 | 7 | 4 |
| C | 2 | 18 | 6 | 13 | 18 | 0 | 13 | 6 | 1 |
| D | 2 | 32 | 6 | 11 | 31 | 0 | 13 | 6 | 3 |
| S | 2 | 9 | 6 | 6 | 6 | 3 | 6 | 6 | 1 |
| B | 2 | 14 | 5 | 10 | 7 | 3 | 11 | 5 | 4 |
| I3 | 2 | 8 | 5 | 4 | 8 | 0 | 7 | 5 | 2 |
| L | 2 | 10 | 5 | 4 | 8 | 0 | 7 | 5 | 5 |
| O3 | 2 | 22 | 4 | 6 | 21 | 1 | 14 | 4 | 1 |
| H | 2 | 4 | 3 | 2 | 4 | 0 | 4 | 3 | 2 |
| J | 2 | 8 | 3 | 4 | 7 | 1 | 6 | 3 | 1 |
| E | 3 | 3 | 2 | 3 | 3 | 0 | 3 | 2 | 2 |
| F | 3 | 2 | 2 | 2 | 2 | 0 | 2 | 2 | 1 |
| N | 3 | 4 | 2 | 2 | 4 | 0 | 3 | 2 | 2 |
| Q | 3 | 2 | 2 | 2 | 2 | 0 | 2 | 2 | 1 |
| U | 3 | 5 | 2 | 4 | 4 | 0 | 3 | 2 | 1 |
| V | 3 | 2 | 2 | 2 | 2 | 0 | 2 | 2 | 1 |
| B2 | 3 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |
| G | 3 | 2 | 1 | 1 | 1 | 0 | 2 | 1 | 1 |
| G2 | 3 | 2 | 1 | 1 | 2 | 0 | 1 | 1 | 1 |
| K | 3 | 5 | 1 | 1 | 5 | 0 | 2 | 1 | 1 |
| NSf-1 | 3 | 2 | 1 | 1 | 2 | 0 | 2 | 1 | 1 |
| R | 3 | 9 | 1 | 2 | 7 | 0 | 4 | 1 | 1 |
| Y | 3 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |

### Tiers under each grouping threshold

| identity | tier 1 | tier 2 | tier 3 |
|---|---|---|---|
| 0.9 | 12 | 10 | 9 |
| 0.7 (current) | 7 | 11 | 13 |
| 0.5 | 6 | 6 | 19 |

## Families within superfamilies

| superfamily | family | tier | entries | entry_groups | species |
|---|---|---|---|---|---|
| C | Consomatin | 2 | 15 | 4 | 12 |
| G2 | 1 | 3 | 2 | 1 | 1 |
| M | Conopeptide Y | 3 | 2 | 1 | 2 |
| O2 | Contryphan | 1 | 45 | 10 | 29 |

## Other entries (not classifier classes)

`reference`: conopeptide families outside `labels.taxon_scope`; `other`: venom proteins kept as future negatives; `unmapped`: awaiting curation.

| superfamily | scope | entries | entry_groups |
|---|---|---|---|
| (unlabelled) | unlabelled | 220 | 157 |
| CRISP | unmapped | 2 | 2 |
| ConoGAY | unmapped | 1 | 1 |
| FARP (FMRFamide related peptide) | unmapped | 9 | 4 |
| I2-like | reference | 2 | 1 |
| NPY | unmapped | 3 | 1 |
| P-like | reference | 11 | 5 |
| Pg turripeptide superfamily | unmapped | 8 | 1 |
| actinoporin | unmapped | 3 | 3 |
| conophan | unmapped | 4 | 3 |
| elevenin | unmapped | 3 | 1 |
| glycoprotein hormones subunit alpha | unmapped | 1 | 1 |
| glycoprotein hormones subunit beta | unmapped | 1 | 1 |
| glycosyl hydrolase 56 | unmapped | 3 | 1 |
| phospholipase A2 | unmapped | 4 | 2 |
| teretoxin A (TA) superfamily | unmapped | 2 | 1 |
| teretoxin C (TC) superfamily | unmapped | 4 | 4 |
| teretoxin H (TH) superfamily | unmapped | 2 | 2 |
| teretoxin J (TJ) superfamily | unmapped | 1 | 1 |
| teretoxin M (TM) superfamily | unmapped | 6 | 6 |
| teretoxin N (TN) superfamily | unmapped | 1 | 1 |
| turripeptide | unmapped | 2 | 2 |
| vasopressin/oxytocin | unmapped | 15 | 9 |
| venom Kunitz-type | unmapped | 6 | 3 |

## Entries without a UniProt family

220 entries: 81 complete precursors, 82 with a signal peptide, 112 mature peptide only, 35 flagged as fragments.

Most frequent first word of the protein name: Conotoxin (115), Turripeptide (28), Augerpeptide (15), Mu-conotoxin-like (14), Teretoxin (7), Alpha-conotoxin (2), Conopeptide (2), Omega-conotoxin (2), A-conotoxin (1), Augerpeptide-s11a (1), Augerpeptide-s6a (1), Augerpeptide-s7a (1)

## Labels needing curation

Raw UniProt family text that no rule or override covers. Add rows to the labels file with scope `conopeptide`, `other` or `exclude`.

| raw_superfamily | raw_family | entries |
|---|---|---|
| vasopressin/oxytocin family |  | 15 |
| FARP (FMRFamide related peptide) family |  | 9 |
| Pg turripeptide superfamily |  | 8 |
| venom Kunitz-type family |  | 6 |
| teretoxin M (TM) superfamily |  | 6 |
| phospholipase A2 family | Group IX subfamily | 4 |
| conophan family |  | 4 |
| teretoxin C (TC) superfamily |  | 4 |
| glycosyl hydrolase 56 family |  | 3 |
| NPY family |  | 3 |
| elevenin family |  | 3 |
| actinoporin family | Conoidea subfamily | 3 |
| CRISP family |  | 2 |
| turripeptide family |  | 2 |
| teretoxin H (TH) superfamily |  | 2 |
| teretoxin A (TA) superfamily |  | 2 |
| ConoGAY family |  | 1 |
| glycoprotein hormones subunit alpha family |  | 1 |
| teretoxin J (TJ) superfamily |  | 1 |
| glycoprotein hormones subunit beta family |  | 1 |
| teretoxin N (TN) superfamily |  | 1 |

## Label conflicts (15)

Sequences >= 0.9 identical (regions >= 10 aa, kinds: FULL, SIG, MAT, PRO) but labelled with different superfamilies, most severe first. FULL = whole sequences; SIG = signal peptides (superfamilies are defined by signal similarity). `review` lists the minority labels. By severity: FULL 4, SIG 2, MAT 6, PRO 3.

| severity | kind | labels | review |
|---|---|---|---|
| FULL | FULL, SIG, MAT, PRO | M x1, A x1 | B2KPN7:M, P0CH39:A |
| FULL | FULL | O1 x4, M x2 | P0DKQ4:M, P0DKQ5:M |
| FULL | FULL | M x1, O1 x1 | P0DKQ3:M, Q9XZK2:O1 |
| FULL | FULL, MAT | L x1, R x1 | P0C6S2:L, P84704:R |
| SIG | SIG | O1 x34, M x2 | P0DKQ3:M, P0DKQ5:M |
| SIG | SIG | T x2, A x1 | P0DPL8:A |
| MAT | MAT | R x2, L x1 | P0C6S3:L |
| MAT | MAT | O1 x2, M x1 | P0DKQ3:M |
| MAT | MAT | O1 x1, M x1 | P0C8V8:O1, P0DKQ4:M |
| MAT | MAT | O1 x3, M x1 | P0DKQ5:M |
| MAT | MAT | T x2, M x1 | P0DMB3:M |
| MAT | MAT | T x3, A x1 | P0DPL8:A |
| PRO | PRO | T x1, A x1 | A0A3G3C7S2:T, P0DPL8:A |
| PRO | PRO | N x1, O2 x1 | A0A125S9F0:N, C8CK74:O2 |
| PRO | PRO | O1 x4, M x2 | P0DKQ3:M, P0DKQ5:M |

## Parsing warnings

| warning | count |
|---|---|
| Propeptide unknown boundary | 1 |
| Signal unknown boundary | 1 |

