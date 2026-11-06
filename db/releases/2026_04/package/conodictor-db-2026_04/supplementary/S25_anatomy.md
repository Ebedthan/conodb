# Region prediction benchmark: UniProt 2026_04

185 held-out conopeptides; profiles, full-precursor models and cleavage-site matrices built from training entries only. Coordinates are scored only where the UniProt boundary is EXACT. These are test-set numbers: the defaults are chosen by cross-validation on training entries (`selection.md`), not from this table.

- full-precursor models: match states = columns occupied >= 0.5; labels = majority region
- transfer/within: unaligned ends extended only within the region of the last aligned state; transfer/cross: across regions; ohN: N residues past where the extension stops keep its label
- sites (S = signal end, MS = mature start, ME = mature end): lightmotif matrices of 5 residues before / 3 after each cleavage; search +/- 6 residues, penalty 0.5 per residue shifted; Cterm = a C-terminal mature peptide loses a trailing G[KR]{0,3} (amidation)
- current settings (framework feature, conodictor): **transfer/within/oh5 + sites S+MS + Cterm**; recommended by cross-validation: **transfer/within/oh5 + sites S+MS + Cterm**

## Family model: true family

### Boundaries and mature peptide

| variant | n | boundary_exact | signal_end_exact | mature_start_exact | mature_end_exact | mature_exact | residue_accuracy |
|---|---|---|---|---|---|---|---|
| envelope | 184 | 0.3811 | 0.6077 | 0.1397 | 0.4581 | 0.1173 | 0.786 |
| envelope + sites S+MS+ME | 184 | 0.4283 | 0.7 | 0.2011 | 0.4581 | 0.1564 | 0.7767 |
| transfer/within/oh0 | 184 | 0.5061 | 0.6231 | 0.3966 | 0.5307 | 0.2346 | 0.9091 |
| transfer/within/oh0 + sites S | 184 | 0.5389 | 0.7308 | 0.4078 | 0.5307 | 0.2402 | 0.9111 |
| transfer/within/oh0 + sites MS | 184 | 0.5881 | 0.6308 | 0.6145 | 0.5307 | 0.3296 | 0.9006 |
| transfer/within/oh0 + sites S+MS | 184 | 0.6127 | 0.7231 | 0.6145 | 0.5307 | 0.3296 | 0.9022 |
| transfer/within/oh0 + sites S+MS+ME | 184 | 0.5656 | 0.7231 | 0.6145 | 0.4022 | 0.2235 | 0.883 |
| transfer/within/oh5 | 184 | 0.5471 | 0.6231 | 0.419 | 0.6201 | 0.3296 | 0.9263 |
| transfer/within/oh5 + sites S | 184 | 0.5799 | 0.7308 | 0.4302 | 0.6201 | 0.3352 | 0.9282 |
| transfer/within/oh5 + sites MS | 184 | 0.6291 | 0.6308 | 0.6369 | 0.6201 | 0.4469 | 0.9217 |
| transfer/within/oh5 + sites S+MS | 184 | 0.6537 | 0.7231 | 0.6369 | 0.6201 | 0.4469 | 0.9233 |
| transfer/within/oh5 + sites S+MS+ME | 184 | 0.6168 | 0.7231 | 0.6369 | 0.5196 | 0.3575 | 0.9144 |
| transfer/cross/oh0 | 184 | 0.5123 | 0.6231 | 0.3966 | 0.5475 | 0.2346 | 0.9128 |
| transfer/cross/oh0 + sites S | 184 | 0.5451 | 0.7308 | 0.4078 | 0.5475 | 0.2402 | 0.9148 |
| transfer/cross/oh0 + sites MS | 184 | 0.5943 | 0.6308 | 0.6145 | 0.5475 | 0.3296 | 0.9039 |
| transfer/cross/oh0 + sites S+MS | 184 | 0.6189 | 0.7231 | 0.6145 | 0.5475 | 0.3296 | 0.9055 |
| transfer/cross/oh0 + sites S+MS+ME | 184 | 0.5717 | 0.7231 | 0.6145 | 0.419 | 0.2235 | 0.8863 |
| transfer/cross/oh5 | 184 | 0.5451 | 0.6231 | 0.3966 | 0.6369 | 0.2905 | 0.9217 |
| transfer/cross/oh5 + sites S | 184 | 0.5779 | 0.7308 | 0.4078 | 0.6369 | 0.2961 | 0.9236 |
| transfer/cross/oh5 + sites MS | 184 | 0.627 | 0.6308 | 0.6145 | 0.6369 | 0.4078 | 0.9127 |
| transfer/cross/oh5 + sites S+MS | 184 | 0.6516 | 0.7231 | 0.6145 | 0.6369 | 0.4078 | 0.9143 |
| transfer/cross/oh5 + sites S+MS+ME | 184 | 0.6066 | 0.7231 | 0.6145 | 0.514 | 0.3073 | 0.8989 |
| envelope + Cterm | 184 | 0.4877 | 0.6077 | 0.1397 | 0.7486 | 0.1341 | 0.7941 |
| envelope + sites S+MS+ME + Cterm | 184 | 0.5328 | 0.7 | 0.2011 | 0.743 | 0.1788 | 0.7849 |
| transfer/within/oh0 + Cterm | 184 | 0.543 | 0.6231 | 0.3966 | 0.6313 | 0.257 | 0.9103 |
| transfer/within/oh0 + sites S + Cterm | 184 | 0.5758 | 0.7308 | 0.4078 | 0.6313 | 0.2626 | 0.9122 |
| transfer/within/oh0 + sites MS + Cterm | 184 | 0.625 | 0.6308 | 0.6145 | 0.6313 | 0.3743 | 0.9018 |
| transfer/within/oh0 + sites S+MS + Cterm | 184 | 0.6496 | 0.7231 | 0.6145 | 0.6313 | 0.3743 | 0.9034 |
| transfer/within/oh0 + sites S+MS+ME + Cterm | 184 | 0.5861 | 0.7231 | 0.6145 | 0.4581 | 0.2737 | 0.8841 |
| transfer/within/oh5 + Cterm | 184 | 0.623 | 0.6231 | 0.419 | 0.8268 | 0.3743 | 0.9319 |
| transfer/within/oh5 + sites S + Cterm | 184 | 0.6557 | 0.7308 | 0.4302 | 0.8268 | 0.3799 | 0.9339 |
| transfer/within/oh5 + sites MS + Cterm | 184 | 0.7049 | 0.6308 | 0.6369 | 0.8268 | 0.5419 | 0.9274 |
| transfer/within/oh5 + sites S+MS + Cterm | 184 | 0.7295 | 0.7231 | 0.6369 | 0.8268 | 0.5419 | 0.929 |
| transfer/within/oh5 + sites S+MS+ME + Cterm | 184 | 0.6906 | 0.7231 | 0.6369 | 0.7207 | 0.4637 | 0.9199 |
| transfer/cross/oh0 + Cterm | 184 | 0.5492 | 0.6231 | 0.3966 | 0.648 | 0.257 | 0.914 |
| transfer/cross/oh0 + sites S + Cterm | 184 | 0.582 | 0.7308 | 0.4078 | 0.648 | 0.2626 | 0.9159 |
| transfer/cross/oh0 + sites MS + Cterm | 184 | 0.6311 | 0.6308 | 0.6145 | 0.648 | 0.3743 | 0.905 |
| transfer/cross/oh0 + sites S+MS + Cterm | 184 | 0.6557 | 0.7231 | 0.6145 | 0.648 | 0.3743 | 0.9066 |
| transfer/cross/oh0 + sites S+MS+ME + Cterm | 184 | 0.5922 | 0.7231 | 0.6145 | 0.4749 | 0.2737 | 0.8874 |
| transfer/cross/oh5 + Cterm | 184 | 0.6004 | 0.6231 | 0.3966 | 0.7877 | 0.324 | 0.9259 |
| transfer/cross/oh5 + sites S + Cterm | 184 | 0.6332 | 0.7308 | 0.4078 | 0.7877 | 0.3296 | 0.9278 |
| transfer/cross/oh5 + sites MS + Cterm | 184 | 0.6824 | 0.6308 | 0.6145 | 0.7877 | 0.486 | 0.9169 |
| transfer/cross/oh5 + sites S+MS + Cterm | 184 | 0.707 | 0.7231 | 0.6145 | 0.7877 | 0.486 | 0.9185 |
| transfer/cross/oh5 + sites S+MS+ME + Cterm | 184 | 0.6619 | 0.7231 | 0.6145 | 0.6648 | 0.4022 | 0.9031 |

### Exact mature peptide by sequence type

| variant | complete precursor | no signal peptide | signal, incomplete |
|---|---|---|---|
| envelope | 0.0 | 0.4667 | 0.0 |
| envelope + sites S+MS+ME | 0.0538 | 0.4667 | 0.0 |
| transfer/within/oh0 | 0.1462 | 0.4667 | 0.5 |
| transfer/within/oh0 + sites S | 0.1538 | 0.4667 | 0.5 |
| transfer/within/oh0 + sites MS | 0.2846 | 0.4667 | 0.25 |
| transfer/within/oh0 + sites S+MS | 0.2846 | 0.4667 | 0.25 |
| transfer/within/oh0 + sites S+MS+ME | 0.1692 | 0.3778 | 0.25 |
| transfer/within/oh5 | 0.2308 | 0.6 | 0.5 |
| transfer/within/oh5 + sites S | 0.2385 | 0.6 | 0.5 |
| transfer/within/oh5 + sites MS | 0.3923 | 0.6 | 0.5 |
| transfer/within/oh5 + sites S+MS | 0.3923 | 0.6 | 0.5 |
| transfer/within/oh5 + sites S+MS+ME | 0.3 | 0.5111 | 0.5 |
| transfer/cross/oh0 | 0.1462 | 0.4667 | 0.5 |
| transfer/cross/oh0 + sites S | 0.1538 | 0.4667 | 0.5 |
| transfer/cross/oh0 + sites MS | 0.2846 | 0.4667 | 0.25 |
| transfer/cross/oh0 + sites S+MS | 0.2846 | 0.4667 | 0.25 |
| transfer/cross/oh0 + sites S+MS+ME | 0.1692 | 0.3778 | 0.25 |
| transfer/cross/oh5 | 0.2231 | 0.4667 | 0.5 |
| transfer/cross/oh5 + sites S | 0.2308 | 0.4667 | 0.5 |
| transfer/cross/oh5 + sites MS | 0.3846 | 0.4667 | 0.5 |
| transfer/cross/oh5 + sites S+MS | 0.3846 | 0.4667 | 0.5 |
| transfer/cross/oh5 + sites S+MS+ME | 0.2769 | 0.3778 | 0.5 |
| envelope + Cterm | 0.0154 | 0.4889 | 0.0 |
| envelope + sites S+MS+ME + Cterm | 0.0769 | 0.4889 | 0.0 |
| transfer/within/oh0 + Cterm | 0.1615 | 0.5111 | 0.5 |
| transfer/within/oh0 + sites S + Cterm | 0.1692 | 0.5111 | 0.5 |
| transfer/within/oh0 + sites MS + Cterm | 0.3385 | 0.4667 | 0.5 |
| transfer/within/oh0 + sites S+MS + Cterm | 0.3385 | 0.4667 | 0.5 |
| transfer/within/oh0 + sites S+MS+ME + Cterm | 0.2231 | 0.4 | 0.5 |
| transfer/within/oh5 + Cterm | 0.2692 | 0.6667 | 0.5 |
| transfer/within/oh5 + sites S + Cterm | 0.2769 | 0.6667 | 0.5 |
| transfer/within/oh5 + sites MS + Cterm | 0.5154 | 0.6 | 0.75 |
| transfer/within/oh5 + sites S+MS + Cterm | 0.5154 | 0.6 | 0.75 |
| transfer/within/oh5 + sites S+MS+ME + Cterm | 0.4308 | 0.5333 | 0.75 |
| transfer/cross/oh0 + Cterm | 0.1615 | 0.5111 | 0.5 |
| transfer/cross/oh0 + sites S + Cterm | 0.1692 | 0.5111 | 0.5 |
| transfer/cross/oh0 + sites MS + Cterm | 0.3385 | 0.4667 | 0.5 |
| transfer/cross/oh0 + sites S+MS + Cterm | 0.3385 | 0.4667 | 0.5 |
| transfer/cross/oh0 + sites S+MS+ME + Cterm | 0.2231 | 0.4 | 0.5 |
| transfer/cross/oh5 + Cterm | 0.2462 | 0.5333 | 0.5 |
| transfer/cross/oh5 + sites S + Cterm | 0.2538 | 0.5333 | 0.5 |
| transfer/cross/oh5 + sites MS + Cterm | 0.4846 | 0.4667 | 0.75 |
| transfer/cross/oh5 + sites S+MS + Cterm | 0.4846 | 0.4667 | 0.75 |
| transfer/cross/oh5 + sites S+MS+ME + Cterm | 0.3923 | 0.4 | 0.75 |

### Per region (current and recommended settings)

| variant | region | n_true | detected | false_region | n_coord | start_exact | start_within2 | start_mae | end_exact | end_within2 | end_mae |
|---|---|---|---|---|---|---|---|---|---|---|---|
| transfer/within/oh5 + sites S+MS + Cterm | SIG | 134 | 1.0 | 0.0 | 130 | 1.0 | 1.0 | 0.0 | 0.7231 | 0.7923 | 1.023 |
| transfer/within/oh5 + sites S+MS + Cterm | PRO | 140 | 0.9857 | 0.2727 | 123 | 0.7236 | 0.7886 | 0.902 | 0.5122 | 0.8943 | 1.024 |
| transfer/within/oh5 + sites S+MS + Cterm | MAT | 184 | 0.9837 |  | 176 | 0.6477 | 0.8693 | 1.324 | 0.8409 | 0.9545 | 0.358 |
| transfer/within/oh5 + sites S+MS + Cterm | POST | 15 | 0.1333 | 0.0 | 2 | 0.5 | 0.5 | 5.0 | 0.5 | 0.5 | 2.0 |

## Family model: predicted family

### Boundaries and mature peptide

| variant | n | boundary_exact | signal_end_exact | mature_start_exact | mature_end_exact | mature_exact | residue_accuracy |
|---|---|---|---|---|---|---|---|
| envelope | 184 | 0.3811 | 0.6077 | 0.1397 | 0.4581 | 0.1173 | 0.7848 |
| envelope + sites S+MS+ME | 184 | 0.4283 | 0.7 | 0.2011 | 0.4581 | 0.1564 | 0.7747 |
| transfer/within/oh0 | 184 | 0.502 | 0.6231 | 0.3855 | 0.5307 | 0.2235 | 0.9078 |
| transfer/within/oh0 + sites S | 184 | 0.5348 | 0.7308 | 0.3966 | 0.5307 | 0.2291 | 0.9098 |
| transfer/within/oh0 + sites MS | 184 | 0.584 | 0.6308 | 0.6034 | 0.5307 | 0.3184 | 0.8941 |
| transfer/within/oh0 + sites S+MS | 184 | 0.6086 | 0.7231 | 0.6034 | 0.5307 | 0.3184 | 0.8957 |
| transfer/within/oh0 + sites S+MS+ME | 184 | 0.5615 | 0.7231 | 0.6034 | 0.4022 | 0.2123 | 0.8765 |
| transfer/within/oh5 | 184 | 0.5471 | 0.6231 | 0.419 | 0.6201 | 0.3296 | 0.9249 |
| transfer/within/oh5 + sites S | 184 | 0.5799 | 0.7308 | 0.4302 | 0.6201 | 0.3352 | 0.9269 |
| transfer/within/oh5 + sites MS | 184 | 0.6291 | 0.6308 | 0.6369 | 0.6201 | 0.4469 | 0.9215 |
| transfer/within/oh5 + sites S+MS | 184 | 0.6537 | 0.7231 | 0.6369 | 0.6201 | 0.4469 | 0.9231 |
| transfer/within/oh5 + sites S+MS+ME | 184 | 0.6168 | 0.7231 | 0.6369 | 0.5196 | 0.3575 | 0.9142 |
| transfer/cross/oh0 | 184 | 0.5082 | 0.6231 | 0.3855 | 0.5475 | 0.2235 | 0.9115 |
| transfer/cross/oh0 + sites S | 184 | 0.541 | 0.7308 | 0.3966 | 0.5475 | 0.2291 | 0.9135 |
| transfer/cross/oh0 + sites MS | 184 | 0.5902 | 0.6308 | 0.6034 | 0.5475 | 0.3184 | 0.8973 |
| transfer/cross/oh0 + sites S+MS | 184 | 0.6148 | 0.7231 | 0.6034 | 0.5475 | 0.3184 | 0.8989 |
| transfer/cross/oh0 + sites S+MS+ME | 184 | 0.5676 | 0.7231 | 0.6034 | 0.419 | 0.2123 | 0.8797 |
| transfer/cross/oh5 | 184 | 0.541 | 0.6231 | 0.3855 | 0.6369 | 0.2793 | 0.9204 |
| transfer/cross/oh5 + sites S | 184 | 0.5738 | 0.7308 | 0.3966 | 0.6369 | 0.2849 | 0.9223 |
| transfer/cross/oh5 + sites MS | 184 | 0.623 | 0.6308 | 0.6034 | 0.6369 | 0.3966 | 0.9062 |
| transfer/cross/oh5 + sites S+MS | 184 | 0.6475 | 0.7231 | 0.6034 | 0.6369 | 0.3966 | 0.9078 |
| transfer/cross/oh5 + sites S+MS+ME | 184 | 0.6025 | 0.7231 | 0.6034 | 0.514 | 0.2961 | 0.8923 |
| envelope + Cterm | 184 | 0.4877 | 0.6077 | 0.1397 | 0.7486 | 0.1341 | 0.793 |
| envelope + sites S+MS+ME + Cterm | 184 | 0.5328 | 0.7 | 0.2011 | 0.743 | 0.1788 | 0.783 |
| transfer/within/oh0 + Cterm | 184 | 0.5389 | 0.6231 | 0.3855 | 0.6313 | 0.2458 | 0.909 |
| transfer/within/oh0 + sites S + Cterm | 184 | 0.5717 | 0.7308 | 0.3966 | 0.6313 | 0.2514 | 0.9109 |
| transfer/within/oh0 + sites MS + Cterm | 184 | 0.6209 | 0.6308 | 0.6034 | 0.6313 | 0.3631 | 0.8952 |
| transfer/within/oh0 + sites S+MS + Cterm | 184 | 0.6455 | 0.7231 | 0.6034 | 0.6313 | 0.3631 | 0.8968 |
| transfer/within/oh0 + sites S+MS+ME + Cterm | 184 | 0.582 | 0.7231 | 0.6034 | 0.4581 | 0.2626 | 0.8776 |
| transfer/within/oh5 + Cterm | 184 | 0.623 | 0.6231 | 0.419 | 0.8268 | 0.3743 | 0.9306 |
| transfer/within/oh5 + sites S + Cterm | 184 | 0.6557 | 0.7308 | 0.4302 | 0.8268 | 0.3799 | 0.9325 |
| transfer/within/oh5 + sites MS + Cterm | 184 | 0.7049 | 0.6308 | 0.6369 | 0.8268 | 0.5419 | 0.9272 |
| transfer/within/oh5 + sites S+MS + Cterm | 184 | 0.7295 | 0.7231 | 0.6369 | 0.8268 | 0.5419 | 0.9288 |
| transfer/within/oh5 + sites S+MS+ME + Cterm | 184 | 0.6906 | 0.7231 | 0.6369 | 0.7207 | 0.4637 | 0.9197 |
| transfer/cross/oh0 + Cterm | 184 | 0.5451 | 0.6231 | 0.3855 | 0.648 | 0.2458 | 0.9127 |
| transfer/cross/oh0 + sites S + Cterm | 184 | 0.5779 | 0.7308 | 0.3966 | 0.648 | 0.2514 | 0.9146 |
| transfer/cross/oh0 + sites MS + Cterm | 184 | 0.627 | 0.6308 | 0.6034 | 0.648 | 0.3631 | 0.8985 |
| transfer/cross/oh0 + sites S+MS + Cterm | 184 | 0.6516 | 0.7231 | 0.6034 | 0.648 | 0.3631 | 0.9001 |
| transfer/cross/oh0 + sites S+MS+ME + Cterm | 184 | 0.5881 | 0.7231 | 0.6034 | 0.4749 | 0.2626 | 0.8809 |
| transfer/cross/oh5 + Cterm | 184 | 0.5963 | 0.6231 | 0.3855 | 0.7877 | 0.3128 | 0.9246 |
| transfer/cross/oh5 + sites S + Cterm | 184 | 0.6291 | 0.7308 | 0.3966 | 0.7877 | 0.3184 | 0.9265 |
| transfer/cross/oh5 + sites MS + Cterm | 184 | 0.6783 | 0.6308 | 0.6034 | 0.7877 | 0.4749 | 0.9104 |
| transfer/cross/oh5 + sites S+MS + Cterm | 184 | 0.7029 | 0.7231 | 0.6034 | 0.7877 | 0.4749 | 0.912 |
| transfer/cross/oh5 + sites S+MS+ME + Cterm | 184 | 0.6578 | 0.7231 | 0.6034 | 0.6648 | 0.3911 | 0.8966 |

### Exact mature peptide by sequence type

| variant | complete precursor | no signal peptide | signal, incomplete |
|---|---|---|---|
| envelope | 0.0 | 0.4667 | 0.0 |
| envelope + sites S+MS+ME | 0.0538 | 0.4667 | 0.0 |
| transfer/within/oh0 | 0.1462 | 0.4222 | 0.5 |
| transfer/within/oh0 + sites S | 0.1538 | 0.4222 | 0.5 |
| transfer/within/oh0 + sites MS | 0.2846 | 0.4222 | 0.25 |
| transfer/within/oh0 + sites S+MS | 0.2846 | 0.4222 | 0.25 |
| transfer/within/oh0 + sites S+MS+ME | 0.1692 | 0.3333 | 0.25 |
| transfer/within/oh5 | 0.2308 | 0.6 | 0.5 |
| transfer/within/oh5 + sites S | 0.2385 | 0.6 | 0.5 |
| transfer/within/oh5 + sites MS | 0.3923 | 0.6 | 0.5 |
| transfer/within/oh5 + sites S+MS | 0.3923 | 0.6 | 0.5 |
| transfer/within/oh5 + sites S+MS+ME | 0.3 | 0.5111 | 0.5 |
| transfer/cross/oh0 | 0.1462 | 0.4222 | 0.5 |
| transfer/cross/oh0 + sites S | 0.1538 | 0.4222 | 0.5 |
| transfer/cross/oh0 + sites MS | 0.2846 | 0.4222 | 0.25 |
| transfer/cross/oh0 + sites S+MS | 0.2846 | 0.4222 | 0.25 |
| transfer/cross/oh0 + sites S+MS+ME | 0.1692 | 0.3333 | 0.25 |
| transfer/cross/oh5 | 0.2231 | 0.4222 | 0.5 |
| transfer/cross/oh5 + sites S | 0.2308 | 0.4222 | 0.5 |
| transfer/cross/oh5 + sites MS | 0.3846 | 0.4222 | 0.5 |
| transfer/cross/oh5 + sites S+MS | 0.3846 | 0.4222 | 0.5 |
| transfer/cross/oh5 + sites S+MS+ME | 0.2769 | 0.3333 | 0.5 |
| envelope + Cterm | 0.0154 | 0.4889 | 0.0 |
| envelope + sites S+MS+ME + Cterm | 0.0769 | 0.4889 | 0.0 |
| transfer/within/oh0 + Cterm | 0.1615 | 0.4667 | 0.5 |
| transfer/within/oh0 + sites S + Cterm | 0.1692 | 0.4667 | 0.5 |
| transfer/within/oh0 + sites MS + Cterm | 0.3385 | 0.4222 | 0.5 |
| transfer/within/oh0 + sites S+MS + Cterm | 0.3385 | 0.4222 | 0.5 |
| transfer/within/oh0 + sites S+MS+ME + Cterm | 0.2231 | 0.3556 | 0.5 |
| transfer/within/oh5 + Cterm | 0.2692 | 0.6667 | 0.5 |
| transfer/within/oh5 + sites S + Cterm | 0.2769 | 0.6667 | 0.5 |
| transfer/within/oh5 + sites MS + Cterm | 0.5154 | 0.6 | 0.75 |
| transfer/within/oh5 + sites S+MS + Cterm | 0.5154 | 0.6 | 0.75 |
| transfer/within/oh5 + sites S+MS+ME + Cterm | 0.4308 | 0.5333 | 0.75 |
| transfer/cross/oh0 + Cterm | 0.1615 | 0.4667 | 0.5 |
| transfer/cross/oh0 + sites S + Cterm | 0.1692 | 0.4667 | 0.5 |
| transfer/cross/oh0 + sites MS + Cterm | 0.3385 | 0.4222 | 0.5 |
| transfer/cross/oh0 + sites S+MS + Cterm | 0.3385 | 0.4222 | 0.5 |
| transfer/cross/oh0 + sites S+MS+ME + Cterm | 0.2231 | 0.3556 | 0.5 |
| transfer/cross/oh5 + Cterm | 0.2462 | 0.4889 | 0.5 |
| transfer/cross/oh5 + sites S + Cterm | 0.2538 | 0.4889 | 0.5 |
| transfer/cross/oh5 + sites MS + Cterm | 0.4846 | 0.4222 | 0.75 |
| transfer/cross/oh5 + sites S+MS + Cterm | 0.4846 | 0.4222 | 0.75 |
| transfer/cross/oh5 + sites S+MS+ME + Cterm | 0.3923 | 0.3556 | 0.75 |

### Per region (current and recommended settings)

| variant | region | n_true | detected | false_region | n_coord | start_exact | start_within2 | start_mae | end_exact | end_within2 | end_mae |
|---|---|---|---|---|---|---|---|---|---|---|---|
| transfer/within/oh5 + sites S+MS + Cterm | SIG | 134 | 1.0 | 0.0 | 130 | 1.0 | 1.0 | 0.0 | 0.7231 | 0.7923 | 1.023 |
| transfer/within/oh5 + sites S+MS + Cterm | PRO | 140 | 0.9857 | 0.2955 | 123 | 0.7236 | 0.7886 | 0.902 | 0.5122 | 0.8943 | 1.024 |
| transfer/within/oh5 + sites S+MS + Cterm | MAT | 184 | 0.9837 |  | 176 | 0.6477 | 0.8693 | 1.33 | 0.8409 | 0.9545 | 0.358 |
| transfer/within/oh5 + sites S+MS + Cterm | POST | 15 | 0.1333 | 0.0 | 2 | 0.5 | 0.5 | 5.0 | 0.5 | 0.5 | 2.0 |

## Boundary errors of the current settings (true family; top 15 of 75)

`residues` lie between the predicted and the true boundary; `extra` = kept in the predicted mature peptide although outside the true one.

| boundary | direction | residues | n |
|---|---|---|---|
| signal_end | shorter | QLDA | 11 |
| mature_start | missing | R | 9 |
| signal_end | shorter | VVS | 5 |
| mature_end | missing | G | 4 |
| mature_start | extra | K | 4 |
| mature_start | extra | R | 4 |
| mature_end | missing | DWRS | 3 |
| mature_end | missing | S | 3 |
| mature_start | missing | 9 residues | 3 |
| mature_start | missing | S | 3 |
| signal_end | shorter | APC | 3 |
| signal_end | shorter | LP | 3 |
| mature_end | extra | LT | 2 |
| mature_end | extra | R | 2 |
| mature_end | missing | PAS | 2 |

All rows, including per-family results and within-1/5 rates: `anatomy_metrics.tsv`, `anatomy_sequences.tsv`; boundary errors of every variant: `anatomy_errors.tsv`; every prediction: `anatomy_predictions.tsv`.
