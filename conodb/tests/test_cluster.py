from conodb.cluster import cluster_levels, greedy_cluster

PRECURSOR = "MKLTCVLIVAVLLLTACQLITAEDSRGTQKHRALRSGCCSDPRCAWRC"
OTHER = "MSKLGVLLTICLLLFPLTALPMDGDQPADRPAERMQDDISSEQHPLFDKRQRCCGKGYSCYGTRC"


def test_fragment_joins_precursor_and_identicals_merge():
    items = [("p", PRECURSOR), ("frag", "GCCSDPRCAWRC"), ("same", PRECURSOR), ("o", OTHER)]
    a = greedy_cluster(items, 0.9)
    assert a["frag"] == a["p"] == a["same"]
    assert a["o"] != a["p"]


def test_levels_are_nested():
    items = [("p", PRECURSOR), ("q", PRECURSOR[:-4] + "AAAA"), ("o", OTHER)]
    levels = cluster_levels(items, [0.9, 0.5])
    hi, lo = levels[0.9], levels[0.5]
    for x in hi:
        for y in hi:
            if hi[x] == hi[y]:
                assert lo[x] == lo[y]
    assert len(set(lo.values())) <= len(set(hi.values()))


def test_deterministic_across_runs_and_threads():
    import random
    rng = random.Random(3)
    base = "MKLTCVLIVAVLLLTACQLITAEDSRGTQKHRALRSGCCSDPRCAWRC"
    items = [(f"s{i}", "".join(c if rng.random() > 0.25 else rng.choice("ADEKLMNQRSTV") for c in base)) for i in range(300)]
    ref = cluster_levels(items, [0.9, 0.7, 0.5], threads=1)
    for threads in (1, 2, 4):
        assert cluster_levels(items, [0.9, 0.7, 0.5], threads=threads) == ref
