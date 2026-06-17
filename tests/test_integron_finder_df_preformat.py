import pandas as pd

from integron_finder_df_preformat import preformat

# preformat reshapes Integron_Finder .integrons element rows to the unified MGE
# schema. Each row is one element of an integron; the integron span is
# min(pos_beg)..max(pos_end) (1-based), the per-integron type is complete/CALIN/In0,
# and the id is "{iso}|{n}|{type}" with n parsed from ID_integron ("integron_01"->1).


def integrons_table(rows, extra=None):
    """rows: [(ID_replicon, ID_integron, pos_beg, pos_end, type)]; extra columns optional."""
    df = pd.DataFrame(
        {
            "ID_replicon": [r[0] for r in rows],
            "ID_integron": [r[1] for r in rows],
            "pos_beg": [r[2] for r in rows],
            "pos_end": [r[3] for r in rows],
            "type": [r[4] for r in rows],
        }
    )
    if extra:
        for k, v in extra.items():
            df[k] = v
    return df


def test_reshapes_to_unified_schema():
    # two elements of one integron collapse to a single row
    df = integrons_table(
        [
            ("iso1", "integron_01", 100, 150, "complete"),
            ("iso1", "integron_01", 900, 980, "complete"),
        ]
    )
    out = preformat([df])
    assert out.index.name == "id"
    assert list(out.columns) == ["iso", "beg", "end", "type"]
    assert len(out) == 1


def test_span_is_min_beg_max_end():
    # span covers every element of the integron
    df = integrons_table(
        [
            ("iso1", "integron_01", 900, 980, "complete"),
            ("iso1", "integron_01", 100, 150, "complete"),
            ("iso1", "integron_01", 400, 460, "complete"),
        ]
    )
    out = preformat([df]).reset_index()
    assert (out["beg"].iloc[0], out["end"].iloc[0]) == (100, 980)


def test_id_is_iso_n_type():
    # n is the integer parsed from ID_integron; type is appended
    df = integrons_table([("NZ_AP022362.1", "integron_01", 209091, 212155, "complete")])
    out = preformat([df])
    assert out.index.tolist() == ["NZ_AP022362.1|1|complete"]
    assert out["type"].iloc[0] == "complete"


def test_type_carried_through():
    df = integrons_table([("iso1", "integron_01", 10, 20, "CALIN")])
    out = preformat([df])
    assert out["type"].iloc[0] == "CALIN"


def test_multiple_integrons_per_isolate_stay_unique():
    df = integrons_table(
        [
            ("iso1", "integron_01", 100, 200, "complete"),
            ("iso1", "integron_02", 5000, 5200, "CALIN"),
        ]
    )
    out = preformat([df])
    assert out.index.tolist() == ["iso1|1|complete", "iso1|2|CALIN"]
    assert out.index.is_unique


def test_concatenates_multiple_isolates():
    a = integrons_table([("isoA", "integron_01", 1, 2, "In0")])
    b = integrons_table([("isoB", "integron_01", 3, 4, "complete")])
    out = preformat([a, b])
    assert len(out) == 2
    assert set(out["iso"]) == {"isoA", "isoB"}


def test_extra_input_columns_dropped():
    # the real table also carries element/strand/evalue/type_elt/annotation/... ;
    # only the unified-schema columns should survive
    df = integrons_table(
        [("iso1", "integron_01", 10, 20, "complete")],
        extra={"type_elt": "attC", "strand": -1, "evalue": 1e-9},
    )
    out = preformat([df])
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_rows_without_integron_id_dropped():
    # a 'no integron found' placeholder row (NA ID_integron) contributes nothing
    df = integrons_table([("iso1", None, None, None, None)])
    out = preformat([df])
    assert out.empty
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_mixed_type_within_one_integron_raises():
    import pytest

    df = integrons_table(
        [
            ("iso1", "integron_01", 10, 20, "complete"),
            ("iso1", "integron_01", 30, 40, "CALIN"),
        ]
    )
    with pytest.raises(AssertionError):
        preformat([df])
