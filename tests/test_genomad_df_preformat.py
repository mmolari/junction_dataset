import pandas as pd
import pytest

from genomad_df_preformat import preformat

# preformat reshapes geNomad virus_summary rows to the unified MGE schema.
# A provirus row is named "{iso}|provirus_{beg}_{end}" with a "{beg}-{end}"
# coordinates field; iso is the seq_name prefix and beg/end come from coordinates.


def virus_summary(rows, extra=None):
    """rows: [(seq_name, coordinates)]; extra: optional dict of additional columns."""
    df = pd.DataFrame(
        {"seq_name": [r[0] for r in rows], "coordinates": [r[1] for r in rows]}
    )
    if extra:
        for k, v in extra.items():
            df[k] = v
    return df


def test_reshapes_to_unified_schema():
    df = virus_summary([("iso1|provirus_100_200", "100-200")])
    out = preformat([df])
    assert out.index.name == "id"
    assert list(out.columns) == ["iso", "beg", "end", "type"]
    assert out.index.tolist() == ["iso1|provirus_100_200"]


def test_iso_extracted_from_seq_name_prefix():
    # iso is the seq_name up to the first "|"
    df = virus_summary([("NZ_CP014495.1|provirus_5_9", "5-9")])
    out = preformat([df])
    assert out["iso"].iloc[0] == "NZ_CP014495.1"


def test_beg_end_parsed_as_ints():
    df = virus_summary([("iso1|provirus_100_200", "100-200")])
    out = preformat([df])
    assert (out["beg"].iloc[0], out["end"].iloc[0]) == (100, 200)
    assert pd.api.types.is_integer_dtype(out["beg"])
    assert pd.api.types.is_integer_dtype(out["end"])


def test_type_is_constant_prophage():
    df = virus_summary([("iso1|provirus_1_2", "1-2"), ("iso1|provirus_3_4", "3-4")])
    out = preformat([df])
    assert (out["type"] == "prophage").all()


def test_extra_input_columns_dropped():
    # geNomad's table carries many columns (topology, taxonomy, ...); only the
    # unified-schema columns should survive
    df = virus_summary(
        [("iso1|provirus_1_2", "1-2")],
        extra={"topology": "Provirus", "taxonomy": "Viruses;...", "virus_score": 0.99},
    )
    out = preformat([df])
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_concatenates_multiple_tables():
    # one table per isolate, concatenated into a single frame
    a = virus_summary([("isoA|provirus_1_2", "1-2")])
    b = virus_summary([("isoB|provirus_3_4", "3-4"), ("isoB|provirus_5_6", "5-6")])
    out = preformat([a, b])
    assert len(out) == 3
    assert set(out["iso"]) == {"isoA", "isoB"}


def test_duplicate_id_rejected():
    # the id index (seq_name) must be unique; verify_integrity guards it
    df = virus_summary([("iso1|provirus_1_2", "1-2"), ("iso1|provirus_1_2", "1-2")])
    with pytest.raises(ValueError):
        preformat([df])
