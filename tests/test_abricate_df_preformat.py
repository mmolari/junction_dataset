import pandas as pd

from abricate_df_preformat import preformat

# preformat reshapes abricate result rows to the unified MGE schema.
# abricate reports one row per resistance-gene hit with SEQUENCE (the isolate),
# 1-based START/END, and GENE; the id is "{iso}|{GENE}|{global row index}" and the
# type is the constant "AMR".


def abricate_table(rows, extra=None):
    """rows: [(SEQUENCE, START, END, GENE)]; extra: optional extra columns."""
    df = pd.DataFrame(
        {
            "SEQUENCE": [r[0] for r in rows],
            "START": [r[1] for r in rows],
            "END": [r[2] for r in rows],
            "GENE": [r[3] for r in rows],
        }
    )
    if extra:
        for k, v in extra.items():
            df[k] = v
    return df


def test_reshapes_to_unified_schema():
    df = abricate_table([("iso1", 100, 200, "blaCTX-M-15")])
    out = preformat([df])
    assert out.index.name == "id"
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_columns_mapped_from_abricate_fields():
    # iso<-SEQUENCE, beg<-START, end<-END
    df = abricate_table([("NZ_CP014495.1", 5, 9, "sul1")])
    out = preformat([df]).reset_index()
    assert out["iso"].iloc[0] == "NZ_CP014495.1"
    assert (out["beg"].iloc[0], out["end"].iloc[0]) == (5, 9)


def test_id_is_iso_gene_index():
    # composite id glues iso, gene and the row index
    df = abricate_table([("iso1", 1, 2, "sul1"), ("iso1", 3, 4, "tetA")])
    out = preformat([df])
    assert out.index.tolist() == ["iso1|sul1|0", "iso1|tetA|1"]


def test_type_is_constant_AMR():
    df = abricate_table([("iso1", 1, 2, "sul1"), ("iso1", 3, 4, "tetA")])
    out = preformat([df])
    assert (out["type"] == "AMR").all()


def test_extra_input_columns_dropped():
    # abricate's table carries many columns (%COVERAGE, %IDENTITY, PRODUCT, ...);
    # only the unified-schema columns should survive
    df = abricate_table(
        [("iso1", 1, 2, "sul1")],
        extra={"%COVERAGE": 100.0, "%IDENTITY": 99.5, "PRODUCT": "sulfonamide"},
    )
    out = preformat([df])
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_concatenates_multiple_tables():
    # one table per isolate, concatenated into a single frame; the global index
    # keeps composite ids unique even when (iso, gene) repeats
    a = abricate_table([("isoA", 1, 2, "sul1")])
    b = abricate_table([("isoB", 3, 4, "sul1"), ("isoB", 5, 6, "tetA")])
    out = preformat([a, b])
    assert len(out) == 3
    assert set(out["iso"]) == {"isoA", "isoB"}
    assert out.index.is_unique
