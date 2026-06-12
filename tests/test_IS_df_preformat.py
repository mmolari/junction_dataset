import pandas as pd

from IS_df_preformat import preformat

# preformat reshapes ISEScan summary rows to the unified MGE schema.
# ISEScan reports one row per IS copy with seqID (the isolate), family, and
# 1-based isBegin/isEnd; the id is "{iso}|{family}|{global row index}".


def is_summary(rows, extra=None):
    """rows: [(seqID, family, isBegin, isEnd)]; extra: optional extra columns."""
    df = pd.DataFrame(
        {
            "seqID": [r[0] for r in rows],
            "family": [r[1] for r in rows],
            "isBegin": [r[2] for r in rows],
            "isEnd": [r[3] for r in rows],
        }
    )
    if extra:
        for k, v in extra.items():
            df[k] = v
    return df


def test_reshapes_to_unified_schema():
    df = is_summary([("iso1", "IS1", 100, 200)])
    out = preformat([df])
    assert out.index.name == "id"
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_columns_mapped_from_isescan_fields():
    # iso<-seqID, beg<-isBegin, end<-isEnd
    df = is_summary([("NZ_CP014495.1", "IS3", 5, 9)])
    out = preformat([df]).reset_index()
    assert out["iso"].iloc[0] == "NZ_CP014495.1"
    assert (out["beg"].iloc[0], out["end"].iloc[0]) == (5, 9)


def test_id_is_iso_family_index():
    # composite id glues iso, family and the row index
    df = is_summary([("iso1", "IS1", 1, 2), ("iso1", "IS6", 3, 4)])
    out = preformat([df])
    assert out.index.tolist() == ["iso1|IS1|0", "iso1|IS6|1"]


def test_type_is_constant_IS():
    df = is_summary([("iso1", "IS1", 1, 2), ("iso1", "IS6", 3, 4)])
    out = preformat([df])
    assert (out["type"] == "IS").all()


def test_extra_input_columns_dropped():
    # ISEScan's table carries many columns (cluster, type, ncopy4is, ...); only
    # the unified-schema columns should survive
    df = is_summary(
        [("iso1", "IS1", 1, 2)],
        extra={"cluster": "IS1_0", "ncopy4is": 3, "score": 0.99},
    )
    out = preformat([df])
    assert list(out.columns) == ["iso", "beg", "end", "type"]


def test_concatenates_multiple_tables():
    # one table per isolate, concatenated into a single frame; the global index
    # keeps composite ids unique even when (iso, family) repeats
    a = is_summary([("isoA", "IS1", 1, 2)])
    b = is_summary([("isoB", "IS1", 3, 4), ("isoB", "IS6", 5, 6)])
    out = preformat([a, b])
    assert len(out) == 3
    assert set(out["iso"]) == {"isoA", "isoB"}
    assert out.index.is_unique
