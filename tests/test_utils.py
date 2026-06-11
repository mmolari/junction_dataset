from utils import load_junction_positions

# Two edges, three rows; covers strand bool parsing and the nested {edge: {iso}} shape.
CSV = """edge,iso,left_start,left_end,right_start,right_end,strand
e1,isoA,100,200,300,400,True
e1,isoB,4841508,4855335,29951,39466,False
e2,isoA,10,20,30,40,True
"""


def _write(tmp_path):
    p = tmp_path / "jp.csv"
    p.write_text(CSV)
    return str(p)


def test_nested_by_edge_then_iso(tmp_path):
    # {edge: {iso: ...}} preserving every (edge, iso) pair
    d = load_junction_positions(_write(tmp_path))
    assert set(d) == {"e1", "e2"}
    assert set(d["e1"]) == {"isoA", "isoB"}


def test_tuple_order_and_strand_type(tmp_path):
    # (left_start, left_end, right_start, right_end, strand) with a real bool
    d = load_junction_positions(_write(tmp_path))
    assert d["e1"]["isoA"] == (100, 200, 300, 400, True)
    assert d["e1"]["isoB"][4] is False  # not the truthy string "False"


def test_index_single_junction(tmp_path):
    # the pattern extract_junctions uses: index the nested dict by edge id
    one = load_junction_positions(_write(tmp_path))["e1"]
    assert set(one) == {"isoA", "isoB"}
