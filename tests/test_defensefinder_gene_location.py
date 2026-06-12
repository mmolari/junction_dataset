import pandas as pd
import pytest

from defensefinder_gene_location import protein_coordinates, add_gene_locations

# protein_coordinates parses prodigal/pyrodigal FASTA headers of the form
#   ">id # begin # end # strand # ...extra"
# add_gene_locations joins those coordinates onto a genes table by hit_id.


def write_prt(tmp_path, text):
    p = tmp_path / "proteins.prt"
    p.write_text(text)
    return str(p)


def test_parses_begin_end_strand_as_ints(tmp_path):
    prt = write_prt(tmp_path, ">g1 # 3 # 119 # 1 # ID=1_1;partial=10\nMSEQ\n")
    coords = protein_coordinates(prt)
    assert coords["g1"] == (3, 119, 1)
    assert all(isinstance(x, int) for x in coords["g1"])


def test_parses_negative_strand(tmp_path):
    prt = write_prt(tmp_path, ">g2 # 200 # 500 # -1 # ID=1_2\nMSEQ\n")
    coords = protein_coordinates(prt)
    assert coords["g2"] == (200, 500, -1)


def test_ignores_extra_hash_in_trailing_field(tmp_path):
    # only the first three numeric fields are read; a '#' in the ID/note tail is harmless
    prt = write_prt(tmp_path, ">g3 # 600 # 700 # 1 # ID=1_3;note=has#hash\nMSEQ\n")
    coords = protein_coordinates(prt)
    assert coords["g3"] == (600, 700, 1)


def test_keyed_by_record_id(tmp_path):
    # the dict key is the FASTA record id (first whitespace token), not the full header
    prt = write_prt(
        tmp_path, ">g1 # 1 # 2 # 1 # ID=a\nMA\n>g2 # 3 # 4 # 1 # ID=b\nMB\n"
    )
    coords = protein_coordinates(prt)
    assert set(coords) == {"g1", "g2"}


def test_add_gene_locations_joins_by_hit_id():
    coords = {"g1": (3, 119, 1), "g2": (200, 500, -1)}
    # rows deliberately out of coords order to confirm per-row lookup
    genes = pd.DataFrame({"hit_id": ["g2", "g1"], "gene_name": ["a", "b"]})
    out = add_gene_locations(genes, coords)
    # existing columns preserved, three location columns appended
    assert list(out.columns) == ["hit_id", "gene_name", "gene_beg", "gene_end", "gene_strand"]
    assert out.loc[0, ["gene_beg", "gene_end", "gene_strand"]].tolist() == [200, 500, -1]
    assert out.loc[1, ["gene_beg", "gene_end", "gene_strand"]].tolist() == [3, 119, 1]


def test_add_gene_locations_missing_hit_id_raises():
    # a hit_id absent from the protein FASTA is a hard error, not a silent NaN
    coords = {"g1": (3, 119, 1)}
    with pytest.raises(KeyError):
        add_gene_locations(pd.DataFrame({"hit_id": ["missing"]}), coords)
