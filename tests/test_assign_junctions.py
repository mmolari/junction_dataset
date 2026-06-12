import pytest

from assign_junctions import within, within_no_loop, within_loop, rotate

# All cases run on a length-100 circle. Two reference targets are reused throughout:
#   no-loop target [B, E] = [30, 60]  (does not wrap the origin)
#   loop    target [B, E] = [60, 30]  (wraps: covers [60, 100] u [0, 30])
L = 100


def wi(B, E, b, e):
    """within() is the public entry; it rotates wrapping queries (b > e) before classifying."""
    return within(B, E, b, e, L)


def test_rotate_wraps_modulo_length():
    assert rotate(95, 10, L) == 5  # past the origin
    assert rotate(5, -10, L) == 95  # negative delta wraps back
    assert rotate(30, 0, L) == 30  # zero delta is identity


# --- no-loop target [30, 60] -------------------------------------------------


def test_no_loop_basic_categories():
    assert wi(30, 60, 35, 45) == "subset"  # query strictly inside
    assert wi(30, 60, 35, 65) == "start"  # left end inside, runs off the right
    assert wi(30, 60, 10, 45) == "end"  # right end inside, starts left of B
    assert wi(30, 60, 25, 65) == "superset"  # query contains the target
    assert wi(30, 60, 80, 85) == "no"  # disjoint, right of E


def test_no_loop_wrapping_query():
    # b > e: the query itself wraps the origin; within() rotates it first
    assert wi(30, 60, 35, 10) == "start"
    assert wi(30, 60, 90, 45) == "end"
    assert wi(30, 60, 90, 70) == "superset"
    assert wi(30, 60, 25, 10) == "superset"
    assert wi(30, 60, 80, 10) == "no"


def test_no_loop_endpoint_touching_counts_as_overlap():
    # touching a boundary is overlap, not "no" -- pins the >=/<= tie-breaks
    assert wi(30, 60, 30, 45) == "subset"  # b == B
    assert wi(30, 60, 45, 60) == "subset"  # e == E
    assert wi(30, 60, 30, 60) == "subset"  # exact match
    assert wi(30, 60, 60, 70) == "start"  # b == E (just touches the right edge)
    assert wi(30, 60, 10, 30) == "end"  # e == B (just touches the left edge)


# --- loop target [60, 30] (wraps origin) -------------------------------------


def test_loop_basic_categories():
    assert wi(60, 30, 65, 90) == "subset"  # inside the [60,100] arm
    assert wi(60, 30, 10, 20) == "subset"  # inside the [0,30] arm
    assert wi(60, 30, 10, 40) == "start"  # left end inside [0,30], runs past E
    assert wi(60, 30, 35, 70) == "end"  # right end inside [60,100], starts in the hole
    assert wi(60, 30, 40, 50) == "no"  # entirely in the hole (30,60)


def test_loop_wrapping_query():
    # b > e: query wraps the origin, rotated before classifying
    assert wi(60, 30, 65, 10) == "subset"
    assert wi(60, 30, 90, 40) == "start"
    assert wi(60, 30, 35, 10) == "end"


def test_loop_endpoint_touching_counts_as_overlap():
    assert wi(60, 30, 10, 30) == "subset"  # e == E
    assert wi(60, 30, 60, 70) == "subset"  # b == B
    assert wi(60, 30, 90, 30) == "subset"  # spans the origin, both ends on the target
    assert wi(60, 30, 30, 40) == "start"  # b == E (touches inner edge of the hole)


# --- degenerate / ambiguous geometries (pin current behavior) ----------------


def test_wraparound_query_raises():
    # a query that enters the target from both sides without being contained is
    # geometrically ambiguous; within_loop refuses to classify it.
    with pytest.raises(AssertionError):
        wi(30, 60, 55, 35)  # rotated no-loop query lands as an ambiguous loop
    with pytest.raises(AssertionError):
        wi(60, 30, 10, 90)  # straddles both arms of the wrapped target
    with pytest.raises(AssertionError):
        wi(60, 30, 25, 65)


def test_query_fully_outside_target_is_no():
    # the "??" cases flagged in the source: query occupies the regions either side
    # of a no-loop target but never the target itself -> no overlap
    assert wi(30, 60, 70, 25) == "no"
    assert wi(30, 60, 61, 29) == "no"
    # and a query sitting entirely in the hole of a wrapped target
    assert wi(60, 30, 31, 59) == "no"


def test_zero_length_query_inside_target_is_subset():
    # b == e (a point) inside the target classifies as subset
    assert wi(30, 60, 45, 45) == "subset"
    assert wi(60, 30, 70, 70) == "subset"


# --- underlying helper contracts ---------------------------------------------


def test_within_no_loop_rejects_wrapped_target():
    # within_no_loop assumes an ordered target; B > E is a caller error
    with pytest.raises(AssertionError):
        within_no_loop(60, 30, 35, 45)


def test_within_loop_rejects_ordered_target():
    # within_loop assumes a wrapped target; B <= E is a caller error
    with pytest.raises(AssertionError):
        within_loop(30, 60, 35, 45)
