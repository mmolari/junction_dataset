import sys
from pathlib import Path

# Put scripts/ and scripts/annotations/ on sys.path so tests can import
# production modules (recombination_filter, path_utils, defensefinder_*, ...)
# without packaging them.
scripts = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(scripts))
sys.path.insert(0, str(scripts / "annotations"))
