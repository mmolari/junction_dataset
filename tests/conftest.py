import sys
from pathlib import Path

# Put scripts/ on sys.path so tests can import production modules
# (recombination_filter, path_utils, ...) without packaging them.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
