"""StrainMake: reproducible hybrid metagenomics with MAG recovery and strain-level resolution."""

from pathlib import Path

__version__ = "0.2.0-pre1"

# Absolute path to the repository root (parent of this package directory).
# Used by CLI commands to locate workflow templates, config files, etc.
REPO_ROOT: Path = Path(__file__).parent.parent
