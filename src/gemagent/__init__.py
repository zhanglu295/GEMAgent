"""Standalone, general-purpose GEMAgent package.

The package depends on the installed official gapseq CLI, not on Codex skills.
"""
from .evaluation import evaluate_files, evaluate_frames
from .workflow import build, run_recipe, run_skill, __version__

__all__ = ["build", "run_recipe", "run_skill", "evaluate_files",
           "evaluate_frames", "__version__"]
