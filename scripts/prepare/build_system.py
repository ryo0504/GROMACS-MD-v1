#!/usr/bin/env python3
"""Thin CLI: generate one molecule; packing/solvation is an explicit separate step."""
import sys
from mdworkbench.cli import main
sys.argv.insert(1, "build")
main()
