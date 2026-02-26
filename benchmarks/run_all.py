"""Run all maptitude vs bms-bio benchmarks.

Usage:
    cd /path/to/maptitude
    PYTHONPATH=python:$PYTHONPATH python benchmarks/run_all.py
"""

import os
import sys

# Ensure the benchmarks directory is on sys.path so helpers can be imported
_bench_dir = os.path.dirname(os.path.abspath(__file__))
if _bench_dir not in sys.path:
    sys.path.insert(0, _bench_dir)


def _banner(title: str):
    print()
    print("#" * 72)
    print(f"#  {title}")
    print("#" * 72)
    print()


def main():
    _banner("BENCHMARK 1: Fc Density Generation")
    from bench_fc_density import run as run_fc
    run_fc()

    _banner("BENCHMARK 2: Scoring Metrics — Speed + Agreement")
    from bench_scoring import run as run_scoring
    run_scoring()

    _banner("BENCHMARK 3: Accuracy vs Published References")
    from bench_accuracy import run as run_accuracy
    run_accuracy()

    print("=" * 72)
    print("All benchmarks complete.")
    print("=" * 72)


if __name__ == "__main__":
    main()
