"""
Invoke tasks for Maptitude project management
"""

import os
import sys
from pathlib import Path

from invoke import task

# Project paths
PROJECT_ROOT = Path(__file__).parent.absolute()
BUILD_DIR = PROJECT_ROOT / "build"


@task
def build(ctx, preset="default"):
    """Build the project using CMake presets.

    :param preset: CMake preset to use (default or debug).
    """
    print(f"Configuring with preset '{preset}'...")
    ctx.run(f"cmake --preset {preset}", pty=True)

    print(f"Building with preset '{preset}'...")
    ctx.run(f"cmake --build --preset {preset}", pty=True)

    print("Build complete!")


@task
def test(ctx, preset="default"):
    """Run C++ tests.

    :param preset: CMake preset to use.
    """
    build_dir = BUILD_DIR if preset == "default" else PROJECT_ROOT / f"build-{preset}"

    print("Running C++ tests...")
    ctx.run(f"ctest --test-dir {build_dir} --output-on-failure", pty=True)


@task
def pytest(ctx):
    """Run Python tests."""
    print("Running Python tests...")
    ctx.run(f"{sys.executable} -m pytest tests/python -v", pty=True)


@task
def clean(ctx):
    """Remove build directories."""
    import shutil

    for d in ["build", "build-debug"]:
        p = PROJECT_ROOT / d
        if p.exists():
            print(f"Removing {p}...")
            shutil.rmtree(p)

    print("Clean complete!")
