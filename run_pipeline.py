"""Validate and execute the Rapalog discovery pipeline."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

NOTEBOOKS = [
    "01_exploratory_analysis.ipynb",
    "02_fragment_and_similarity.ipynb",
    "03_unsupervised_analysis.ipynb",
    "04_similarity_prediction.ipynb",
    "05_docking_pipeline.ipynb",
    "06_binding_affinity_prediction.ipynb",
]

REQUIRED_PACKAGES = [
    "Bio",
    "matplotlib",
    "meeko",
    "numpy",
    "pandas",
    "pyarrow",
    "rdkit",
    "shap",
    "sklearn",
    "umap",
    "vina",
]

REQUIRED_SOURCE_FILES = [
    "src/config.py",
    "src/data_processing.py",
    "src/fragment_extraction.py",
    "src/visualization.py",
    "src/docking.py",
]


def report_file(relative_path: str, problems: list[str]) -> None:
    """Report whether a required project file exists."""
    path = ROOT / relative_path
    status = "OK" if path.exists() else "MISSING"
    print(f"[{status}] {relative_path}")

    if not path.exists():
        problems.append(relative_path)


def validate_setup(require_coconut: bool = False) -> int:
    """Validate inputs, source files, dependencies, and executables."""
    problems: list[str] = []

    print(f"Project root: {ROOT}")
    print("\nReference structures:")

    for relative_path in [
        "data/raw/1FKB.cif",
        "data/raw/1FKB.pdb",
    ]:
        report_file(relative_path, problems)

    print("\nExternal dataset:")

    coconut_path = ROOT / "data/raw/coconut_csv-04-2026.csv"
    coconut_exists = coconut_path.exists()

    if coconut_exists:
        print(f"[OK] {coconut_path.relative_to(ROOT)}")
    else:
        status = "MISSING" if require_coconut else "OPTIONAL-MISSING"
        print(f"[{status}] data/raw/coconut_csv-04-2026.csv")

        if require_coconut:
            problems.append("data/raw/coconut_csv-04-2026.csv")

    print("\nSource files:")

    for relative_path in REQUIRED_SOURCE_FILES:
        report_file(relative_path, problems)

    print("\nNotebooks:")

    for notebook in NOTEBOOKS:
        report_file(f"notebooks/{notebook}", problems)

    print("\nPython packages:")

    for package in REQUIRED_PACKAGES:
        available = importlib.util.find_spec(package) is not None
        status = "OK" if available else "MISSING"
        print(f"[{status}] {package}")

        if not available:
            problems.append(f"Python package: {package}")

    print("\nExternal executables:")

    for executable in ["obabel", "vina"]:
        executable_path = shutil.which(executable)

        if executable_path:
            print(f"[OK] {executable}: {executable_path}")
        else:
            print(f"[MISSING] {executable}")
            problems.append(f"Executable: {executable}")

    if problems:
        print("\nSetup validation failed. Missing requirements:")

        for problem in problems:
            print(f"  - {problem}")

        return 1

    print("\nSetup validation passed.")
    return 0


def execute_notebooks(stop_after: int | None = None) -> int:
    """Execute the notebooks sequentially and stop on the first error."""
    selected = NOTEBOOKS if stop_after is None else NOTEBOOKS[:stop_after]

    for notebook in selected:
        notebook_path = ROOT / "notebooks" / notebook
        print(f"\nExecuting notebooks/{notebook}")

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--inplace",
                "--ExecutePreprocessor.timeout=-1",
                str(notebook_path),
            ],
            cwd=ROOT,
            check=False,
        )

        if result.returncode != 0:
            print(f"Notebook execution failed: {notebook}")
            return result.returncode

    print("\nAll selected notebooks completed successfully.")
    return 0


def main() -> int:
    """Parse command-line arguments and run the requested operation."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--check",
        action="store_true",
        help="validate project files, dependencies, and executables",
    )
    parser.add_argument(
        "--require-coconut",
        action="store_true",
        help="treat the external COCONUT CSV as required",
    )
    parser.add_argument(
        "--run-notebooks",
        action="store_true",
        help="execute notebooks 01 through 06 in place",
    )
    parser.add_argument(
        "--stop-after",
        type=int,
        choices=range(1, 7),
        metavar="N",
        help="stop after notebook N",
    )

    args = parser.parse_args()

    if not args.check and not args.run_notebooks:
        parser.print_help()
        return 0

    if args.check:
        status = validate_setup(require_coconut=args.require_coconut)

        if status != 0:
            return status

    if args.run_notebooks:
        return execute_notebooks(stop_after=args.stop_after)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())