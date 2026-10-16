"""Clean saved notebook outputs without altering code or Markdown cells."""

from __future__ import annotations

import sys
from pathlib import Path

import nbformat


PROJECT_ROOT = Path.cwd().resolve()
NOTEBOOK_DIR = PROJECT_ROOT / "notebooks"

REPLACEMENTS = (
    (str(PROJECT_ROOT), "<PROJECT_ROOT>"),
    (str(Path(sys.prefix).resolve()), "<CONDA_ENV>"),
    (str(Path.home().resolve()), "<HOME>"),
)


def sanitize_text(value):
    """Replace machine-specific paths in strings or lists of strings."""
    if isinstance(value, str):
        for original, replacement in REPLACEMENTS:
            value = value.replace(original, replacement)
        return value

    if isinstance(value, list):
        return [sanitize_text(item) for item in value]

    return value


def clean_notebook(path: Path) -> tuple[int, int]:
    """Remove saved stderr streams and sanitize paths in remaining outputs."""
    notebook = nbformat.read(path, as_version=4)
    removed_stderr = 0
    changed_outputs = 0

    for cell in notebook.cells:
        if cell.get("cell_type") != "code":
            cell.pop("outputs", None)
            cell.pop("execution_count", None)
            continue

        cleaned_outputs = []

        for output in cell.get("outputs", []):
            if output.get("output_type") == "stream" and output.get("name") == "stderr":
                removed_stderr += 1
                continue

            before = repr(output)

            if "text" in output:
                output["text"] = sanitize_text(output["text"])

            if "traceback" in output:
                output["traceback"] = sanitize_text(output["traceback"])

            data = output.get("data", {})
            for mime_type in ("text/plain", "text/markdown", "text/html"):
                if mime_type in data:
                    data[mime_type] = sanitize_text(data[mime_type])

            if repr(output) != before:
                changed_outputs += 1

            cleaned_outputs.append(output)

        cell["outputs"] = cleaned_outputs

    nbformat.write(notebook, path)
    return removed_stderr, changed_outputs


def main() -> None:
    for path in sorted(NOTEBOOK_DIR.glob("*.ipynb")):
        removed, changed = clean_notebook(path)
        print(
            f"{path.name}: removed {removed} stderr stream(s), "
            f"sanitized {changed} output(s)"
        )


if __name__ == "__main__":
    main()