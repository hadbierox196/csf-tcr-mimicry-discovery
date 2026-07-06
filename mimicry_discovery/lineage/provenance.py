"""Per-run data/code/model provenance manifests.

Generic (not TCR-biology-specific) provenance tracking: hashes inputs,
captures the current git commit, and writes a JSON manifest alongside
pipeline outputs so any result can be traced back to the exact
inputs/code/model that produced it (see the architecture rationale on
data lineage). Unlike the other modules in this "stubs" tier, this one
is fully implemented -- it has no domain-specific unknowns, only
generic file/subprocess handling.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def hash_file(path: Path | str, algorithm: str = "sha256") -> str:
    """Compute a hex digest of a file's contents.

    Args:
        path: Path to the file to hash.
        algorithm: Any algorithm name accepted by :func:`hashlib.new`.

    Returns:
        The hex digest string.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Cannot hash missing file: {path}")
    hasher = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def current_git_sha(repo_dir: Path | str = ".") -> str | None:
    """Return the current git commit SHA, or None if unavailable.

    Args:
        repo_dir: Directory inside the git repository to check.

    Returns:
        The full commit SHA, or ``None`` if ``repo_dir`` is not inside a
        git repository (e.g. when running from an extracted archive) or
        git is not installed.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return result.stdout.strip()


def write_provenance_manifest(
    output_path: Path | str,
    input_paths: list[Path | str],
    tool_versions: dict[str, str],
    extra: dict[str, Any] | None = None,
) -> Path:
    """Write a JSON provenance manifest for one pipeline run.

    Args:
        output_path: Where to write the manifest JSON.
        input_paths: Input files whose hashes should be recorded.
        tool_versions: e.g. ``{"tcrdock": "1.2.0", "torch": "2.3.0"}``.
        extra: Any additional key-value metadata to include verbatim.

    Returns:
        The path the manifest was written to.

    Raises:
        FileNotFoundError: If any path in ``input_paths`` does not
            exist.
    """
    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_sha": current_git_sha(),
        "input_hashes": {str(p): hash_file(p) for p in input_paths},
        "tool_versions": tool_versions,
        "extra": extra or {},
    }
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2))
    return output_path
