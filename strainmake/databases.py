"""
Reference database bookkeeping shared by the CLI and the workflow.

The pipeline reaches out to the internet in two situations: when it
creates its conda environments, and when it downloads a reference database.
Clusters whose compute nodes have no route to the internet therefore need both
to happen somewhere else, ahead of time. This module backs the commands that
make that possible (`strainmake fetch-databases`, and the preflight check that
`strainmake run` performs), by answering three questions about a configuration:

* where does each database live,
* which of them will this run actually read,
* which of those are missing.

The path logic itself is *not* reimplemented here. It is imported from
`workflow/rules/utils.py`, the same module Snakemake loads, so the CLI and the
rules cannot disagree about where a database belongs.
"""

from __future__ import annotations

import importlib.util
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from strainmake import REPO_ROOT

# ---------------------------------------------------------------------------
# Loading the workflow's own helpers
# ---------------------------------------------------------------------------

_UTILS_PATH: Path = REPO_ROOT / "workflow" / "rules" / "utils.py"


def _load_workflow_utils():
    """
    Import `workflow/rules/utils.py` as a module.
    """

    if not _UTILS_PATH.exists():
        raise FileNotFoundError(
            f"Workflow helpers not found: {_UTILS_PATH}\n"
            "The StrainMake package must sit next to the `workflow/` directory."
        )

    spec = importlib.util.spec_from_file_location(
        "strainmake_workflow_utils", _UTILS_PATH
    )
    module = importlib.util.module_from_spec(spec)
    # registered so that repeated calls reuse one module object
    sys.modules.setdefault("strainmake_workflow_utils", module)
    spec.loader.exec_module(module)

    return module


_utils = _load_workflow_utils()

DATABASE_KEYS: List[str] = list(_utils.DATABASE_DEFAULT_PATHS)


# ---------------------------------------------------------------------------
# Describing a database
# ---------------------------------------------------------------------------


class ConfigError(Exception):
    """A configuration file that cannot be used, described for a human."""


def _describe_yaml_error(
    config_path: Path, text: str, error: yaml.YAMLError
) -> str:
    """
    Turn a YAML parser error into something a user can act on.

    The parser reports where it gave up, which is often not where the mistake
    is: a `key:value` written without a space is not a mapping entry at all but
    a plain scalar, so it swallows the following lines and the complaint lands
    on whichever later line first contains a colon. Quoting the reported line
    and naming that cause is usually enough to spot it.

    Parameters:
    config_path (Path): the file being read
    text (str): its contents
    error (yaml.YAMLError): what the parser raised

    Returns:
    str: a message naming the file, the line, and the likely cause
    """

    message = [f"{config_path} is not valid YAML."]

    mark = getattr(error, "problem_mark", None)
    problem = getattr(error, "problem", None)

    if problem:
        message.append(f"  {problem}")

    if mark is not None:
        lines = text.splitlines()
        message.append(f"  at line {mark.line + 1}, column {mark.column + 1}:")
        if 0 <= mark.line < len(lines):
            message.append(f"    {lines[mark.line]}")
            message.append("    " + " " * mark.column + "^")

        culprits = [
            (number, line)
            for number, line in enumerate(lines[: mark.line + 1], start=1)
            if re.match(r"^\s*[A-Za-z_][\w-]*:(?!//)\S", line)
        ]
        if culprits:
            number, line = culprits[-1]
            message.append(
                f"\n  Line {number} is missing a space after the colon, which makes"
                "\n  it a single string rather than a setting:"
                f"\n    {line.strip()}"
                f"\n  YAML needs `key: value`, not `key:value`."
            )

    return "\n".join(message)


@dataclass(frozen=True)
class DatabaseStatus:
    """The state of one reference database for a given configuration."""

    key: str
    label: str
    path: str
    target: str
    required: bool
    present: bool
    downloadable: bool
    manual_source: Optional[str] = None
    approx_size: Optional[str] = None

    @property
    def blocking(self) -> bool:
        """True when this database is needed by the run but is not there."""

        return self.required and not self.present


def load_config(config_path: Path) -> dict:
    """
    Read a StrainMake configuration YAML.

    A syntax error here is a typo in a file the user wrote, so it is reported as
    such (with the line, the offending text and the most common cause) rather
    than as a traceback through the YAML parser's internals.

    Parameters:
    config_path (Path): path to the configuration file

    Returns:
    dict: the parsed configuration

    Raises:
    ConfigError: when the file is not valid YAML, or is not a mapping
    """

    with open(config_path, encoding="utf-8") as fh:
        text = fh.read()

    try:
        config = yaml.safe_load(text)
    except yaml.YAMLError as error:
        raise ConfigError(_describe_yaml_error(config_path, text, error)) from None

    if config is None:
        return {}

    if not isinstance(config, dict):
        raise ConfigError(
            f"{config_path} does not describe a configuration.\n"
            f"  Expected a mapping of settings, found {type(config).__name__}."
        )

    return config


def database_status(config: dict, key: str) -> DatabaseStatus:
    """
    Describe one database: where it should be, and whether it is there.

    Parameters:
    config (dict): the pipeline configuration
    key (str): a database key, e.g. 'checkm2'

    Returns:
    DatabaseStatus: the resolved state of that database
    """

    target = _utils.database_target(config, key)

    return DatabaseStatus(
        key=key,
        label=_utils.DATABASE_LABELS.get(key, key),
        path=_utils.database_path(config, key),
        target=target,
        required=_utils.database_is_required(config, key),
        present=os.path.exists(target),
        downloadable=key in _utils.DOWNLOADABLE_DATABASES,
        manual_source=_utils.MANUAL_DATABASE_SOURCES.get(key),
        approx_size=_utils.DATABASE_APPROX_SIZES.get(key),
    )


def inspect_databases(config: dict) -> Dict[str, DatabaseStatus]:
    """
    Describe every database StrainMake knows about.

    Parameters:
    config (dict): the pipeline configuration

    Returns:
    dict: database key -> DatabaseStatus, in declaration order
    """

    return {key: database_status(config, key) for key in DATABASE_KEYS}


def missing_databases(config: dict) -> List[DatabaseStatus]:
    """
    The databases this run needs but cannot find.

    Returns:
    list: DatabaseStatus entries that would make the run fail
    """

    return [status for status in inspect_databases(config).values() if status.blocking]


def is_offline(config: dict) -> bool:
    """
    Whether the configuration declares that no network is available.

    Returns:
    bool: the value of `deployment: offline:`
    """

    return _utils.is_offline(config)


def conda_prefix(config: dict) -> Optional[str]:
    """
    The directory holding the pipeline's conda environments, if one is set.

    A shared prefix is what allows the environments to be built once, where
    there is internet access, and reused by runs that have none.

    Returns:
    str or None: the configured directory, or None to use Snakemake's default
    """

    deployment = config.get("deployment") or {}
    prefix = deployment.get("conda_prefix")

    return str(prefix) if prefix else None
