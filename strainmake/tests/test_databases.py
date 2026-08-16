"""
Tests for reference database resolution and the offline deployment path.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from strainmake.databases import (
    DATABASE_KEYS,
    ConfigError,
    conda_prefix,
    database_status,
    inspect_databases,
    is_offline,
    load_config,
    missing_databases,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_CONFIG = REPO_ROOT / "config" / "template_config.yaml"


@pytest.fixture(name="config")
def _fixture_config() -> dict:
    """The shipped template configuration."""
    with open(TEMPLATE_CONFIG, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------


def test_each_database_is_located_independently(config: dict):
    """
    Reference databases are not kept together on a real system, so pointing one
    somewhere else must not move any of the others.
    """
    config = copy.deepcopy(config)
    config["databases"]["gtdbtk"] = "/opt/gtdb/release220"
    config["databases"]["bakta"] = "/scratch/bakta"

    assert database_status(config, "gtdbtk").path == "/opt/gtdb/release220"
    assert database_status(config, "bakta").path == "/scratch/bakta"
    # untouched entries keep their own defaults
    assert database_status(config, "checkm2").path == "data/databases/checkm2"


def test_unset_databases_fall_back_to_their_own_defaults(config: dict):
    """A project that configures nothing still resolves every database."""
    for key in DATABASE_KEYS:
        assert database_status(config, key).path.startswith("data/databases/")


def test_missing_databases_section_falls_back_to_defaults():
    """A config written before the `databases:` section still resolves."""
    status = database_status({}, "checkm2")

    assert status.path == "data/databases/checkm2"


def test_meteor_target_follows_the_configured_catalogue(config: dict):
    """The Meteor reference is the catalogue subdirectory, not the download dir."""
    config = copy.deepcopy(config)
    config["databases"]["meteor"] = "/dbs/meteor"
    config["databases"]["meteor_catalogue"] = "hs_8_4_oral"

    assert database_status(config, "meteor").target == "/dbs/meteor/hs_8_4_oral"


def test_bakta_target_points_inside_the_db_subdirectory(config: dict):
    """`bakta_db download` unpacks into a `db/` subdirectory."""
    config = copy.deepcopy(config)
    config["databases"]["bakta"] = "/dbs/bakta"

    assert database_status(config, "bakta").target == "/dbs/bakta/db/version.json"


# ---------------------------------------------------------------------------
# Requirement detection
# ---------------------------------------------------------------------------


def test_removing_a_tool_section_drops_its_database(config: dict):
    """Deleting a tool's config section is how a user opts out of its database."""
    config = copy.deepcopy(config)
    assert database_status(config, "metaphlan").required

    del config["taxonomic_profiling"]["metaphlan"]
    assert not database_status(config, "metaphlan").required


def test_host_index_is_required_by_default(config: dict):
    """Host decontamination feeds nearly every downstream step."""
    assert database_status(config, "bowtie2_index").required


def test_a_database_can_be_declared_unused(config: dict):
    """
    `bowtie2:` holds settings the preprocessing rules read at parse time, so it
    cannot be deleted from the config the way a tool section can. A run that
    starts from already-preprocessed reads never reads the host index though,
    and must be able to say so, otherwise it is blocked offline by a database
    it does not use.
    """
    config = copy.deepcopy(config)
    config["databases"]["bowtie2_index"] = False

    assert not database_status(config, "bowtie2_index").required


def test_declaring_a_database_unused_does_not_affect_others(config: dict):
    config = copy.deepcopy(config)
    config["databases"]["bowtie2_index"] = False

    assert database_status(config, "checkm2").required


def test_strainscan_is_declared_but_not_downloadable(config: dict):
    """
    StrainScan ships only per-species Google Drive links.

    It must therefore be reported as a manual step rather than attempted, while
    still being checked for, since the shipped config profiles with it.
    """
    status = database_status(config, "strainscan")

    assert status.required
    assert not status.downloadable
    assert "StrainScan" in (status.manual_source or "")


def test_every_other_database_is_downloadable(config: dict):
    """Everything except StrainScan can be fetched unattended."""
    downloadable = {
        key for key, status in inspect_databases(config).items() if status.downloadable
    }

    assert downloadable == set(DATABASE_KEYS) - {"strainscan"}


# ---------------------------------------------------------------------------
# Presence detection
# ---------------------------------------------------------------------------


def test_absent_database_is_reported_as_blocking(config: dict, tmp_path: Path):
    """A required database that is not on disk stops the run."""
    config = copy.deepcopy(config)
    for key in DATABASE_KEYS:
        config["databases"][key] = str(tmp_path / "empty" / key)

    blocking = {status.key for status in missing_databases(config)}

    assert blocking == set(DATABASE_KEYS)


def test_present_database_is_not_blocking(config: dict, tmp_path: Path):
    """Creating the marker file is enough to satisfy the check."""
    config = copy.deepcopy(config)
    config["databases"]["checkm2"] = str(tmp_path / "dbs" / "checkm2")

    target = Path(database_status(config, "checkm2").target)
    target.parent.mkdir(parents=True)
    target.touch()

    blocking = {status.key for status in missing_databases(config)}

    assert "checkm2" not in blocking


# ---------------------------------------------------------------------------
# Deployment settings
# ---------------------------------------------------------------------------


def test_offline_defaults_to_false(config: dict):
    """A normal run is unaffected by the offline machinery."""
    assert not is_offline(config)


def test_conda_prefix_defaults_to_unset(config: dict):
    """Without a shared prefix, Snakemake's per-project default is used."""
    assert conda_prefix(config) is None


def test_conda_prefix_is_configurable(config: dict):
    """A shared prefix is what lets environments be built once and reused."""
    config = copy.deepcopy(config)
    config["deployment"]["conda_prefix"] = "/shared/strainmake/envs"

    assert conda_prefix(config) == "/shared/strainmake/envs"


def test_conda_prefix_treats_empty_string_as_unset(config: dict):
    """The shipped template ships the key empty rather than absent."""
    config = copy.deepcopy(config)
    config["deployment"]["conda_prefix"] = ""

    assert conda_prefix(config) is None


def test_load_config_reads_the_template():
    """The shipped template parses and carries the new sections."""
    config = load_config(TEMPLATE_CONFIG)

    assert "databases" in config
    assert "deployment" in config


# ---------------------------------------------------------------------------
# Conda environments
# ---------------------------------------------------------------------------


def test_environments_do_not_declare_variables():
    """
    `variables:` in an environment file are baked in when the environment is
    built, so with a shared `conda_prefix` they would freeze one project's
    database paths for every other project using those environments. Database
    locations belong in the `databases:` config section instead.
    """
    envs_dir = REPO_ROOT / "workflow" / "envs"
    offenders = [
        env.name
        for env in sorted(envs_dir.glob("*.yaml"))
        if "\nvariables:" in env.read_text()
    ]

    assert not offenders, (
        f"{offenders} declare `variables:`; move those settings to the "
        "`databases:` section of the config instead"
    )


# ---------------------------------------------------------------------------
# The documented walkthrough
# ---------------------------------------------------------------------------


def _walkthrough_config(tmp_path: Path) -> dict:
    """
    The configuration of the step-by-step example in the documentation.

    It is the most demanding real scenario we have: it starts from
    already-preprocessed reads (so the host index is never read), skips
    taxonomic profiling entirely, and post-processes with GTDB-Tk, CheckM1,
    CarveMe and Bakta. Only the parts that decide which databases are needed
    are reproduced here.
    """

    return {
        "samples": str(tmp_path / "samples.tsv"),
        # kept in the config because the preprocessing rules read it at parse
        # time, even though this run starts from preprocessed reads
        "bowtie2": {"index_name": "chm13v2.0", "threads": 4},
        "assembly": {"hybrid_assembler": ["hybridspades"]},
        "binning": {"binner": ["metabat2", "vamb"]},
        "checkm2": {"threads": 4},
        "bins_postprocessing": {
            "gtdbtk": {"threads": 1},
            "drep": {"ani": [95, 97]},
            "carveme": {"threads": 5, "solver": "gurobi"},
            "bakta": {"threads": 10, "parallel_jobs": 2},
        },
        # the walkthrough performs no taxonomic profiling: the section is absent
        "strains_profiling": {"instrain": {}},
        "databases": {
            "checkm2": str(tmp_path / "dbs" / "checkm2"),
            "bakta": str(tmp_path / "dbs" / "bakta"),
            "gtdbtk": str(tmp_path / "dbs" / "gtdb_tk"),
            # already-preprocessed reads: the host index is never read
            "bowtie2_index": False,
        },
    }


def test_walkthrough_needs_exactly_three_databases(tmp_path: Path):
    """
    The documented example must not be asked for databases it never reads.

    Getting this wrong is not cosmetic: offline, a database judged required but
    absent blocks the run outright, so an over-broad answer here makes the
    documented walkthrough impossible to follow on a cluster.
    """
    config = _walkthrough_config(tmp_path)

    required = {
        key for key, status in inspect_databases(config).items() if status.required
    }

    assert required == {"checkm2", "bakta", "gtdbtk"}


def test_walkthrough_databases_are_all_downloadable(tmp_path: Path):
    """`strainmake fetch-databases` can prepare the walkthrough unattended."""
    config = _walkthrough_config(tmp_path)

    manual = [
        status.label
        for status in inspect_databases(config).values()
        if status.required and not status.downloadable
    ]

    assert not manual


# ---------------------------------------------------------------------------
# Reading the configuration
# ---------------------------------------------------------------------------


def test_malformed_yaml_is_reported_not_raised(tmp_path: Path):
    """
    A typo in a hand-written config is a user error, not a crash.

    Left alone, PyYAML's exception surfaces as a long traceback through the
    parser's internals, which says nothing about what to change.
    """
    config = tmp_path / "config.yaml"
    config.write_text("samples:/data/samples.tsv\nlr_seq_format: fastq\n")

    with pytest.raises(ConfigError) as raised:
        load_config(config)

    assert "not valid YAML" in str(raised.value)


def test_a_key_glued_to_its_value_is_named(tmp_path: Path):
    """
    The parser blames the line it gave up on, not the line at fault.

    `key:value` without a space is a plain scalar, so it absorbs the lines that
    follow and the error is reported against the first of those containing a
    colon. The message has to point back at the real culprit.
    """
    config = tmp_path / "config.yaml"
    config.write_text("samples:/data/samples.tsv\nlr_seq_format: fastq\n")

    with pytest.raises(ConfigError) as raised:
        load_config(config)

    message = str(raised.value)
    assert "Line 1" in message
    assert "missing a space after the colon" in message


def test_a_bare_url_is_not_mistaken_for_a_glued_key(tmp_path: Path):
    """`https://...` has the same shape as `key:value` but is not the problem."""
    config = tmp_path / "config.yaml"
    config.write_text("a: 1\nhttps://data.gtdb.org/x.tar.gz\nb: [unclosed\n")

    with pytest.raises(ConfigError) as raised:
        load_config(config)

    assert "missing a space after the colon" not in str(raised.value)


def test_urls_in_a_valid_config_are_untouched(tmp_path: Path):
    """The shipped template contains both https:// and docker:// values."""
    config = tmp_path / "config.yaml"
    config.write_text(
        'gtdbtk_release_url: "https://data.gtdb.org/r220.tar.gz"\n'
        "deployment:\n  container: docker://x/y:1\n"
    )

    assert load_config(config)["gtdbtk_release_url"].startswith("https://")


def test_an_empty_config_is_empty_not_broken(tmp_path: Path):
    config = tmp_path / "config.yaml"
    config.write_text("")

    assert load_config(config) == {}


def test_a_config_that_is_not_a_mapping_is_rejected(tmp_path: Path):
    """A list parses fine as YAML but is not a configuration."""
    config = tmp_path / "config.yaml"
    config.write_text("- one\n- two\n")

    with pytest.raises(ConfigError):
        load_config(config)


# ---------------------------------------------------------------------------
# Rerun triggers
# ---------------------------------------------------------------------------


def _workflow_utils():
    """The workflow helper module, as Snakemake loads it."""
    from strainmake import databases

    return databases._utils


def test_offline_guard_is_identical_online_and_offline(tmp_path: Path):
    """
    The guard is a rule parameter, and `params` is one of Snakemake's default
    rerun triggers.

    If its text differed between an online and an offline run, going offline
    would mark every download rule out of date, and re-running a download rule
    offline fails, even when the database is already on disk. A run that had
    just fetched its databases would break the moment it went offline.
    """
    utils = _workflow_utils()
    base = {"databases": {"checkm2": str(tmp_path / "checkm2")}}

    online = dict(base, deployment={"offline": False})
    offline = dict(base, deployment={"offline": True})

    assert utils.offline_guard(online, "CheckM2", "checkm2") == utils.offline_guard(
        offline, "CheckM2", "checkm2"
    )


def test_offline_guard_tests_the_environment_at_run_time():
    """Offline-ness has to reach the rule without going through its params."""
    guard = _workflow_utils().offline_guard({}, "CheckM2", "checkm2")

    assert "STRAINMAKE_OFFLINE" in guard
    assert "exit 1" in guard


def test_offline_guard_still_changes_when_the_database_moves(tmp_path: Path):
    """
    Relocating a database *should* invalidate a previous download, so the guard
    is allowedNand expected, to differ then.
    """
    utils = _workflow_utils()
    here = {"databases": {"checkm2": str(tmp_path / "a")}}
    there = {"databases": {"checkm2": str(tmp_path / "b")}}

    assert utils.offline_guard(here, "CheckM2", "checkm2") != utils.offline_guard(
        there, "CheckM2", "checkm2"
    )
