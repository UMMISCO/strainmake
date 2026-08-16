"""
Tests for `strainmake run` - Snakemake wrapper.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from strainmake import REPO_ROOT
from strainmake.cli.main import app
from strainmake.cli.run_cmd import _DEFAULT_SNAKEFILE, _PREPROCESSED_SNAKEFILE

runner = CliRunner()


def _invoke_run(args: list, returncode: int = 0, skip_checks: bool = True):
    """
    Invoke the run command with subprocess.run mocked.

    The preflight is skipped by default so that tests about how the Snakemake
    command line is assembled do not depend on which reference databases happen
    to exist on the machine running them. Tests that are about the preflight
    itself pass `skip_checks=False`.
    """
    mock_result = MagicMock()
    mock_result.returncode = returncode
    # the Snakemake invocation lives in the shared _snakemake helper, which is
    # where subprocess is actually called from
    if skip_checks and "--skip-checks" not in args:
        # anything after a bare `--` belongs to Snakemake, so the flag has to go
        # before it to reach StrainMake itself
        cut = args.index("--") if "--" in args else len(args)
        args = args[:cut] + ["--skip-checks"] + args[cut:]
    with patch("strainmake.cli._snakemake.subprocess.run", return_value=mock_result) as mock_sub:
        result = runner.invoke(app, ["run"] + args, catch_exceptions=False)
    return result, mock_sub


class TestDefaultSnakefile:
    def test_default_snakefile_points_to_workflow(self):
        """The default Snakefile should be the one in the workflow directory."""
        assert _DEFAULT_SNAKEFILE == REPO_ROOT / "workflow" / "Snakefile"

    def test_preprocessed_snakefile_exists(self):
        """The Snakefile for already preprocessed sequences should exist."""
        assert _PREPROCESSED_SNAKEFILE == REPO_ROOT / "workflow" / "Snakefile.already_preprocessed_seq.smk"


class TestRunCommand:
    def test_calls_snakemake(self):
        """Running `strainmake run` should call Snakemake."""
        result, mock_sub = _invoke_run(["--cores", "2"])
        assert mock_sub.called
        cmd = mock_sub.call_args[0][0]
        assert "snakemake" in cmd

    def test_use_conda_always_present(self):
        """The `--use-conda` flag should always be included when calling Snakemake, regardless of other options."""
        _, mock_sub = _invoke_run(["--cores", "4"])
        cmd = mock_sub.call_args[0][0]
        assert "--use-conda" in cmd

    def test_cores_forwarded(self):
        """The number of cores specified with `--cores` should be forwarded to Snakemake."""
        _, mock_sub = _invoke_run(["--cores", "16"])
        cmd = mock_sub.call_args[0][0]
        assert "--cores" in cmd
        assert "16" in cmd

    def test_dry_run_flag(self):
        """The `--dry-run` flag should be forwarded to Snakemake when specified."""
        _, mock_sub = _invoke_run(["--dry-run"])
        cmd = mock_sub.call_args[0][0]
        assert "--dry-run" in cmd

    def test_configfile_forwarded(self, tmp_path):
        """The `--configfile` option should be forwarded to Snakemake with the correct path."""
        cfg = tmp_path / "config.yaml"
        cfg.write_text("samples: data/metadata.tsv\n")
        _, mock_sub = _invoke_run(["--configfile", str(cfg)])
        cmd = mock_sub.call_args[0][0]
        assert "--configfile" in cmd
        assert str(cfg) in cmd

    def test_already_preprocessed_uses_correct_snakefile(self):
        """The `--already-preprocessed` flag should use the correct Snakefile."""
        _, mock_sub = _invoke_run(["--already-preprocessed"])
        cmd = mock_sub.call_args[0][0]
        assert "already_preprocessed" in " ".join(cmd)

    def test_extra_args_forwarded(self):
        """Extra arguments after `--` should be forwarded to Snakemake."""
        _, mock_sub = _invoke_run(["--cores", "4", "--", "--rerun-incomplete"])
        cmd = mock_sub.call_args[0][0]
        assert "--rerun-incomplete" in cmd

    def test_exit_code_propagated(self):
        """If the Snakemake subprocess returns a non-zero exit code, `strainmake run` should also exit with a non-zero code."""
        result, _ = _invoke_run([], returncode=1)
        assert result.exit_code == 1

    def test_default_cores_is_4(self):
        """If no `--cores` option is provided, the default should be 4 cores."""
        _, mock_sub = _invoke_run([])
        cmd = mock_sub.call_args[0][0]
        idx = cmd.index("--cores")
        assert cmd[idx + 1] == "4"


class TestOfflineMode:
    """
    The no-internet workflow splits into a preparation step, run where there is
    network access, and the run itself, which must not need any. `run` has to
    reach for the prepared environments and refuse to start rather than fail
    hours later against an unreachable host.
    """

    @staticmethod
    def _offline_config(
        tmp_path: Path, databases_present: bool, envs_built: bool = True
    ) -> Path:
        """A config with a shared env prefix and, optionally, everything staged."""
        import yaml

        from strainmake.databases import DATABASE_KEYS, database_status

        envs = tmp_path / "envs"
        if envs_built:
            envs.mkdir()

        cfg = {
            "samples": "data/metadata.tsv",
            "databases": {
                key: str(tmp_path / "dbs" / key) for key in DATABASE_KEYS
            },
            "deployment": {"conda_prefix": str(envs), "offline": False},
        }

        if databases_present:
            for key in DATABASE_KEYS:
                target = Path(database_status(cfg, key).target)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.touch()

        path = tmp_path / "offline_config.yaml"
        path.write_text(yaml.dump(cfg, sort_keys=False))
        return path

    def test_conda_prefix_is_forwarded(self, tmp_path):
        """The shared environment directory has to reach Snakemake."""
        cfg = self._offline_config(tmp_path, databases_present=True)
        _, mock_sub = _invoke_run(["--config", str(cfg)])
        cmd = mock_sub.call_args[0][0]
        assert "--conda-prefix" in cmd
        assert cmd[cmd.index("--conda-prefix") + 1] == str(tmp_path / "envs")

    def test_no_conda_prefix_when_unconfigured(self):
        """Without the setting, Snakemake's own default is left alone."""
        _, mock_sub = _invoke_run(["--cores", "2"])
        cmd = mock_sub.call_args[0][0]
        assert "--conda-prefix" not in cmd

    def test_create_envs_only_forwards_the_snakemake_flag(self, tmp_path):
        """This is the login-node preparation step."""
        cfg = self._offline_config(tmp_path, databases_present=True)
        _, mock_sub = _invoke_run(["--config", str(cfg), "--create-envs-only"])
        cmd = mock_sub.call_args[0][0]
        assert "--conda-create-envs-only" in cmd

    def test_create_envs_only_ignores_missing_databases(self, tmp_path):
        """
        Building environments reads no database, and it is the first command a
        user runs, before anything has been fetched. Checking for databases
        here would block the whole workflow at step one.
        """
        cfg = self._offline_config(tmp_path, databases_present=False)
        result, mock_sub = _invoke_run(skip_checks=False, args=["--config", str(cfg), "--create-envs-only"])
        assert result.exit_code == 0
        assert mock_sub.called

    def test_offline_tells_the_workflow_it_is_offline(self, tmp_path):
        """The workflow needs to know, so its download rules can fail fast."""
        cfg = self._offline_config(tmp_path, databases_present=True)
        _, mock_sub = _invoke_run(["--config", str(cfg), "--offline"])
        cmd = mock_sub.call_args[0][0]
        override = cmd[cmd.index("--config") + 1]
        assert '"offline": true' in override

    def test_offline_override_preserves_other_deployment_settings(self, tmp_path):
        """--config replaces a whole top-level key, so the rest must be resent."""
        cfg = self._offline_config(tmp_path, databases_present=True)
        _, mock_sub = _invoke_run(["--config", str(cfg), "--offline"])
        cmd = mock_sub.call_args[0][0]
        override = cmd[cmd.index("--config") + 1]
        assert "conda_prefix" in override

    def test_offline_refuses_to_start_with_missing_databases(self, tmp_path):
        """Offline, a missing database can never be resolved by running."""
        cfg = self._offline_config(tmp_path, databases_present=False)
        result, mock_sub = _invoke_run(skip_checks=False, args=["--config", str(cfg), "--offline"])
        assert result.exit_code == 1
        assert not mock_sub.called

    def test_offline_refuses_to_start_with_unbuilt_environments(self, tmp_path):
        """The environments cannot be solved without network access either."""
        cfg = self._offline_config(
            tmp_path, databases_present=True, envs_built=False
        )
        result, mock_sub = _invoke_run(skip_checks=False, args=["--config", str(cfg), "--offline"])
        assert result.exit_code == 1
        assert not mock_sub.called

    def test_offline_refuses_without_a_shared_prefix(self, tmp_path):
        """A per-project environment directory defeats the whole approach."""
        import yaml

        cfg_path = self._offline_config(tmp_path, databases_present=True)
        cfg = yaml.safe_load(cfg_path.read_text())
        cfg["deployment"]["conda_prefix"] = ""
        cfg_path.write_text(yaml.dump(cfg, sort_keys=False))

        result, mock_sub = _invoke_run(skip_checks=False, args=["--config", str(cfg_path), "--offline"])
        assert result.exit_code == 1
        assert not mock_sub.called

    def test_online_only_warns_about_missing_databases(self, tmp_path):
        """Online, the pipeline downloads what it needs, as it always has."""
        cfg = self._offline_config(tmp_path, databases_present=False)
        result, mock_sub = _invoke_run(skip_checks=False, args=["--config", str(cfg)])
        assert result.exit_code == 0
        assert mock_sub.called

    def test_offline_falls_back_to_the_snakefile_configuration(self, tmp_path):
        """
        Every Snakefile StrainMake runs names its own config, so the user should
        not have to pass it twice, including for --offline, which needs it to
        locate the prebuilt environments.
        """
        cfg = self._offline_config(tmp_path, databases_present=True)
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text(f'configfile: "{cfg}"\n')

        _, mock_sub = _invoke_run(
            skip_checks=False,
            args=["--offline", "--snakefile", str(snakefile), "--cores", "2"],
        )
        cmd = mock_sub.call_args[0][0]
        assert cmd[cmd.index("--conda-prefix") + 1] == str(tmp_path / "envs")

    def test_offline_refuses_when_no_configuration_can_be_found(self, tmp_path):
        """
        Without a configuration the environment directory is unknown, and the
        offline override would replace the whole `deployment` mapping, so the
        run could not work. Refuse rather than start it.
        """
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text("rule all:\n    input: []\n")

        result, mock_sub = _invoke_run(
            skip_checks=False,
            args=["--offline", "--snakefile", str(snakefile), "--cores", "2"],
        )
        assert result.exit_code == 1
        assert not mock_sub.called

    def test_skip_checks_bypasses_preflight(self, tmp_path):
        """The escape hatch must actually reach Snakemake."""
        cfg = self._offline_config(tmp_path, databases_present=False)
        _, mock_sub = _invoke_run(["--config", str(cfg), "--offline", "--skip-checks"])
        assert mock_sub.called


class TestHistoricalUsage:
    """
    The pipeline was used online long before any of the offline machinery
    existed, with commands that pass no configuration at all. Those must keep
    working, and must not silently acquire new behaviour.
    """

    def test_plain_run_still_works(self):
        """`strainmake run --cores N`, the command the README has always shown."""
        _, mock_sub = _invoke_run(["--cores", "16"])
        cmd = mock_sub.call_args[0][0]
        assert cmd[:2] == ["snakemake", "--snakefile"]
        assert "--use-conda" in cmd
        assert cmd[cmd.index("--cores") + 1] == "16"

    def test_plain_run_does_not_pin_a_conda_prefix(self):
        """
        Without `deployment: conda_prefix:` Snakemake keeps managing the
        environments itself, exactly as before.
        """
        _, mock_sub = _invoke_run(["--cores", "16"])
        assert "--conda-prefix" not in mock_sub.call_args[0][0]

    def test_plain_run_does_not_declare_itself_offline(self):
        """The workflow's download rules stay enabled for an ordinary run."""
        _, mock_sub = _invoke_run(["--cores", "16"])
        assert "--config" not in mock_sub.call_args[0][0]

    def test_plain_run_is_still_checked(self, tmp_path):
        """
        A run that passes no configuration still gets the database checks, by
        reading the configuration the Snakefile itself names. Otherwise the
        most common invocation of all would be the one invocation with no
        warning before a multi-hundred-gigabyte download.
        """
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            "samples: data/metadata.tsv\n"
            f"databases: {{checkm2: {tmp_path / 'dbs' / 'checkm2'}}}\n"
            "binning: {binner: [metabat2]}\n"
        )
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text(f'configfile: "{cfg}"\n')

        result, _ = _invoke_run(
            skip_checks=False,
            args=["--snakefile", str(snakefile), "--cores", "4"],
        )
        assert "CheckM2" in result.stderr
        assert "will be downloaded" in result.stderr
