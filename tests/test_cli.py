from typer.testing import CliRunner
from smelt.cli import app

runner = CliRunner()


def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "site" in result.stdout


def test_site_help():
    result = runner.invoke(app, ["site", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.stdout


def test_window_help():
    result = runner.invoke(app, ["window", "--help"])
    assert result.exit_code == 0
    assert "--window" in result.stdout


def test_element_help():
    result = runner.invoke(app, ["element", "--help"])
    assert result.exit_code == 0
    assert "--gtf" in result.stdout


def test_metaplot_help():
    result = runner.invoke(app, ["metaplot", "--help"])
    assert result.exit_code == 0
    assert "--body-bins" in result.stdout


def test_custom_help():
    result = runner.invoke(app, ["custom", "--help"])
    assert result.exit_code == 0
    assert "--bed" in result.stdout


def test_dmr_help():
    result = runner.invoke(app, ["dmr", "--help"])
    assert result.exit_code == 0
    assert "--group1" in result.stdout


def test_stats_help():
    result = runner.invoke(app, ["stats", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.stdout


def test_site_without_input_fails():
    result = runner.invoke(app, ["site"])
    assert result.exit_code != 0


def test_site_no_fasta_fails():
    """Without FASTA, context can't be determined - should error."""
    result = runner.invoke(app, [
        "site", "--input", "tests/data/sample.cov",
    ])
    assert result.exit_code != 0
