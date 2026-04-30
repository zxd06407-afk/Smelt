import pandas as pd
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


def test_site_unsupported_format():
    """Unsupported file extension gives a clear error."""
    result = runner.invoke(app, [
        "site", "--input", "tests/data/sample.bed",
    ])
    assert result.exit_code != 0


def test_site_cli_pipeline():
    """Run site via CLI with CX_report input, verify output file."""
    import tempfile, os
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        out_path = f.name
    try:
        result = runner.invoke(app, [
            "site", "--input", "tests/data/test.CX_report.txt",
            "--sample-name", "test",
            "-o", out_path,
        ])
        assert result.exit_code == 0
        df = pd.read_csv(out_path, sep="\t")
        assert len(df) > 0
        assert "sample" in df.columns
        assert all(df["sample"] == "test")
    finally:
        os.unlink(out_path)


def test_cli_error_message():
    """Nonexistent input gives non-zero exit."""
    result = runner.invoke(app, [
        "site", "--input", "nonexistent_file.bam",
    ])
    assert result.exit_code != 0
