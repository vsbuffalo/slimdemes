import pytest
import os
import yaml
from pathlib import Path
from snakemake.api import DAGSettings
from snakemake.api import ResourceSettings
from snakemake.api import SnakemakeApi
from snakemake.settings.types import OutputSettings
from snakemake.settings.types import Quietness
import polars as pl
from scipy import stats

# Significance level for statistical tests
ALPHA = 0.05


def get_workflow_rules(snakefile: Path) -> list[str]:
    """
    Get list of rules from a Snakefile using Snakemake API.

    Args:
        snakefile: Path to the Snakefile

    Returns:
        List of rule names
    """
    with SnakemakeApi(OutputSettings(quiet={Quietness.ALL})) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            snakefile=snakefile,
            resource_settings=ResourceSettings(),
            workdir=snakefile.parent,
        )

        # Get all rules from the workflow object
        workflow = workflow_api._workflow
        return workflow.rules


def get_rule_all_targets(snakefile: Path) -> list[str]:
    """
    Get target files specifically from the 'all' rule in a Snakefile.

    Args:
        snakefile: Path to the Snakefile

    Returns:
        List of target files specified in the 'all' rule
    """
    with SnakemakeApi(OutputSettings(quiet={Quietness.ALL})) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            snakefile=snakefile,
            resource_settings=ResourceSettings(),
            workdir=snakefile.parent,
        )

        # Get workflow object and find 'all' rule
        workflow = workflow_api._workflow
        rules = workflow.rules

        all_rule = [r for r in rules if r.name == "all"]

        if len(all_rule) == 0:
            return []

        assert len(all_rule) == 1
        all_rule = all_rule[0]

        # Get input files from 'all' rule
        return all_rule.input


@pytest.fixture(scope="session")
def resource_settings():
    """
    Configure resource settings based on environment.

    Priorities:
    1. Environment variable SNAKEMAKE_CORES
    2. GitHub Actions environment (use max available cores)
    3. Default to 1 core as fallback
    """
    # Check for explicit core count from environment
    cores_env = os.getenv("SNAKEMAKE_CORES")

    if cores_env is not None:
        cores = int(cores_env)
    # Check if running in GitHub Actions
    elif os.getenv("GITHUB_ACTIONS") == "true":
        # GitHub Actions provides NUMBER_OF_PROCESSORS
        cores = int(os.getenv("NUMBER_OF_PROCESSORS", "2"))
    else:
        # Default fallback
        cores = 1

    return ResourceSettings(
        cores=cores,
    )


@pytest.fixture(scope="session", autouse=True)
def run_snakemake(resource_settings):
    """Session-scoped fixture to run Snakemake workflow before any tests."""
    workflow_dir = Path("workflows/msprime_validation")
    snakefile = workflow_dir / "Snakefile"

    if not workflow_dir.exists():
        raise RuntimeError(
            f"Workflow directory not found at {workflow_dir}. "
            "Are you running tests from the project root?"
        )

    with SnakemakeApi(OutputSettings(quiet={Quietness.ALL})) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            snakefile=snakefile,
            resource_settings=resource_settings,
            # config_settings=ConfigSettings(config=config_args),
            workdir=snakefile.parent,
        )
        all_targets = get_rule_all_targets(Path("Snakefile"))
        dag_api = workflow_api.dag(
            dag_settings=DAGSettings(targets=all_targets),
        )
        dag_api.unlock()
        dag_api.execute_workflow()


@pytest.fixture()
def stats_files():
    """Create all the expected output files."""
    return get_rule_all_targets(snakefile)


# ----- main tests ------
@pytest.fixture()
def output_files(stats_files):
    """Convert stats files to parameters for testing."""
    return [pytest.param(Path(file), id=Path(file).name) for file in stats_files]


@pytest.mark.parametrize("output_file", "output_files", indirect=True)
def test_file_exists(output_file):
    """Test if each individual output file exists."""
    assert output_file.exists(), f"Output file {output_file} was not created"


if __name__ == "__main__":
    # For testing.
    snakefile = Path("workflows/msprime_validation/Snakefile")

    all_rules = get_workflow_rules(snakefile)
    print("All rules:", all_rules)

    # Get specific targets from 'all' rule
    all_rule_targets = get_rule_all_targets(snakefile)
    print("Rule 'all' targets:", all_rule_targets)

    run_snakemake()
