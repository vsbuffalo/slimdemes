import json
from pathlib import Path
from typing import Optional

from moments.Demes.DemesUtil import rescale as rescale_demes
import msprime
import defopt
from slimdemes.demes import load_demes
from slimdemes.utilities import safe_float


def run_msprime(
    input_file: Path,
    *,
    ignore_gene_flow: Optional[bool] = False,
    rescale_q: Optional[float] = None,  # Added parameter
    num_samples: Optional[int] = 10,
    sequence_length: Optional[int] = 100_000,
    random_seed: Optional[int] = None,
    out: Optional[Path] = None,
) -> None:
    """Load the model into msprime and run simulation.

    Args:
        input_file: Input YAML file path
        ignore_gene_flow: Remove gene flow blocks (pulses and migrations)
        rescale_q: The scaling factor for population sizes (N' = N/q)
        num_samples: The number of samples to draw from each present-day deme.
        sequence_length: Length of sequence to simulate
        random_seed: Random seed for reproducibility
        out: Output tree sequence file path (.trees)
    """
    graph = load_demes(input_file, ignore_gene_flow, True)

    # Apply rescaling if specified
    if rescale_q is not None:
        assert rescale_q >= 1.0, "q in rescale_q must be >= 1.0"
        graph = rescale_demes(graph, rescale_q)

    # Convert demes graph to msprime demography
    demography = msprime.Demography.from_demes(graph)

    # Set up sampling
    num_samples = {
        deme.name: num_samples for deme in graph.demes if deme.end_time == 0}

    samples = []
    for deme_name, n in num_samples.items():
        samples.extend([msprime.SampleSet(n, population=deme_name, time=0)])

    # Run simulation
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        # msprime seeds must be >= 1
        random_seed=random_seed + 1,
    )

    # Save tree sequence
    if out:
        ts.dump(out)
    else:
        # Print summary statistics if no output file specified
        print("Tree sequence statistics:")
        print(f"Number of trees: {ts.num_trees}")
        print(f"Sample size: {ts.num_samples}")
        print(f"Sequence length: {ts.sequence_length}")


def convert(
    input_file: Path,
    *,
    ignore_gene_flow: Optional[bool] = False,
    rescale_q: Optional[float] = None,
    out: Optional[Path] = None,
) -> None:
    """Convert YAML file to JSON.
    Args:
        input_file: Input YAML file path
        ignore_gene_flow: Remove gene flow blocks (pulses and migrations).
        rescale_q: The scaling factor, e.g. setting q = 10 rescales the population size such that N' = N/10.
        out: Optional output JSON file path. If not provided, writes to stdout
    """
    graph = load_demes(
        input_file,
        ignore_gene_flow=ignore_gene_flow,
        convert_to_generations=True,
    )

    # rescale if necessary
    if rescale_q is not None:
        rescale_q = float(rescale_q)
        assert rescale_q >= 1.0, "q in --rescale-q must be >= 1.0"
        graph = rescale_demes(graph, rescale_q)

    data = graph.asdict()

    # Convert float('inf') to "Infinity" string
    def convert_infinity(obj):
        if isinstance(obj, dict):
            return {k: convert_infinity(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_infinity(v) for v in obj]
        elif isinstance(obj, float) and obj == float("inf"):
            return "Infinity"  # This will be properly quoted in JSON output
        return obj

    data = convert_infinity(data)

    # Include the rescaling factor so it can be used downstream
    data["metadata"] = dict(rescale_q=safe_float(rescale_q))

    # Use ensure_ascii=False to handle any special characters
    json_str = json.dumps(data, indent=2, ensure_ascii=False)

    if out:
        with open(out, "w") as f:
            f.write(json_str)
    else:
        print(json_str)


def main() -> None:
    """CLI entry point."""
    defopt.run([convert, run_msprime])


if __name__ == "__main__":
    main()
