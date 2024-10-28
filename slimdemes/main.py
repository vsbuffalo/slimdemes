import json
import sys
from pathlib import Path
from typing import Optional

import demes
import defopt
import yaml


def remove_gene_flow(graph):
    msg = "Removing '{block}' from demes graph, " "since --ignore-gene-flow is set."
    if "pulses" in graph:
        sys.stderr.write(msg.format(block="pulses"))
        graph.pop("pulses")
    if "migrations" in graph:
        sys.stderr.write(msg.format(block="migrations"))
        graph.pop("migrations")
    return graph


def load_demes(demes_file, ignore_gene_flow, convert_to_generations):
    with open(demes_file) as f:
        data = yaml.safe_load(f)

    if ignore_gene_flow:
        data = remove_gene_flow(data)

    graph = demes.Graph.fromdict(data)

    if convert_to_generations:
        graph = graph.in_generations()

    return graph


def msprime(
    input_file: Path,
    *,
    ignore_gene_flow: Optional[bool] = False,
    out: Optional[Path] = None,
) -> None:
    """Load the model into msprime.

    Args:
        input_file: Input YAML file path
        ignore_gene_flow: Remove gene flow blocks (pulses and migrations).
    """
    graph = load_demes(input_file, ignore_gene_flow, True)
    print(graph)


def convert(
    input_file: Path,
    *,
    ignore_gene_flow: Optional[bool] = False,
    out: Optional[Path] = None,
) -> None:
    """Convert YAML file to JSON.
    Args:
        input_file: Input YAML file path
        ignore_gene_flow: Remove gene flow blocks (pulses and migrations).
        out: Optional output JSON file path. If not provided, writes to stdout
    """
    graph = load_demes(
        input_file,
        ignore_gene_flow=ignore_gene_flow,
        convert_to_generations=True,
    )
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
    # Use ensure_ascii=False to handle any special characters
    json_str = json.dumps(data, indent=2, ensure_ascii=False)

    if out:
        with open(out, "w") as f:
            f.write(json_str)
    else:
        print(json_str)


def main() -> None:
    """CLI entry point."""
    defopt.run([convert, msprime])


if __name__ == "__main__":
    main()
