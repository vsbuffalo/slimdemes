import sys
import json
import demes
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


def load_demes_json(demes_json_file):
    """
    Load a demes graph from a JSON file.
    Parameters
    ----------
    demes_json_file : str
        Path to JSON file containing demes graph specification
    Returns
    -------
    demes.Graph
        The loaded demes graph
    """

    def convert_infinity_strings(obj):
        if isinstance(obj, dict):
            return {k: convert_infinity_strings(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_infinity_strings(v) for v in obj]
        elif isinstance(obj, str) and obj == "Infinity":
            return float("inf")
        return obj

    with open(demes_json_file) as f:
        data = json.load(f)

    # Convert any "Infinity" strings to float('inf')
    data = convert_infinity_strings(data)

    graph = demes.Graph.fromdict(data)
    return graph
