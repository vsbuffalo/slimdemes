import random
import pyslim
from collections import defaultdict
from warnings import warn


def safe_float(value):
    """
    Return None if value is None/"None" otherwise
    convert to float.
    """
    if isinstance(value, str) and value.lower() == "none":
        return None
    if value is None:
        return None
    return float(value)


def subsample_by_population(ts, num_samples, seed=None):
    """
    Subsample individuals from each population.

    NOTE: to pass to simplify, you need to get the nodes
    associated with each individual.

    Args:
        ts: Tree sequence
        alive_inds: List of alive individual IDs
        num_samples: Number of samples to take from each population
        seed: Random seed for reproducibility

    Returns:
        List of sampled individual IDs
    """
    if seed is not None:
        random.seed(seed)

    alive_inds = pyslim.individuals_alive_at(ts, 0)

    # Group by population
    grouped_by_pop = defaultdict(list)
    for ind_id in alive_inds:
        ind = ts.individual(ind_id)
        grouped_by_pop[ind.population].append(ind_id)

    # Sample from each population
    sampled_ids = []
    for pop, ind_ids in grouped_by_pop.items():
        if len(ind_ids) < num_samples:
            warn(
                f"Population {pop} has fewer individuals ({len(ind_ids)}) "
                f"than requested samples ({num_samples}): {ind_ids}"
            )
        sampled_ids.extend(random.sample(ind_ids, min(num_samples, len(ind_ids))))

    return sampled_ids


def individuals_to_nodes(ts, individuals):
    nodes = []
    for i in individuals:
        nodes.extend(ts.individual(i).nodes)
    return nodes
