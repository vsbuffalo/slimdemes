import itertools
import scipy
import polars as pl
import numpy as np
import tskit
import msprime
import pyslim
from pathlib import Path
from typing import Optional
import matplotlib
from slimdemes.demes import load_demes_json
from slimdemes.utilities import subsample_by_population, individuals_to_nodes

matplotlib.use("Agg")


def check_coalescence(ts):
    """
    Check if tree sequence has fully coalesced.

    Args:
        ts: tskit.TreeSequence

    Returns:
        bool: True if fully coalesced, False otherwise

    Raises:
        ValueError: If any trees haven't coalesced
    """
    for tree in ts.trees():
        if tree.num_roots > 1:
            raise ValueError(
                f"Tree at position {tree.interval} has not fully coalesced "
                f"(num_roots={tree.num_roots})"
            )
    return True


def recapitate_if_needed(ts, random_seed=None):
    """
    Check if tree sequence has fully coalesced and recapitate if not.

    Args:
        ts: tskit.TreeSequence
        random_seed: Random seed for recapitation

    Returns:
        tskit.TreeSequence: Original or recapitated tree sequence
    """
    needs_recapitation = False
    num_uncoalesced = 0
    for tree in ts.trees():
        if tree.num_roots > 1:
            needs_recapitation = True
            num_uncoalesced += 1

    if needs_recapitation:
        print(f"Recapitating {num_uncoalesced} uncoalesced trees...")

        # Get recombination rate from metadata
        try:
            slim_md = ts.metadata["SLiM"]["user_metadata"]
        except IndexError:
            raise IndexError(
                "Metadata entry 'SLiM' > 'user_metadata' not found, "
                "is this a SLiM tree sequence?"
            )
        recmap = slim_md["recombination_map"][0]
        recomb_rates = recmap["rates"]
        recomb_ends = recmap["ends"]

        # If there's just one rate, use that
        if len(recomb_rates) == 1:
            recomb_rate = recomb_rates[0]
        else:
            # Create recombination map if multiple rates
            positions = [0] + recomb_ends
            rates = recomb_rates + [recomb_rates[-1]]  # repeat last rate
            recomb_rate = msprime.RateMap(position=positions, rate=rates)

        # Get initial Ne from first population in demes graph
        if "demes_json" in slim_md:
            demes_graph = load_demes_json(slim_md["demes_json"][0])
            # Use the ancestral population size
            demography = msprime.Demography.from_demes(demes_graph)
        else:
            raise ValueError("'demes_json' must be in SLiM metadata.")

        ts = pyslim.recapitate(
            ts,
            demography=demography,
            recombination_rate=recomb_rate,
            random_seed=random_seed,
        )

        print("Recapitation complete")

    return ts


def mimic_slimsim_with_msprime(
    slim_ts_path: Path,
    *,
    num_samples: Optional[int] = 10,
    random_seed: Optional[int] = None,
    out: Optional[Path] = None,
) -> None:
    """Run msprime simulation matching parameters from a SLiM tree sequence.

    Args:
        slim_ts_path: Path to SLiM tree sequence file to match parameters from
        num_samples: The number of samples to draw from each present-day deme
        random_seed: Random seed for reproducibility
        out: Output tree sequence file path (.trees)
    """
    # Load SLiM tree sequence to get metadata
    slim_ts = tskit.load(slim_ts_path)
    slim_md = slim_ts.metadata["SLiM"]["user_metadata"]

    # Load demographic model from metadata
    if "demes_json" not in slim_md:
        raise ValueError("'demes_json' must be in SLiM metadata")

    graph = load_demes_json(slim_md["demes_json"][0])
    demography = msprime.Demography.from_demes(graph)

    # Get recombination map from metadata
    if "recombination_map" not in slim_md:
        raise ValueError("'recombination_map' must be in SLiM metadata")

    recmap = slim_md["recombination_map"][0]
    recomb_rates = recmap["rates"]
    recomb_ends = recmap["ends"]

    # If there's just one rate, use that
    if len(recomb_rates) == 1:
        recomb_rate = recomb_rates[0]
    else:
        # Create recombination map if multiple rates
        positions = [0] + recomb_ends
        rates = recomb_rates + [recomb_rates[-1]]  # repeat last rate
        recomb_rate = msprime.RateMap(position=positions, rate=rates)

    # Get sequence length (note that slim returns arrays...)
    sequence_length = slim_md["contig_length"][0]

    # Set up sampling - similar to run_msprime but using max time
    samples = []
    for deme in graph.demes:
        if deme.end_time == 0:
            # sample from extent lineages only
            samples.extend(
                [msprime.SampleSet(num_samples, population=deme.name, time=0)]
            )

    # Run simulation
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recomb_rate,
        random_seed=None if random_seed is None else random_seed + 1,
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


def analyze_trees(ts_path, subsample_size=None, random_seed=None):
    """
    Analyze a tree sequence file and compute population genetic statistics

    Args:
        ts_path: Path to .trees file
        subsample_size: Optional int specifying number of samples to draw from each present-day population.
        random_seed: Random seed for recapitation and subsampling if needed

    Returns:
        dict of computed statistics including FST between all population pairs
    """
    ts = tskit.load(ts_path)

    # Recapitate if needed
    ts = recapitate_if_needed(ts, random_seed=random_seed)

    # Subsample if specified
    if subsample_size is not None:
        sample_individuals = subsample_by_population(
            ts, subsample_size, seed=random_seed
        )
        ts = ts.simplify(
            samples=individuals_to_nodes(ts, sample_individuals),
            filter_populations=False,
        )

    # Basic diversity statistics
    stats = {
        "nucleotide_diversity": ts.diversity(mode="branch"),
        "tajimas_d": ts.Tajimas_D(mode="branch"),
        "total_branch_length": sum(tree.total_branch_length for tree in ts.trees()),
    }

    # Compute site frequency spectrum
    afs = ts.allele_frequency_spectrum(
        polarised=True, span_normalise=False, mode="branch"
    )
    stats["afs"] = afs

    # Compute per-population diversity and AFS
    populations = list(set(node.population for node in ts.nodes() if node.time == 0))
    for pop in populations:
        pop_samples = [
            node.id for node in ts.nodes() if node.time == 0 and node.population == pop
        ]
        pop_ts = ts.simplify(pop_samples, filter_populations=False)
        stats[f"diversity_pop{pop}"] = pop_ts.diversity(mode="branch")

        afs = pop_ts.allele_frequency_spectrum(
            polarised=True, span_normalise=False, mode="branch"
        )
        stats[f"afs_pop{pop}"] = afs

    # Compute FST between all population pairs
    for pop1, pop2 in itertools.combinations(populations, 2):
        samples1 = [
            node.id for node in ts.nodes() if node.time == 0 and node.population == pop1
        ]
        samples2 = [
            node.id for node in ts.nodes() if node.time == 0 and node.population == pop2
        ]
        fst = ts.Fst([samples1, samples2], mode="branch")
        stats[f"fst_pop{pop1}_pop{pop2}"] = fst

    return stats


def compare_sims(sim_dir, model, rescale_q, num_samples, replicates):
    """
    Compare SLiM and msprime simulation results

    Args:
        sim_dir: Simulation results directory (for both msprime and slim)
        model: Model name
        rescale_q: Rescaling factor used in simulations
        num_samples: Number of samples used for msprime (SLiM sims sampled accordingly).
        replicates: Number of replicates to analyze
    """
    sim_engines = ("slim", "msprime")

    rows = []
    for rep in range(replicates):
        for sim_engine in sim_engines:
            if sim_engine == "msprime":
                tree_path = f"{sim_dir}/{model}_q{rescale_q}/{sim_engine}/n_{num_samples}/rep{rep}.trees"
                stats = analyze_trees(tree_path, random_seed=rep + 1)
            else:
                # SLiM tree of full population - we need to down-subsample
                tree_path = (
                    f"{sim_dir}/{model}_q{rescale_q}/{sim_engine}/rep{rep}.trees"
                )
                stats = analyze_trees(
                    tree_path, subsample_size=num_samples, random_seed=rep + 1
                )
            stats["rep"] = rep
            stats["sim_engine"] = sim_engine
            rows.append(stats)

    return pl.DataFrame(rows)


def compare_afs_distributions(afs1, afs2):
    """Compare two allele frequency spectra using Kolmogorov-Smirnov test.

    Args:
        afs1, afs2: Arrays containing allele frequency spectra

    Returns:
        tuple: (KS statistic, p-value)
    """
    from scipy import stats

    # Normalize the spectra so we're comparing distributions
    afs1_norm = afs1 / afs1.sum()
    afs2_norm = afs2 / afs2.sum()

    # Perform KS test
    ks_stat, p_value = stats.ks_2samp(afs1_norm, afs2_norm)

    return ks_stat, p_value


def compute_comparison_stats(stats_df: pl.DataFrame):
    """Compute statistical comparisons between SLiM and msprime results.

    Args:
        stats_df: DataFrame with statistics from compare_sims

    Returns:
        dict: Statistical test results for both regular stats and AFS
    """
    from scipy import stats
    import numpy as np

    results = {}

    # Separate AFS columns from other statistics
    afs_cols = [col for col in stats_df.columns if col.startswith("afs_pop")]
    other_stat_cols = [
        col
        for col in stats_df.columns
        if col not in ["rep", "sim_engine", "afs"] + afs_cols
    ]

    # Compute stats for regular statistics
    for stat in other_stat_cols:
        slim_values = stats_df.filter(pl.col("sim_engine") == "slim")[stat].to_numpy()
        msp_values = stats_df.filter(pl.col("sim_engine") == "msprime")[stat].to_numpy()

        # Perform t-test
        t_stat, p_val = stats.ttest_ind(slim_values, msp_values)

        results[stat] = {
            "test_type": "t-test",
            "test_stat": t_stat,
            "p_value": p_val,
            "slim_mean": np.mean(slim_values),
            "msp_mean": np.mean(msp_values),
            "relative_diff": abs(np.mean(slim_values) - np.mean(msp_values))
            / np.mean(msp_values),
            "significant": p_val < 0.05,
        }

    # Compute stats for AFS
    for afs_col in afs_cols:
        # Get mean AFS for each simulator
        slim_afs_list = stats_df.filter(pl.col("sim_engine") == "slim")[
            afs_col
        ].to_list()
        msp_afs_list = stats_df.filter(pl.col("sim_engine") == "msprime")[
            afs_col
        ].to_list()

        slim_afs = np.mean(
            np.array([np.array(afs[1:-1]) for afs in slim_afs_list]), axis=0
        )
        msp_afs = np.mean(
            np.array([np.array(afs[1:-1]) for afs in msp_afs_list]), axis=0
        )

        # Chi-square test for AFS
        chi2_stat, p_val = compare_afs_distributions(slim_afs, msp_afs)

        results[afs_col] = {
            "test_type": "chi-square",
            "test_stat": chi2_stat,
            "p_value": p_val,
            "slim_mean": np.mean([np.sum(afs) for afs in slim_afs_list]),
            "msp_mean": np.mean([np.sum(afs) for afs in msp_afs_list]),
            "relative_diff": abs(np.sum(slim_afs) - np.sum(msp_afs)) / np.sum(msp_afs),
            "significant": p_val < 0.05,
        }

    return results


def plot_stats_comparison(stats_df: pl.DataFrame, test_results: dict, output_path=None):
    """Plot comparison of statistics between SLiM and msprime using box plots.

    Args:
        stats_df: DataFrame with statistics from compare_sims
        test_results: Dict of statistical test results from compute_comparison_stats
        output_path: Optional path to save figure to
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Separate AFS columns from other statistics
    afs_cols = [col for col in stats_df.columns if col.startswith("afs_pop")]
    other_stat_cols = [
        col
        for col in stats_df.columns
        if col not in ["rep", "sim_engine", "afs"] + afs_cols
    ]

    # Calculate plot layout
    n_stats = len(other_stat_cols)
    n_afs = len(afs_cols)
    total_plots = n_stats + n_afs
    n_rows = (total_plots + 1) // 2

    fig, axes = plt.subplots(n_rows, 2, figsize=(12, 4 * n_rows))
    axes = axes.flatten()

    # Plot regular statistics
    for i, stat in enumerate(other_stat_cols):
        ax = axes[i]
        result = test_results[stat]

        # Get values for each simulation engine
        slim_values = stats_df.filter(pl.col("sim_engine") == "slim")[stat].to_numpy()
        msp_values = stats_df.filter(pl.col("sim_engine") == "msprime")[stat].to_numpy()

        # Create box plots
        box_data = [slim_values, msp_values]
        ax.boxplot(box_data, labels=["SLiM", "msprime"])

        # Add individual points
        x_slim = np.random.normal(1, 0.04, size=len(slim_values))
        x_msp = np.random.normal(2, 0.04, size=len(msp_values))
        ax.scatter(x_slim, slim_values, alpha=0.4, color="blue")
        ax.scatter(x_msp, msp_values, alpha=0.4, color="orange")

        # Add test results to plot title - color red if significant
        title_color = "red" if result["significant"] else "black"
        title = f"{stat}\n" f"t={result['test_stat']:.3f}, p={result['p_value']:.3f}"
        ax.set_title(title, color=title_color)
        ax.grid(True, linestyle="--", alpha=0.7)

    # Plot AFS for each population
    for j, afs_col in enumerate(afs_cols):
        ax = axes[n_stats + j]
        result = test_results[afs_col]

        # Plot AFS for each simulation type
        for sim_engine, color in zip(["slim", "msprime"], ["blue", "orange"]):
            afs_list = stats_df.filter(pl.col("sim_engine") == sim_engine)[
                afs_col
            ].to_list()
            afs_array = np.array([afs[1:-1] for afs in afs_list])
            mean_afs = np.mean(afs_array, axis=0)

            x = np.arange(1, len(mean_afs) + 1)
            ax.plot(x, mean_afs, "-o", label=sim_engine, color=color, markersize=3)

        ax.semilogx()
        ax.set_xlabel("Allele count")
        ax.set_ylabel("Number of variants")

        # Add test results to plot title
        title_color = "red" if result["significant"] else "black"
        title = (
            f"AFS Pop {afs_col.split('pop')[1]}\n"
            f"χ²={result['test_stat']:.3f}, p={result['p_value']:.3f}"
        )
        ax.set_title(title, color=title_color)
        ax.legend()
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.set_yscale("log")

    # Hide unused subplots
    for k in range(n_stats + len(afs_cols), len(axes)):
        axes[k].set_visible(False)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def check_significant_differences(
    test_results: dict, significance_level: float = 0.05
) -> bool:
    """Check if any tests show significant differences.

    Args:
        test_results: Dict of statistical test results from compute_comparison_stats
        significance_level: P-value threshold for significance

    Returns:
        bool: True if all tests pass (no significant differences), False otherwise
    """
    failed_tests = []

    for stat, result in test_results.items():
        if result["p_value"] < significance_level:
            failed_tests.append(f"{stat}: p={result['p_value']:.4f}")

    if failed_tests:
        print("\nSignificant differences found:")
        for test in failed_tests:
            print(f"- {test}")
        return False
    return True
