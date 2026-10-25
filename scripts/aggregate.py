#!/usr/bin/env python3
"""Aggregate the raw output of the experiment drivers into the data files of the paper figures.

Usage:
    python3 scripts/aggregate.py <output directory> [<output directory> ...] [--out <directory>]

Each output directory is one directory written by a driver, e.g. `equilateral_45nodes_30seeds_30distanceOfCostZero_4sigma`.
The kind of experiment is recognized from the directory name. The aggregated files are written to `--out`
(default: `aggregated`) and are named like the data files of the paper figures.

For every configuration, the statistics are taken over the seeds:
- the percentage of fixed variables (of fixed triples for `tripletCut`) is 100 * (initial - final) / initial,
  where final is the number of remaining variables (triples) after the last step of a run,
- mean, std (ddof=1), median, q25, q75 (linear interpolation), min and max of these percentages,
- the runtime of a run is the sum of the runTime column [ms], as in the paper. Note that the drivers write the
  accumulated runtime of each criterion, so runs whose join phase is repeated count earlier passes more than once.
  For the runtimes, std uses ddof=0, as in the paper.
"""
import argparse
import os
import re
import sys

import pandas as pd

VARIABLES = "remainingVariables"
TRIPLES = "remainingTriples"


def fixedPercentage(df, key, column):
    """percentage of fixed variables or triples per configuration and seed"""
    runs = df.groupby(key + ["seed"], sort=True)
    numberOfRuns = df[df.criterionName == "initial"].groupby(key + ["seed"]).size()
    if (numberOfRuns != 1).any():
        raise RuntimeError("Some seeds were run more than once. The drivers append to existing files, "
                           "so delete the output directory before running a driver again.")
    initial = runs[column].first()
    final = runs[column].last()
    return (100 * (initial - final) / initial).rename("value").reset_index()


def runtime(df, key):
    return df[df.criterionName != "initial"].groupby(key + ["seed"], sort=True).runTime.sum().rename("value").reset_index()


def statistics(values, key, ddof=1, withMinMax=True):
    groups = values.groupby(key, sort=True).value
    stats = pd.DataFrame({
        "mean": groups.mean(),
        "std": groups.std(ddof=ddof),
        "median": groups.median(),
        "q25": groups.quantile(0.25),
        "q75": groups.quantile(0.75),
    })
    if withMinMax:
        stats["min"] = groups.min()
        stats["max"] = groups.max()
    return stats.reset_index()


def name(value):
    """parameter value as used in the file names of the paper data, e.g. 0.5, 0.01, 0.0"""
    return repr(float(value))


def write(df, outDir, fileName):
    path = os.path.join(outDir, fileName)
    df.to_csv(path, index=False)
    print(f"wrote {path}")


def readCsv(path):
    return pd.read_csv(path)


def aggregateTrianglesJoint(directory, outDir):
    df = readCsv(os.path.join(directory, "experiment.csv"))
    key = ["numberOfNodes", "sigma"]
    write(statistics(fixedPercentage(df, key, VARIABLES), key), outDir,
          "data-triangles-experiment_equilateral_triangles_all_vars.csv")
    write(statistics(runtime(df, ["sigma"]), ["sigma"], ddof=0, withMinMax=False), outDir,
          "data-triangles-experiment_runtimes.csv")


def criterionFiles(directory, prefix):
    pattern = re.compile(re.escape(prefix) + r"_([A-Za-z]+)\.csv$")
    for fileName in sorted(os.listdir(directory)):
        match = pattern.match(fileName)
        if match:
            yield match.group(1), os.path.join(directory, fileName)


def aggregateTrianglesIndividual(directory, outDir):
    key = ["numberOfNodes", "sigma"]
    for criterion, path in criterionFiles(directory, "experiment"):
        df = readCsv(path)
        write(statistics(fixedPercentage(df, key, VARIABLES), key), outDir,
              f"data-triangles-individual-experiment_{criterion}_vars.csv")
        if criterion == "tripletCut":
            write(statistics(fixedPercentage(df, key, TRIPLES), key), outDir,
                  f"data-triangles-individual-experiment_{criterion}_triplets.csv")


def withNumberOfPoints(df):
    df = df.copy()
    df["numberOfPoints"] = 8 * df["n"]
    return df


def aggregatePartitionJoint(directory, outDir):
    df = withNumberOfPoints(readCsv(os.path.join(directory, "experiment_all.csv")))
    key = ["numberOfPoints", "alpha"]
    for beta, part in df.groupby("beta"):
        write(statistics(fixedPercentage(part, key, VARIABLES), key), outDir,
              f"data-partition-experiment_vars_beta={name(beta)}.csv")
        write(statistics(runtime(part, ["alpha"]), ["alpha"], ddof=0, withMinMax=False), outDir,
              f"data-partition-experiment_runtimes_beta={name(beta)}.csv")


def aggregatePartitionIndividual(directory, outDir):
    key = ["numberOfPoints", "alpha"]
    for fileName in sorted(os.listdir(directory)):
        match = re.match(r"(experiment_full_n=\d+)_([A-Za-z]+)\.csv$", fileName)
        if not match:
            continue
        base, criterion = match.groups()
        df = withNumberOfPoints(readCsv(os.path.join(directory, fileName)))
        for beta, part in df.groupby("beta"):
            write(statistics(fixedPercentage(part, key, VARIABLES), key), outDir,
                  f"data-partition-individual-{base}_{criterion}_vars_beta={name(beta)}.csv")
            if criterion == "tripletCut":
                write(statistics(fixedPercentage(part, key, TRIPLES), key), outDir,
                      f"data-partition-individual-{base}_{criterion}_triples_beta={name(beta)}.csv")


def readInstanceSizes(directory):
    files = [f for f in os.listdir(directory) if re.match(r"\d+nodes(_all)?\.csv$", f)]
    return pd.concat([readCsv(os.path.join(directory, f)) for f in files], ignore_index=True)


def aggregateTrianglesInstanceSizes(directory, outDir):
    df = readInstanceSizes(directory).rename(columns={"numberOfNodes": "n"})
    for sigma, part in df.groupby("sigma"):
        write(statistics(fixedPercentage(part, ["n"], VARIABLES), ["n"], withMinMax=False), outDir,
              f"data-triangles-variables_sigma_{name(sigma)}.csv")
        write(statistics(runtime(part, ["n"]), ["n"], ddof=0, withMinMax=False), outDir,
              f"data-triangles-runtimes_sigma_{name(sigma)}.csv")


def aggregatePartitionInstanceSizes(directory, outDir):
    df = readInstanceSizes(directory)
    df["n"] = 8 * df["n"]
    for (alpha, beta), part in df.groupby(["alpha", "beta"]):
        suffix = f"alpha_{name(alpha)}_beta_{name(beta)}"
        write(statistics(fixedPercentage(part, ["n"], VARIABLES), ["n"], withMinMax=False), outDir,
              f"data-partition-variables_{suffix}.csv")
        write(statistics(runtime(part, ["n"]), ["n"], ddof=0, withMinMax=False), outDir,
              f"data-partition-runtimes_{suffix}.csv")


# the order matters: the first pattern that matches the directory name is used
EXPERIMENTS = [
    (r"equilateral_instance_sizes_", aggregateTrianglesInstanceSizes),
    (r"partition_instance_sizes_", aggregatePartitionInstanceSizes),
    (r"equilateral_individual_", aggregateTrianglesIndividual),
    (r"partition_individual_", aggregatePartitionIndividual),
    (r"equilateral_\d+nodes_", aggregateTrianglesJoint),
    (r"partition_\d+nodes_", aggregatePartitionJoint),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directories", nargs="+", help="output directories written by the drivers")
    parser.add_argument("--out", default="aggregated", help="directory for the aggregated files")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    for directory in args.directories:
        baseName = os.path.basename(os.path.normpath(directory))
        for pattern, aggregate in EXPERIMENTS:
            if re.match(pattern, baseName):
                aggregate(directory, args.out)
                break
        else:
            sys.exit(f"Cannot tell which experiment wrote {directory}.")


if __name__ == "__main__":
    main()
