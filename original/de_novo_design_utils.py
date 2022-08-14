"""Some utility functions for de novo design"""

import math
import random
from typing import Dict, List, Tuple

import numpy as np
from rdkit import Chem, RDLogger


def generate_smiles(
    known_values: Dict[str, Tuple[float, Dict[str, float]]],
    k: int = 10_000,
    qed_change_threshold=0.05,
):
    """
    Randomly choose one of the top K SMILES, add 'C',
    and return if QED doesn't change much and it is still valid.
    """

    # Turn off logging (a lot of noise is generated here)
    RDLogger.DisableLog("rdApp.error")

    # Sort all known values which are not NaN
    dataset = [(v, k) for k, v in known_values.items() if not np.isnan(v[0])]
    dataset.sort(key=lambda x: x[0][0])

    # Try to make a new SMILES
    while True:
        chosen_entry = random.choice(dataset[:k])
        new_smiles = "C" + chosen_entry[1]
        new_mol = Chem.MolFromSmiles(new_smiles)

        # Only return if mol is ok
        if new_mol is not None:
            if math.isclose(Chem.QED.qed(new_mol), chosen_entry[0][1]["QED"], abs_tol=qed_change_threshold):
                return Chem.MolToSmiles(new_mol)


def topk_mols_over_time(
    results: Dict[str, Tuple[Tuple[float, Dict[str, float]], int]],
    times: List[int],
    k: int,
) -> Dict[int, List[Tuple[str, float, Dict[str, float]]]]:
    """Helper function to take a list of results and return the best molecules over time."""

    output = dict()
    for curr_time in times:
        assert curr_time >= 1

        # Sort all results at this time, excluding NaNs
        results_by_time_t = [(smiles, r, t) for smiles, (r, t) in results.items()
                             if t <= curr_time and not np.isnan(r[0])]
        results_by_time_t.sort(key=lambda x: x[1][0])

        # Add top k results from this time to the output
        current_output = []
        for i in range(k):
            if i < len(results_by_time_t):
                current_output.append((results_by_time_t[i][0], results_by_time_t[i][1][0], results_by_time_t[i][1][1]))
        output[curr_time] = current_output
    return output


def print_results_over_time(results_over_time):
    """Pretty printing the output of the function above."""
    for t, l in sorted(results_over_time.items()):
        print(f"\tTop molecules at time={t} function evaluations.")
        for rank, li in enumerate(l):
            print(f"\t\t#{rank+1}: obj={li[1]:.2g}, smiles={li[0]}")
