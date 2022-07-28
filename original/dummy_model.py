from typing import List
from sklearn.linear_model import LinearRegression
from rdkit import Chem

ATOM_LIST = ["C", "O", "H", "S", "N", "F", "Cl", "Br", "I", "P"]


def smiles_to_atom_counts(smiles: str) -> List[int]:
    """
    Given a SMILES string, produce a vector of its atom counts (including hydrogens).

    Examples:

    "CCC" is C3H8
    "CCC" -> [3, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0]

    "CC[OH]" is C2O1H6
    "CC[OH]" -> [2, 1, 6, 0, 0, 0, 0, 0, 0, 0, 0]
    """

    # List the atom types which we are counting explicitly
    atom_to_idx = {symbol: i for i, symbol in enumerate(ATOM_LIST)}

    # Make the mol object (including implicit hydrogens)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Count the atoms
    atom_counts = [0 for _ in ATOM_LIST] + [0]  # extra atom
    for atom in mol.GetAtoms():
        # if it isn't in the dictionary, it is counted
        # as an "extra" atom at the end
        idx = atom_to_idx.get(atom.GetSymbol(), -1)
        atom_counts[idx] += 1

    return atom_counts


def get_model() -> LinearRegression:
    """Creates and returns a new linear regression model"""
    return LinearRegression()
