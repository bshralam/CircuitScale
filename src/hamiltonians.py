from __future__ import annotations
from typing import Tuple

from openfermion import QubitOperator, MolecularData
from openfermion.transforms import jordan_wigner
from openfermion.utils import count_qubits
from openfermionpyscf import run_pyscf

from .systems import MoleculeSpec


def molecular_qubit_hamiltonian(spec: MoleculeSpec) -> Tuple[QubitOperator, int]:
    """
    Robust across openfermionpyscf versions:
    - Build openfermion.MolecularData
    - run_pyscf(molecule, ...)
    - get molecular Hamiltonian
    - Jordan–Wigner map to qubits
    """
    # OpenFermion expects geometry as: [("H", (x,y,z)), ...]
    geometry = [(sym, tuple(xyz)) for sym, xyz in spec.geometry]

    multiplicity = spec.spin + 1  # spec.spin is 2S (PySCF convention)

    mol = MolecularData(
        geometry=geometry,
        basis=spec.basis,
        multiplicity=multiplicity,
        charge=spec.charge,
        filename=None,  # keeps it in-memory
    )

    # Depending on version, these flags are accepted; keep minimal.
    mol = run_pyscf(mol, run_scf=True, run_fci=False)

    fermion_ham = mol.get_molecular_hamiltonian()
    qubit_ham = jordan_wigner(fermion_ham)
    n_qubits = count_qubits(qubit_ham)
    return qubit_ham, n_qubits

