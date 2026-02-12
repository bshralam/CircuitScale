from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List

import numpy as np
from scipy.optimize import minimize

from openfermion import QubitOperator

from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit import QuantumCircuit


def qubitop_to_sparsepauliop(qubit_op: QubitOperator, n_qubits: int) -> SparsePauliOp:
    """
    Convert OpenFermion QubitOperator -> Qiskit SparsePauliOp.
    """
    pauli_strs: List[str] = []
    coeffs: List[complex] = []

    # Qiskit Pauli strings are ordered from highest qubit index to 0
    # We build a string like "IXYZ..." length n_qubits.
    for term, coeff in qubit_op.terms.items():
        s = ["I"] * n_qubits
        for (q, p) in term:  # p in {'X','Y','Z'}
            s[q] = p
        pauli_str = "".join(reversed(s))
        pauli_strs.append(pauli_str)
        coeffs.append(coeff)

    return SparsePauliOp(pauli_strs, coeffs=np.array(coeffs, dtype=complex))


def expectation_statevector(op: SparsePauliOp, qc: QuantumCircuit, x: np.ndarray) -> float:
    bound = qc.assign_parameters(x, inplace=False)
    psi = Statevector.from_instruction(bound)
    val = psi.expectation_value(op)
    return float(np.real(val))


@dataclass
class VQEResult:
    energy: float
    x_opt: np.ndarray
    nfev: int
    success: bool
    message: str


def run_vqe(
    op: SparsePauliOp,
    ansatz: QuantumCircuit,
    x0: Optional[np.ndarray] = None,
    maxiter: int = 200,
    seed: int = 7,
) -> VQEResult:
    n_params = ansatz.num_parameters
    rng = np.random.default_rng(seed)

    if x0 is None:
        # small random init; HEA is periodic so wide init not necessary
        x0 = 0.01 * rng.standard_normal(n_params)

    fun = lambda x: expectation_statevector(op, ansatz, x)

    res = minimize(
        fun,
        x0,
        method="COBYLA",
        options={"maxiter": maxiter, "disp": False},
    )

    return VQEResult(
        energy=float(res.fun),
        x_opt=np.array(res.x, dtype=float),
        nfev=int(getattr(res, "nfev", -1)),
        success=bool(res.success),
        message=str(res.message),
    )

