#Build the two ansätze (HEA + UCCSD) + count resources
from __future__ import annotations
from typing import Tuple, List
import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

def hea(n_qubits: int, layers: int = 1, entangle: str = "linear") -> QuantumCircuit:
    """
    Simple HEA: Ry on all qubits + entanglers each layer.
    """
    theta = ParameterVector("theta", length=n_qubits * layers)
    qc = QuantumCircuit(n_qubits)
    k = 0
    for _ in range(layers):
        for q in range(n_qubits):
            qc.ry(theta[k], q)
            k += 1

        if entangle == "linear":
            for q in range(n_qubits - 1):
                qc.cx(q, q + 1)
        elif entangle == "full":
            for i in range(n_qubits):
                for j in range(i + 1, n_qubits):
                    qc.cx(i, j)
        else:
            raise ValueError("entangle must be 'linear' or 'full'")
    return qc

def uccsd_placeholder(n_qubits: int) -> QuantumCircuit:
    """
    Minimal placeholder so the pipeline runs before you wire in real UCCSD.
    We'll replace this with Qiskit Nature UCCSD *after* the plotting pipeline works.
    """
    qc = QuantumCircuit(n_qubits)
    # no parameters; it's a placeholder
    return qc



def hf_bitstring(n_qubits: int, n_electrons: int) -> list[int]:
    """
    Return occupied spin-orbital indices for a simple RHF reference
    under Jordan–Wigner ordering: occupy the lowest n_electrons orbitals.
    """
    return list(range(n_electrons))

def hea_with_hf_init(
    n_qubits: int,
    n_electrons: int,
    layers: int = 1,
    entangle: str = "linear",
) -> QuantumCircuit:
    """
    Prepare |HF> (fill lowest n_electrons qubits with X) then apply HEA layers.
    """
    theta = ParameterVector("theta", length=n_qubits * layers)
    qc = QuantumCircuit(n_qubits)

    # HF reference (JW ordering, lowest orbitals occupied)
    for q in hf_bitstring(n_qubits, n_electrons):
        qc.x(q)

    k = 0
    for _ in range(layers):
        for q in range(n_qubits):
            qc.ry(theta[k], q)
            k += 1

        if entangle == "linear":
            for q in range(n_qubits - 1):
                qc.cx(q, q + 1)
        elif entangle == "full":
            for i in range(n_qubits):
                for j in range(i + 1, n_qubits):
                    qc.cx(i, j)
        else:
            raise ValueError("entangle must be 'linear' or 'full'")

    return qc

