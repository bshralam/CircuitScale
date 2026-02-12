#Resource metrics (depth, 2Q gates, parameter count)
from __future__ import annotations
from dataclasses import dataclass
from qiskit import QuantumCircuit

@dataclass
class CircuitResources:
    n_qubits: int
    n_params: int
    depth: int
    n_2q: int

def count_2q_gates(qc: QuantumCircuit) -> int:
    n = 0
    for inst, qargs, _ in qc.data:
        if len(qargs) == 2:
            n += 1
    return n

def circuit_resources(qc: QuantumCircuit) -> CircuitResources:
    return CircuitResources(
        n_qubits=qc.num_qubits,
        n_params=qc.num_parameters,
        depth=qc.depth(),
        n_2q=count_2q_gates(qc),
    )

