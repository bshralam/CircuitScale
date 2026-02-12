from __future__ import annotations
from typing import List, Dict, Any, Tuple
import re
import pandas as pd

from .systems import MoleculeSpec
from .references import compute_reference
from .hamiltonians import molecular_qubit_hamiltonian
from .ansatze import hea, uccsd_placeholder
from .resources import circuit_resources


def _parse_h4_delta_from_name(name: str) -> float:
    """
    Expected name like: H4_sq_d+0.10_a1.00
    """
    m = re.search(r"_d([+-]?\d+\.\d+)", name)
    if not m:
        raise ValueError(f"Could not parse delta from name: {name}")
    return float(m.group(1))


def run_sweep_h4(
    specs: List[MoleculeSpec],
    hea_layers=(1, 2, 3),
    ncas: int = 4,
    nelecas: Tuple[int, int] = (2, 2),
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    for spec in specs:
        delta = _parse_h4_delta_from_name(spec.name)

        # ----- Reference / MR metric
        ref = compute_reference(
            spec,
            do_fci=True,        # still feasible for H4 STO-3G
            do_casscf=True,
            ncas=ncas,
            nelecas=nelecas,
        )

        # ----- Hamiltonian / qubits
        qham, n_qubits = molecular_qubit_hamiltonian(spec)

        # ----- HEA resources at different depths
        for L in hea_layers:
            qc = hea(n_qubits, layers=L, entangle="linear")
            res = circuit_resources(qc)
            rows.append({
                "system": "H4",
                "name": spec.name,
                "param": delta,          # distortion parameter
                "mr_score": ref.mr_score,
                "ansatz": "HEA",
                "layers": L,
                "n_qubits": res.n_qubits,
                "n_params": res.n_params,
                "depth": res.depth,
                "n_2q": res.n_2q,
                "E_HF": ref.ehf,
                "E_FCI": ref.e_fci,
                # store NOONs (up to 4 for the active space)
                "noon_0": float(ref.noon[0]),
                "noon_1": float(ref.noon[1]),
                "noon_2": float(ref.noon[2]),
                "noon_3": float(ref.noon[3]),
            })

        # ----- UCCSD placeholder (swap later)
        qc_u = uccsd_placeholder(n_qubits)
        res_u = circuit_resources(qc_u)
        rows.append({
            "system": "H4",
            "name": spec.name,
            "param": delta,
            "mr_score": ref.mr_score,
            "ansatz": "UCCSD",
            "layers": 0,
            "n_qubits": res_u.n_qubits,
            "n_params": res_u.n_params,
            "depth": res_u.depth,
            "n_2q": res_u.n_2q,
            "E_HF": ref.ehf,
            "E_FCI": ref.e_fci,
            "noon_0": float(ref.noon[0]),
            "noon_1": float(ref.noon[1]),
            "noon_2": float(ref.noon[2]),
            "noon_3": float(ref.noon[3]),
        })

    return pd.DataFrame(rows)

