#Runner: sweep geometry → CSV (this is what makes it feel “real”)
from __future__ import annotations
from typing import List, Dict, Any
import pandas as pd

from .systems import MoleculeSpec
from .references import compute_reference
from .hamiltonians import molecular_qubit_hamiltonian
from .ansatze import hea, uccsd_placeholder
from .resources import circuit_resources

def run_sweep_h2(specs: List[MoleculeSpec], hea_layers=(1,2,3)) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    for spec in specs:
        # ----- Reference / MR metric
        ref = compute_reference(
            spec,
            do_fci=True,
            do_casscf=True,
            ncas=2,
            nelecas=(1, 1),
        )

        # ----- Hamiltonian / qubits
        qham, n_qubits = molecular_qubit_hamiltonian(spec)

        # ----- HEA resources at different depths
        for L in hea_layers:
            qc = hea(n_qubits, layers=L, entangle="linear")
            res = circuit_resources(qc)
            rows.append({
                "system": "H2",
                "name": spec.name,
                "param": float(spec.name.split("R")[1]),  # quick parse
                "mr_score": ref.mr_score,
                "noon_0": float(ref.noon[0]),
                "noon_1": float(ref.noon[1]),
                "ansatz": "HEA",
                "layers": L,
                "n_qubits": res.n_qubits,
                "n_params": res.n_params,
                "depth": res.depth,
                "n_2q": res.n_2q,
                "E_HF": ref.ehf,
                "E_FCI": ref.e_fci,
            })

        # ----- UCCSD placeholder (swap later)
        qc_u = uccsd_placeholder(n_qubits)
        res_u = circuit_resources(qc_u)
        rows.append({
            "system": "H2",
            "name": spec.name,
            "param": float(spec.name.split("R")[1]),
            "mr_score": ref.mr_score,
            "noon_0": float(ref.noon[0]),
            "noon_1": float(ref.noon[1]),
            "ansatz": "UCCSD",
            "layers": 0,
            "n_qubits": res_u.n_qubits,
            "n_params": res_u.n_params,
            "depth": res_u.depth,
            "n_2q": res_u.n_2q,
            "E_HF": ref.ehf,
            "E_FCI": ref.e_fci,
        })

    return pd.DataFrame(rows)

