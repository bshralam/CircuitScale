from __future__ import annotations
from typing import List, Dict, Any, Tuple
import re
import pandas as pd

from .systems import MoleculeSpec
from .references import compute_reference
from .hamiltonians import molecular_qubit_hamiltonian
#from .ansatze import hea, uccsd_placeholder
from .ansatze import hea_with_hf_init, uccsd_placeholder
from .vqe_utils import qubitop_to_sparsepauliop, run_vqe
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
    hea_layers=(1, 2, 3, 4, 5, 6),
    target_error: float = 1.6e-3,
    maxiter: int = 300,
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

        # ----- Qubit Hamiltonian for VQE
        op = qubitop_to_sparsepauliop(qham, n_qubits)

        # H4 has 4 electrons in this setup
        n_electrons = 4

        # ----- HEA resources + VQE error at different depths
        best_layer = None
        best_depth = None
        best_n_2q = None
        best_energy = None
        best_error = None

        for L in hea_layers:
            qc = hea_with_hf_init(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                layers=L,
                entangle="linear",
            )
            res = circuit_resources(qc)

            vqe_res = run_vqe(op, qc, maxiter=maxiter, seed=7)
            err = abs(vqe_res.energy - ref.e_fci)

            rows.append({
                "system": "H4",
                "name": spec.name,
                "param": delta,
                "mr_score": ref.mr_score,
                "ansatz": "HEA",
                "layers": L,
                "n_qubits": res.n_qubits,
                "n_params": res.n_params,
                "depth": res.depth,
                "n_2q": res.n_2q,
                "E_HF": ref.ehf,
                "E_FCI": ref.e_fci,
                "E_VQE": vqe_res.energy,
                "vqe_error": err,
                "target_error": target_error,
                "meets_target": err <= target_error,
                "noon_0": float(ref.noon[0]),
                "noon_1": float(ref.noon[1]),
                "noon_2": float(ref.noon[2]),
                "noon_3": float(ref.noon[3]),
            })

            if best_layer is None and err <= target_error:
                best_layer = L
                best_depth = res.depth
                best_n_2q = res.n_2q
                best_energy = vqe_res.energy
                best_error = err

        rows.append({
            "system": "H4",
            "name": spec.name,
            "param": delta,
            "mr_score": ref.mr_score,
            "ansatz": "HEA_required",
            "layers": best_layer,
            "n_qubits": n_qubits,
            "n_params": None,
            "depth": best_depth,
            "n_2q": best_n_2q,
            "E_HF": ref.ehf,
            "E_FCI": ref.e_fci,
            "E_VQE": best_energy,
            "vqe_error": best_error,
            "target_error": target_error,
            "meets_target": best_layer is not None,
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

