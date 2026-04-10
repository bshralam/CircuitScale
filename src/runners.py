#Runner: sweep geometry → CSV (this is what makes it feel “real”)
from __future__ import annotations
from typing import List, Dict, Any
import pandas as pd

from .systems import MoleculeSpec
from .references import compute_reference
from .hamiltonians import molecular_qubit_hamiltonian
#from .ansatze import hea, uccsd_placeholder
from .ansatze import hea_with_hf_init, uccsd_placeholder
from .vqe_utils import qubitop_to_sparsepauliop, run_vqe
from .resources import circuit_resources

#def run_sweep_h2(specs: List[MoleculeSpec], hea_layers=(1,2,3)) -> pd.DataFrame:
def run_sweep_h2(specs: List[MoleculeSpec],hea_layers=(1, 2, 3, 4, 5, 6),target_error: float = 1.6e-3,maxiter: int = 200) -> pd.DataFrame: 
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

        """# ----- HEA resources at different depths
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
            })"""

                # ----- Qubit Hamiltonian for VQE
        op = qubitop_to_sparsepauliop(qham, n_qubits)

        # H2 has 2 electrons
        n_electrons = 2

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
                "system": "H2",
                "name": spec.name,
                "param": float(spec.name.split("R")[1]),
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
                "E_VQE": vqe_res.energy,
                "vqe_error": err,
                "target_error": target_error,
                "meets_target": err <= target_error,
            })

            if best_layer is None and err <= target_error:
                best_layer = L
                best_depth = res.depth
                best_n_2q = res.n_2q
                best_energy = vqe_res.energy
                best_error = err

        rows.append({
            "system": "H2",
            "name": spec.name,
            "param": float(spec.name.split("R")[1]),
            "mr_score": ref.mr_score,
            "noon_0": float(ref.noon[0]),
            "noon_1": float(ref.noon[1]),
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

