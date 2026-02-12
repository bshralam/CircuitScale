#Classical references + NOONs + MR score (the backbone)
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Any, Tuple, Optional
import numpy as np
from pyscf import gto, scf, mcscf, fci

from .systems import MoleculeSpec

@dataclass
class ReferenceResult:
    ehf: float
    e_fci: Optional[float]
    noon: Optional[np.ndarray]
    mr_score: Optional[float]
    meta: Dict[str, Any]

def build_mol(spec: MoleculeSpec) -> gto.Mole:
    mol = gto.Mole()
    mol.atom = [(sym, *xyz) for sym, xyz in spec.geometry]
    mol.basis = spec.basis
    mol.charge = spec.charge
    mol.spin = spec.spin
    mol.unit = spec.unit
    mol.build()
    return mol

def run_rhf(mol: gto.Mole) -> scf.hf.RHF:
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10
    mf.kernel()
    return mf

def run_fci_energy(mf: scf.hf.RHF) -> float:
    # For tiny systems (H2/H4 STO-3G), FCI is feasible.
    cisolver = fci.FCI(mf.mol, mf.mo_coeff)
    e, _ = cisolver.kernel()
    return float(e)

def run_casscf_noons(
    mf: scf.hf.RHF,
    ncas: int,
    nelecas: Tuple[int, int] | int,
) -> Tuple[np.ndarray, float]:
    """
    Returns (NOONs, E_CASSCF). NOONs are eigenvalues of 1-RDM in active space basis.
    """
    mc = mcscf.CASSCF(mf, ncas=ncas, nelecas=nelecas)
    mc.conv_tol = 1e-10
    mc.kernel()

    # 1-RDM in MO basis (full), then take active block.
    rdm1 = mc.make_rdm1()
    # Active orbitals indices in MO basis: mc.ncore ... mc.ncore+mc.ncas-1
    a0 = mc.ncore
    a1 = mc.ncore + mc.ncas
    rdm1_act = rdm1[a0:a1, a0:a1]
    noons = np.linalg.eigvalsh(rdm1_act)[::-1]  # descending
    return noons, float(mc.e_tot)

def mr_score_from_noons(noons: np.ndarray) -> float:
    """
    Simple, interpretable MR score:
    sum_i min(n_i, 2 - n_i) over active orbitals (closed-shell convention).
    - If occupations are (2,0,2,0...) => score ~ 0 (single reference)
    - If occupations are ~1,1,... => score grows (strong correlation)
    """
    noons = np.array(noons, dtype=float)
    return float(np.sum(np.minimum(noons, 2.0 - noons)))

def compute_reference(
    spec: MoleculeSpec,
    do_fci: bool = True,
    do_casscf: bool = True,
    ncas: int = 2,
    nelecas: Tuple[int, int] | int = (1, 1),
) -> ReferenceResult:
    mol = build_mol(spec)
    mf = run_rhf(mol)
    ehf = float(mf.e_tot)

    e_fci = None
    if do_fci:
        e_fci = run_fci_energy(mf)

    noon = None
    mr_score = None
    if do_casscf:
        noon, e_cas = run_casscf_noons(mf, ncas=ncas, nelecas=nelecas)
        mr_score = mr_score_from_noons(noon)
        meta = {"e_casscf": e_cas, "ncas": ncas, "nelecas": nelecas}
    else:
        meta = {}

    meta.update({"nao": mol.nao_nr(), "nelec": mol.nelectron})
    return ReferenceResult(ehf=ehf, e_fci=e_fci, noon=noon, mr_score=mr_score, meta=meta)


