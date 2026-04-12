# CircuitScale: Multireference Electronic Structure and VQE Resource Scaling

## Objective

The goal of this project is to quantify how multireference (MR) character in molecular electronic structure impacts:

- VQE accuracy (energy error relative to FCI)
- optimization behavior
- quantum circuit resource requirements (depth, two-qubit gates, ansatz layers)

We use minimal hydrogen systems (H₂ and H₄) as controlled testbeds to isolate the relationship between:
\[
\text{electronic structure complexity} \;\longleftrightarrow\; \text{quantum algorithm difficulty}
\]

---

## Theoretical Background

### Multireference Diagnostic

We quantify correlation using natural orbital occupation numbers (NOONs):

\[
MR = \sum_i \min(n_i, 2 - n_i)
\]

- \( n_i \approx 2, 0 \): single-reference (weak correlation)
- \( n_i \approx 1 \): strong multireference (near-degeneracy)

This metric captures deviation from idempotency of the one-particle density matrix.

---

### Systems

#### H₂ (bond stretching)
- Correlation increases **monotonically** with bond length
- Classic single → multireference transition

#### H₄ (rectangular distortion, parameter δ)
- Correlation is **non-monotonic**
- Maximum MR occurs at **δ = 0 (symmetric geometry)** due to orbital degeneracy

---

## Results and Figure-by-Figure Analysis

---

### 1. H₂: NOONs vs bond distance

- NOON₁ decreases from ~2 → ~1
- NOON₂ increases from ~0 → ~1

This reflects gradual emergence of static correlation as the bond is stretched.

**Key point:**  
H₂ exhibits a *smooth, monotonic increase* in multireference character.

---

### 2. H₂: MR score vs bond distance

- MR increases with bond length
- Peaks at large R (dissociation limit)

**Interpretation:**  
Correlation is driven purely by bond stretching → no symmetry-induced effects.

---

### 3. H₂: VQE error vs MR score

- Error increases with MR
- Strong correlation → larger deviation from FCI

**Observation:**
- Behavior is relatively smooth
- Optimization landscape remains manageable

---

### 4. H₂: Required HEA layers / depth / 2Q gates vs MR

Two regimes:

#### (a) Fixed ansatz (fixed L)
- Circuit depth and gate counts are constant
- Independent of MR

#### (b) Minimum required ansatz
- Required layers increase with MR
- Reflects increased expressibility needed

**Key point:**  
MR directly correlates with required circuit expressibility.

---

## Contrast: H₂ vs H₄

| Property | H₂ | H₄ |
|------|----|----|
| Correlation driver | Bond stretching | Orbital degeneracy (symmetry) |
| MR behavior | Monotonic | Non-monotonic (peaks at δ = 0) |
| Electronic structure | Smooth evolution | Abrupt changes near symmetry |
| VQE behavior | Smooth scaling | Highly irregular |

---

## 5. H₄: MR score vs distortion δ

- MR peaks sharply at δ = 0
- Drops away from symmetry

**Interpretation:**
- Symmetric geometry → near-degenerate orbitals
- Leads to maximal multireference character

---

### 6. H₄: NOONs vs distortion δ

At δ = 0:
- Multiple orbitals have occupation ~1
- Strong configuration mixing

Away from δ = 0:
- Occupations return toward ~2 / ~0
- System becomes single-reference

**Key point:**  
Correlation in H₄ is **symmetry-driven**, not geometry-driven.

---

### 7. H₄: Best achievable VQE error vs distortion

- Error is largest near δ = 0
- Decreases away from symmetry

**Interpretation:**
- Strong correlation → difficult optimization landscape
- HEA struggles to represent multi-configurational wavefunctions

---

### 8. H₄: VQE error vs MR score

- Non-monotonic relationship
- Large scatter in error for similar MR values

**Key observation:**
- MR alone does not fully determine difficulty
- Optimization landscape effects are significant

---

### 9. H₄: Required HEA layers vs MR

- Layers increase near high MR region
- Sharp transitions (non-smooth scaling)

**Interpretation:**
- Required expressibility increases abruptly near degeneracy
- Indicates sensitivity to electronic structure topology

---

### 10. H₄: Fixed HEA (layers vs gates/depth)

- Number of gates and depth are constant for fixed L
- Independent of MR

**Important clarification:**
These plots reflect **ansatz structure**, not problem difficulty.

---

## Key Technical Insights

1. **Two distinct sources of correlation**
   - H₂: bond stretching (continuous)
   - H₄: symmetry-induced degeneracy (discontinuous)

2. **MR is necessary but not sufficient**
   - Captures correlation strength
   - Does not capture optimization difficulty

3. **Symmetry points are the hardest cases**
   - Maximal degeneracy
   - Flat/complex energy landscapes
   - Largest VQE errors

4. **HEA limitations**
   - Fixed structure → does not adapt to correlation
   - Inefficient for multireference systems

5. **Resource scaling depends on perspective**
   - Fixed ansatz → constant cost
   - Required ansatz → grows with MR

---

## Takeaway

This study shows that:

\[
\text{Electronic structure (NOONs, degeneracy)} 
\;\Rightarrow\; 
\text{MR character}
\;\Rightarrow\;
\text{VQE difficulty + resource requirements}
\]

However, **symmetry-induced multireference effects (H₄)** introduce additional complexity beyond simple MR scaling, making them a more stringent test for quantum algorithms than bond stretching systems like H₂.

---

## Tools

- PySCF (FCI, NOONs)
- Qiskit / PennyLane (VQE)
- NumPy / Matplotlib
