# Simulation_IVPFISM
# IVPF-ISM vs. FISM Similarity Simulation

This repository contains a Python-based simulation study designed to validate the structural consistency of the proposed **Interval-Valued Picture Fuzzy Interpretive Structural Modeling (IVPF-ISM)** method by comparing it to the classical **Fuzzy ISM (FISM)**.

The comparison is performed through repeated random simulations using **Dice–Sørensen Similarity (DSS)** as the structural similarity metric.

---

## 🎯 Objective

To evaluate whether IVPF-ISM produces hierarchical structures similar to those of the classical FISM method under identical conditions, and to statistically assess the level of structural agreement between the two.

---

## 🧪 Methodology

- Randomly generate decision matrices of varying sizes (5–20 factors) and expert counts (10–30).
- For each matrix:
  - Generate both fuzzy triangular judgments and IVPF values from the same linguistic inputs.
  - Aggregate expert opinions using:
    - Triangular fuzzy averaging for FISM
    - **IVPFOWIA operator** for IVPF-ISM
  - Convert to crisp matrices using respective score functions.
  - Derive Initial Reachability Matrices (IRM) and apply transitive closure for both methods.
- Compute **Dice–Sørensen Similarity (DSS)** between the final binary matrices.

---

## 📈 Output

- A **boxplot** visualization of DSS scores across all simulations.
- A report showing:
  - Mean similarity
  - Standard deviation
  - Number of replications

The final plots are saved as:
- `IVPF_FISM_DSS_Boxplot.png`
- `IVPF_FISM_DSS_Boxplot.pdf`

---

## 🛠️ Requirements

- Python 3.x
- Required packages:
  - `numpy`
  - `matplotlib`

You can install dependencies using:

```bash
pip install numpy matplotlib
