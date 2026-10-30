# Movable-Fin Attitude Stabilization of a Falcon 9-Inspired Rocket

**AE 370 — Aerospace Numerical Methods | Group Project 1 | Fall 2026**

**Group:** Bharat Rajesh, Thomas Theodor, Jorgen Nelson, Carlos Avila

This project models the rotational attitude recovery of a Falcon 9-inspired rocket equipped with four movable fins. We are investigating how **fin deflection limits** and **actuator response time** affect both recovery time and fin effort. The physical parameters and control model are simplified assumptions for a numerical-methods study, **not** a validated or flight-ready Falcon 9 model.

The ten-state nonlinear initial value problem is advanced using a **self-written classical fourth-order Runge–Kutta (RK4)** integrator, with forward Euler used for comparison. **No external ODE integrator is used.**

## Repository contents

| Location | Purpose |
|---|---|
| `notebooks/falcon9_fin_stabilization.ipynb` | **Complete original executable notebook:** model, method, verification checks, convergence studies, figures, 49-case parameter sweep, sensitivity analyses, and report text. This is the primary reproduction artifact. |
| `src/model.py` | Physical and actuator model, initial state and parameters. |
| `src/integrators.py` | From-scratch Euler and RK4 integration plus rotation-matrix utility. |
| `src/stability.py` | Linearization, RK4 stability polynomial and local time-step limit. |
| `src/metrics.py` | Attitude recovery and fin-effort calculations. |
| `scripts/run_verification.py` | Independent checks, stability analysis, convergence and cost studies; creates Figures 1–3. |
| `scripts/run_results.py` | Nominal response, 49-case fin/actuator sweep, sensitivity study; creates Figures 4–6. |
| `figures/` | Images generated from the code. |
| `requirements.txt` | Python dependencies. |

The `src/` and `scripts/` files reorganize the original notebook code so that the numerical method can be inspected independently from the verification tests and the results. The code in the notebook has been retained.

## Running the notebook (simplest method)

1. Open `notebooks/falcon9_fin_stabilization.ipynb` in **Google Colab** or **Jupyter**.
2. Make sure **NumPy** and **Matplotlib** are installed (both are available in common Colab environments).
3. Select **Runtime → Run all** in Colab, or **Kernel → Restart & Run All** in Jupyter.
4. The notebook generates numerical verification outputs, all six figures, convergence/error and cost results, and the full parameter sweep. Run the cells top-to-bottom; later cells depend on definitions in earlier cells.

## Running the separate Python scripts

Run these commands **from the repository root**:

```bash
python -m pip install -r requirements.txt
python -m scripts.run_verification
python -m scripts.run_results
```

The scripts print the numerical results and save `figures/figure-01.png` through `figures/figure-06.png`. Depending on your system, use `python3` instead of `python`.

The verification script includes: comparison to exact special-case solutions; conservation tests for a torque-free body; finite-difference confirmation of the linearized Jacobian; RK4 absolute-stability and convergence studies; and comparisons against a refined time-step reference. The results script computes nominal settling time, zero-command comparison, the 49-setting fin-limit/actuator-lag study with refinement checks, and sensitivity to assumed force/inertia.

**Expected nominal settling time:** about **7.27 s** for the model and chosen initial conditions. Results are model predictions, not measured Falcon 9 flight performance.

## Report and assignment submission

The **course assignment requires a separately submitted LaTeX-typeset PDF report**. This repository provides the source code and figures to be linked in that report; a notebook alone does not replace the PDF. In the report's Reproducibility section, link to this GitHub repository (or to the full notebook within it). Make sure the repository is public or accessible to the instructor.

> Before submitting, the group should confirm the actual member contribution percentages. The uploaded notebook and report-only version contain different provisional contribution allocations. Also replace any `REPLACE_WITH_YOUR_ACCESSIBLE_NOTEBOOK_OR_REPOSITORY_LINK` placeholder in the submitted report.

## Acknowledgments

The accompanying project notebook describes substantive use of ChatGPT (OpenAI) for ideation, modelling, implementation, tests, figures, and report preparation. Refer to the notebook's Acknowledgments section for the group's disclosure and retain an accurate disclosure in the submitted report.
