# Publish and link this project

1. Sign in to [GitHub](https://github.com/) and open [Create a new repository](https://github.com/new).
2. **Repository name:** `ae370-falcon9-fin-stabilization`. Select **Public** (or give your instructor access if private). Optional description: `AE 370 Group Project 1 — Falcon 9-inspired movable-fin attitude stabilization using a custom RK4 solver`.
3. **Do not** initialize it with a new README, .gitignore, or license: this upload already contains README and .gitignore.
4. Create the repository. On its empty-repository page click **uploading an existing file**, or navigate to **Add file → Upload files**.
5. **Unzip** the package first. Drag the *contents* of the `ae370-falcon9-fin-stabilization` folder (README.md, requirements.txt, folders, etc.) into the GitHub upload area. Do **not** upload only the ZIP archive; the instructor should be able to browse code online.
6. Choose a commit message such as `Add AE 370 Falcon 9 numerical methods project` and click **Commit changes**.
7. Open the code notebook on GitHub and confirm you can see the code cells. In the root README, check that the six output figures and scripts are accessible.
8. Copy your repository URL (in the form shown below) into your LaTeX report.

If GitHub's web upload does not accept your folder structure, upload it in smaller batches or use the GitHub Desktop application to publish the unzipped folder.

## Paste into LaTeX (Overleaf)

In your document preamble, include `\usepackage{hyperref}` if it is not already loaded. Then add this near the Reproducibility section:

```latex
\section{Reproducibility}
The full Python implementation, including the custom Runge--Kutta solver,
independent verification tests, convergence studies, parameter sweep, and
code used to generate all reported figures, is available at
\url{https://github.com/YOUR-USERNAME/ae370-falcon9-fin-stabilization}.
The repository contains a runnable Jupyter notebook, separated source modules,
and instructions for reproducing the results.
```

Replace `YOUR-USERNAME` with the exact account or organization name you used.

**Important:** Replace any old code-link placeholder in the report. Submit the PDF generated from your LaTeX file, as required by the assignment; merely uploading a notebook does not satisfy that part of the instructions.
