# Complex Networks | projects

In this repository you will find the simulations and analyses carried out for the final exam of the *Physics of Complex Networks: Structure and Dynamics* course in the *Physics of Data* master’s programme (University of Padova, AY 2023/24).

The repo contains the code, processed data and report for two selected tasks:

- **Voter model (#28)**  
  Theoretical and numerical study of the classical voter model, the noisy voter model and the q-voter model on various network topologies to explore ordering, consensus time and finite-size effects.  
  **Score:** 0.5

- **European transportation networks I (#45)**  
  Construction and analysis of multi-layer European transport networks (rail, air and road); computation of centrality measures and community structure, followed by comparative discussion.  
  **Score:** 1.0

---

## About reproducibility

The simulations depend on the C++/Python library **graph-tool**, which is easiest to install via *conda-forge* (PyPI wheels are not available).

```bash
# Recommended environment
conda create -n pocn python=3.11 graph-tool jupyterlab numpy scipy matplotlib
conda activate pocn
