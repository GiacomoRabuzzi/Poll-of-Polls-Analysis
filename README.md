# Poll-of-Polls-Analysis

This repository contains the code and selected files associated with my MSc thesis on the aggregation of political opinion polls in Italy from 2018 to 2024.  
The project focuses on improving the accuracy, transparency, and interpretability of *poll of polls* estimates using advanced state-space modeling techniques.

## Abstract

This thesis studies the aggregation of political opinion polls in Italy (2018–2024), focusing on improving the accuracy and interpretability of *poll of polls* estimates. Existing approaches often neglect systematic differences across polling institutes and rarely disentangle variance from bias in polling errors.

To address these issues, I develop a family of state-space models that extend the Pelagatti–Monti framework by incorporating:

- Institute-specific **house effects** (bias)
- **Heteroskedastic variances** linked to sample size
- **Binomial measurement structures**
- A **covariance matrix** capturing structural dependencies among polling institutes

Model estimation relies on **Kalman filtering**, **Durbin–Koopman approximations**, and **importance sampling**.

The results show that variance is the main driver of discrepancies across polls, while institute-specific biases are generally modest (around ±1 percentage point) and party-specific. A small number of institutes display persistent deviations, whereas others differ mainly due to higher sampling variance. Structural correlations suggest limited clusters of methodological affinity rather than broad polarization. Smoothed vote-share trends closely align with well-known Italian electoral dynamics.

Overall, the thesis demonstrates that effective poll aggregation requires bias correction, variance-based weighting, and structural diagnostics, and it proposes a replicable methodology with applications beyond the Italian case.


## Methodology Overview

The core contribution of this project is a flexible state-space framework for poll aggregation that:

- Separates **bias** from **sampling variance**
- Accounts for **institute heterogeneity**
- Allows for **correlated polling errors**
- Produces **smoothed latent vote-share trajectories**

The models are estimated using simulation-based methods suitable for non-Gaussian and high-dimensional settings.

## Requirements

The code is written primarily in **R** (with possible extensions in other languages).  
Required packages are listed within the scripts and should be installed before running the analysis.

## Reproducibility

The repository is designed to be as reproducible as possible given data-sharing constraints.  
Where raw polling data cannot be redistributed, instructions are provided to reconstruct the datasets from public sources.

## License

This repository is intended for **academic and research use**.  
Please cite the thesis if you use or adapt the code.

## Author

Giacomo Rabuzzi  
Università degli Studi di Milano Bicocca
2026


The repository is organized as follows:

