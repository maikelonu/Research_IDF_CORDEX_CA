# Research_IDF_CORDEX_CA 
# IDF_CORDEX

HBV_OPI is an R program designed to perform optimization and sensitivity analysis of hydrological models. To do so, it uses several local optimization methods, including: (1) Nelder-Mead, (2) Broyden-Fletcher-Goldfarb-Shanno, (3) Hooke-Jeeves, (4) Variable Nonlinear Minimization, (5) Bound Optimization BY Quadratic Approximation, (6) Spectral Projected Gradient, (7) PORT Gradient Algorithm and (8) Levenberg-Marquardt Algorithm, as well as global optimization methods, including: (1) Generalized Simulated Annealing, (2) Differential Evolution Optimization, (3) Genetic Algorithms, (4) Shuffled Complex Evolution, (5) Enhanced Particle Swarm Optimization, (6) DIviding RECTangles for Global Optimization, (7) Controlled Random Search Local Mutation and (8) Augmented Lagrangian Minimization Algorithm.

## Installation

Use Git to clone and install the program

## Usage

Run HBV_OPTI.R script accordingly

## Contributing

Maikel Mendez Morales. Escuela de Ingeniería en Construcción, Instituto Tecnológico de Costa Rica. email: mamendez@itcr.ac.cr

## Publications

Mendez, M.; Calvo-Valverde, L.A. Comparison of global and local optimization methods for the calibration and sensitivity analysis of a conceptual hydrological model. TM 2019.
https://doi.org/10.18845/tm.v32i3.4477

![alt test](/edp01.png)

Graphical Abstract

![alt test](/FIG_TM_02.png)

Abstract: 

This work aims to examine the effect of bias correction (BC) methods on the development of Intensity-Duration-Frequency (IDF) curves under climate change at multiple temporal scales.     Daily outputs from a 9-member CORDEX-CA GCM-RCM multi-model ensemble (MME) under RCP 8.5 were used to represent future precipitation. Two stationary BC methods; empirical quantile mapping (EQM) and gamma-pareto quantile mapping (GPM), along with three non-stationary BC methods; detrended quantile mapping (DQM), quantile delta mapping (QDM) and robust quantile mapping (RQM), were selected to adjust daily biases between MME members and observations from the SJO weather station located in Costa Rica. The equidistant quantile-matching (EDQM) temporal disaggregation method was applied to obtain future sub-daily annual maximum precipitation series (AMPs) based on daily projections from the bias-corrected ensemble members. Both historical and future IDF curves were developed based on 5-min temporal resolution AMPs series using the Generalized Extreme Value (GEV) distribution. Results indicate that projected future precipitation intensities (2020-2100) vary significantly from historical IDF curves (1970-2020), depending on individual GCM-RCMs, BC methods, durations, and return periods. Regardless of stationarity, the ensemble-spread increases steadily with return period, as uncertainties are further amplified with increasing return periods. Stationary BC methods show a wide variety of trends depending on individual GCM-RCM models, many of which are unrealistic and physically improbable. In contrast, non-stationary BC methods generally show a tendency towards higher precipitation intensities as the return period increases for individual GCM-RCMs, despite differences in the magnitude of changes. Precipitation intensities based on ensemble-means are found to increase with a change-factor (CF) ranging between 2 and 25% depending on temporal scale, return period and non-stationary BC method, with moderately

![alt test](/Opti02.png)

![alt test](/Opti03.png)

![alt test](/Opti04.png)

## License

[MIT](https://choosealicense.com/licenses/mit/)
