# Research_IDF_CORDEX_CA 
# IDF_CORDEX

IDF_CORDEX_CA is an R program designed to generate Intensity-Duration-Frequency (IDF) curves under Climate Change at multiple temporal scales using precipitation projections from the Coordinated Regional Climate Downscaling Experiment (CORDEX) for the Central America domain (CA) (https://esg-dn1.nsc.liu.se/projects/cordex/). IDF_CORDEX_CA covers a 9-member multi-model ensemble (MME) combining various GCMs and RCMs with a spatial resolution of 0.22° x 0.22° (~25 km). To do so, it uses several bias correction (BC) techniques, including two stationary BC methods; empirical quantile mapping (EQM) and gamma-pareto quantile mapping (GPM), along with three non-stationary BC methods; detrended quantile mapping (DQM), quantile delta mapping (QDM) and robust quantile mapping (RQM). IDF_CORDEX_CA uses the equidistant quantile-matching (EDQM) temporal disaggregation to obtain future sub-daily annual maximum precipitation series (AMPs) based on daily projections from the bias-corrected ensemble members. Extreme Value Analysis (EVA) can be executed using  (1) Generalized Extreme Value (GEV), (2) 2-parameter Gumbel (EV1), or (3) 3-parameter Log-Pearson type 3 (LP3).

## Installation

Use Git to clone and install the program

## Usage

Run sources.R script accordingly

## Contributing

Maikel Mendez Morales. Escuela de Ingeniería en Construcción, Instituto Tecnológico de Costa Rica. email: mamendez@itcr.ac.cr

## Publications

To be defined (in press)

![alt test](/fig01.png)

Graphical Abstract

To be defined (in press)

![alt test](/fig02.png)

Abstract: 

To be defined (in press)

![alt test](/fig03.png)

![alt test](/fig04.png)

![alt test](/fig05.png)

## License

[MIT](https://choosealicense.com/licenses/mit/)
