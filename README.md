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

Assessing the Effect of Bias Correction Methods on the Development of Intensity–Duration–Frequency Curves Based on Projections from the CORDEX Central America GCM-RCM Multimodel-Ensemble 

![alt test](/fig01.png)

Graphical Abstract

https://doi.org/10.3390/w16233473

![alt test](/Graphical_Abstract.png)

Abstract: 

This work aims to examine the effect of bias correction (BC) methods on the development of Intensity–Duration–Frequency (IDF) curves under climate change at multiple temporal scales. Daily outputs from a 9-member CORDEX-CA GCM-RCM multi-model ensemble (MME) under RCP 8.5 were used to represent future precipitation. Two stationary BC methods, empirical quantile mapping (EQM) and gamma-pareto quantile mapping (GPM), along with three non-stationary BC methods, detrended quantile mapping (DQM), quantile delta mapping (QDM), and robust quantile mapping (RQM), were selected to adjust daily biases between MME members and observations from the SJO weather station located in Costa Rica. The equidistant quantile-matching (EDQM) temporal disaggregation method was applied to obtain future sub-daily annual maximum precipitation series (AMPs) based on daily projections from the bias-corrected ensemble members. Both historical and future IDF curves were developed based on 5 min temporal resolution AMP series using the Generalized Extreme Value (GEV) distribution. The results indicate that projected future precipitation intensities (2020–2100) vary significantly from historical IDF curves (1970–2020), depending on individual GCM-RCMs, BC methods, durations, and return periods. Regardless of stationarity, the ensemble spread increases steadily with the return period, as uncertainties are further amplified with increasing return periods. Stationary BC methods show a wide variety of trends depending on individual GCM-RCM models, many of which are unrealistic and physically improbable. In contrast, non-stationary BC methods generally show a tendency towards higher precipitation intensities as the return period increases for individual GCM-RCMs, despite differences in the magnitude of changes. Precipitation intensities based on ensemble means are found to increase with the change factor (CF), ranging between 2 and 25% depending on the temporal scale, return period, and non-stationary BC method, with moderately smaller increases for short-durations and long-durations, and slightly higher for mid-durations. In summary, it can be concluded that stationary BC methods underperform compared to non-stationary BC methods. DQM and RQM are the most suitable BC methods for generating future IDF curves, recommending the use of ensemble means over ensemble medians or individual GCM-RCM outcomes.

![alt test](/FIG01.png)

![alt test](/FIG10.png)

![alt test](/FIG11.png)

![alt test](/FIG12.png)

## License

[MIT](https://choosealicense.com/licenses/mit/)
