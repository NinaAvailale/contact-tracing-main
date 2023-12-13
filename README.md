# Contact-Tracing-Main
Code for: Revisiting the effectiveness of digital contact tracing for cospace-time interactions



## Citation
Revisiting the effectiveness of digital contact tracing for cospace-time interactions



## Abstract
While many countries employed digital contact tracing to contain the spread of SARS-CoV-2, the effectiveness of tracing numerous individuals involved in cospace-time interaction (i.e., individuals who shared the same space and time) in real world has seldom been revisited due to the lack of contacts who tested negative. To address this issue, we utilized data from 2,230 infections and 220,878 contacts with detailed epidemiological information during the Omicron outbreak in Beijing in 2022. We observed that daily contact numbers for individuals in dwelling, workplace, cospace-time interaction, and community settings could be described by gamma distribution with distinct parameters. Our findings revealed that 37.99% of infections occurred through cospace-time interactions with control measures. However, using a mathematical model to incorporate contacts under different locations, we found that cospace-time interactions contributed to only 11.09% of transmission and 3.85% of super-spreading risk without control measures, the lowest among all locations studied. These results suggest that public health measures should be adapted to strike a balance between the benefits of digital contact tracing for cospace-time interactions and the involvement of a substantial number of individuals in this location.

## System requirements
### Software Version Information:
This code executes COVID-19 mathematical analyses using R version 4.2.1. Specifically, it addresses the infectiousness model p(tau) for various locations with input parameters provided. The code also calculates contributions to transmission and the probability of super spreading events (SSE) across different locations.

The analyses take approximately 3 seconds to run on a machine that meets the recommended specifications.


### R dependencies
Users should install the following packages prior to performing the COVID-19 mathematical analysis, from an R terminal:
```
install.packages(c('ggplot2', 'tidyverse', 'rriskDistributions', 'gridExtra', 'RColorBrewer', 'R.utils'))
```

## Notes on the code
 

About the infectiousness model: Our infectiousness model is a modified version based on Ferretti, Wymant et al., Science 2020. If you use or refer to this code, kindly cite the aforementioned publication:

Ferretti, L. et al. Quantifying SARS-CoV-2 transmission suggests epidemic control with digital contact tracing. Science 368, eabb6936 (2020).

About the code folder: We used the "model" folder to estimate the infectiousness, i.e., p(tau), under different locations for specified input parameters and to estimate the contributions to transmission. The "tranmission_and_SSE_risk_across_locations" folder was used to calculate the probability of super-spreading events (SSE) under different locations. Detailed information on the contact tracing and epidemiological data of the Omicron outbreak in the spring of 2022 in Beijing, relevant to the input data, is stored in the "data_processed" folder. For configuration and computational methodology for parameter settings, please refer to the supplementary materials.

