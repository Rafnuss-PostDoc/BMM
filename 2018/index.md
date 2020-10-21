---
title: Modelling the flow of nocturnal bird migration with year-round European weather radar network.
---

# Modelling the flow of nocturnal bird migration with year-round European weather radar network.
[<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-8185-1020){:target="_blank"}Raphaël Nussbaumer<sup>1,2</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-8182-0152){:target="_blank"}Lionel Benoit<sup>2</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-8820-2808){:target="_blank"}Grégoire Mariethoz<sup>2</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0001-9473-0837){:target="_blank"}Felix Liechti<sup>1</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-0844-164X){:target="_blank"}Silke
Bauer<sup>1</sup> and [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-7736-7527){:target="_blank"}Baptiste Schmid<sup>1</sup>
- <sup>1</sup>Swiss Ornithological Institute, Sempach, Switzerland
- <sup>2</sup>Institute of Earth Surface Dynamics, University of Lausanne, Lausanne, Switzerland

### Links:
[<i class="ai ai-biorxiv"></i> Biorxiv](https://www.biorxiv.org/content/),  [<i class="ai ai-researchgate"></i> Researchgate](https://www.researchgate.net/project/Bird-Migration-Modelling-BMM), [Demo](https://bmm.raphaelnussbaumer.com/2018), [](https://bmm.raphaelnussbaumer.com/2018).



## Dataset selection and pre-processing

1. [`script_download_cleaning.m`](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m)
- [L3-L10](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L3-L10): Download all data available on enram repository for the entire year of 2018.
- [L50-L75](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L50-L75): Delete the data of all radar with aparrent error in them
- [L96-L117](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L96-L117): Reshape data to be on a regular grid: 5min interval, same altitudinal bins.
- [L119-L157](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L119-L157): Compute sunrise and sunset at each radar location
- [L256-L277](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L256-L277): Manually clean the data (use a MATLAB app `CleaningV.mlapp`)
- [L281-L345](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L281-L345): More automatic post-cleaning processing: removal of single data point isolated by 30min. Noise removal (replace value exeeding twice the average of 8 surounding values by the mean of the surounding values), interpolate to 5000m.
- [L449-L498](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L449-L498): Cleaning of flight speed componenet and `sd_vvp` based on the clearning performed on density.
- [L502-L771](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_download_cleaning.m#L502-L771): Insect removal procedure 

 2. [`script_elevation_correction.m`](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_elevation_correction.m)
- [L91-L181](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_elevation_correction.m#L91-L181): MPS simulation for the lower elevation
- [L212-L362](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_elevation_correction.m#L212-L362): Compute volume available for flight
- [L365-L400](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2018/script_elevation_correction.m#L365-L400): Export data for Zenodo

The final interpolated spatio-temporal map can also be downloaded from zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243397.svg)](https://doi.org/10.5281/zenodo.3243397).

## Interpolation
See each livescript for more information.
1. Inference of Bird density [`Density_inference.mlx`](https://rafnuss-postdoc.github.io/BMM/2018/LiveScript/Density_inference)
2. Interpolation of Bird density [`Density_estimationMap.mlx`](https://rafnuss-postdoc.github.io/BMM/2018/LiveScript/Density_estimationMap)
3. Inference of flight speed [`Flight_inference.mlx`](https://rafnuss-postdoc.github.io/BMM/2018/LiveScript/Flight_inference)
4. Interpolation of flight speed [`Flight_estimationMap.mlx`](https://rafnuss-postdoc.github.io/BMM/2018/LiveScript/Flight_estimationMap)

## Flow Model
[`SinkSource.mlx`](https://rafnuss-postdoc.github.io/BMM/2018/LiveScript/SinkSource)

This file contains both the code to compute the fluxes and also produce the figures of the paper. 
