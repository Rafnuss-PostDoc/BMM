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

1. `script_download_cleaning.m`
- Download all data available on enram repository for the entire year of 2018.
- Delete the data of all radar with aparrent error in them
- Compute sunrise and sunset at each radar location
- Manually clean the data (use a MATLAB app `CleaningV.mlapp`)
- More post-cleaning automatic processing (see for detail)

 2. `script_elevation_correction.m`
- MPS simulation for the lower elevation
- Compute volume available for flight
- Export data for Zenodo

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