---
title: 2016
---

# A geostatistical approach to estimate high resolution nocturnal bird migration densities from a weather radar network
[<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-8185-1020){:target="_blank"}Raphaël Nussbaumer<sup>1,2</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-8182-0152){:target="_blank"}Lionel Benoit<sup>2</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-8820-2808){:target="_blank"}Grégoire Mariethoz<sup>2</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0001-9473-0837){:target="_blank"}Felix Liechti<sup>1</sup> , [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-0844-164X){:target="_blank"}Silke
Bauer<sup>1</sup> and [<i class="ai ai-orcid"></i>](https://orcid.org/0000-0002-7736-7527){:target="_blank"}Baptiste Schmid<sup>1</sup>
- <sup>1</sup>Swiss Ornithological Institute, Sempach, Switzerland
- <sup>2</sup>Institute of Earth Surface Dynamics, University of Lausanne, Lausanne, Switzerland

### Links:
[<i class="ai ai-biorxiv"></i> Biorxiv](https://www.biorxiv.org/content/10.1101/690065), [<i class="ai ai-doi"></i> Publication in Remote sensing](https://doi.org/10.3390/rs11192233),  [<i class="ai ai-researchgate"></i> Researchgate](https://www.researchgate.net/project/Bird-Migration-Modelling-BMM), [Demo](https://bmm.raphaelnussbaumer.com/2016), 
[Presentation](https://docs.google.com/viewer?url=https://www.researchgate.net/profile/Raphael_Nussbaumer/publication/332028742_Space-time_interpolation_of_nocturnal_bird_migration/links/5c9b85cda6fdccd4603f1120/Space-time-interpolation-of-nocturnal-bird-migration.pdf).

## Abstract
1. Quantifying nocturnal bird migration at high resolution is essential for (1) understanding the phenology of migration and its drivers, (2) identifying critical spatio-temporal protection zones for migratory birds, and (3) assessing the risk of collision with man-made structures.
2. We propose a tailored geostatistical model to interpolate migration intensity monitored by a network of weather radars. The model is applied to data collected in autumn 2016 from 69 European weather radars. To validate the model, we performed a cross-validation and also compared our interpolation results with independent measurements of two bird radars.
3. Our model estimated bird densities at high resolution (0.2° latitude-longitude, 15min) and assessed the associated uncertainty. Within the area covered by the radar network, we estimated that around 120 million birds were simultaneously in flight (10-90 quantiles: 107-134). Local estimations can be easily visualized and retrieved from a dedicated interactive website: [birdmigrationmap.vogelwarte.ch](https://birdmigrationmap.vogelwarte.ch/).
4. This proof-of-concept study demonstrates that a network of weather radar is able to quantify bird migration at high resolution and accuracy. The model presented has the ability to monitor population of migratory birds at scales ranging from regional to continental in space and daily to yearly in time. Near-real-time estimation should soon be possible with an update of the infrastructure and processing software.

## Dataset
Our dataset originates from measurements of 69 European weather radars covering a large part of the Western-European flyway during fall migration 2016.
The raw bird density data used in this study are found on the repository of [European Network for the Radar surveillance of Animal Movement (ENRAM)](http://enram.github.io/data-repository/) and were generated with [vol2bird](https://github.com/adokter/vol2bird).
These data were cleaned manually into vertical profile of reflectivity (see Appendix A in paper for details). These data are available on zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243397.svg)](https://doi.org/10.5281/zenodo.3243397)

<img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/10-paper/figure/Figure2.png">
<span style="font-size:0.8em;">Illustration of the spatio-temporal variability of bird densities measured by a weather radar network. (a) Average bird densities measured by each radar over the whole study period. (b-d) Time series of bird densities measured by the radar with the corresponding outer ring colour in panel a. (e-g) Zoom on a two-days period.
A strong continental trend appears in (a) as well as a correlation at the multi-night scale when comparing (b), (c) and (d). These spatial correlations are even stronger at the regional scale when comparing within a subplot (e.g. (f)).The intra-night scale shows an obvious bell-shape curve pattern during each night (e.g. (g)).</span>


## Model description

The strong spatio-temporal correlations  of bird migration motivated the use of a Gaussian process regression to interpolate bird densities measured by weather radars at high temporal resolution. Because of the multi-scale temporal structure of bird migration, we consider here an additive model combining two temporal scales: first, a multi-night process that models bird density averaged over the night, and second, an intra-nightly process that models variations within each night. Subsequently, each scale-specific process is further split into two terms: a smoothly-varying (in space and/or time) deterministic trend, and a stationary Gaussian process.

The bird density $$B(\mathbf{s},t)$$ observed at location $$\mathbf{s}$$ and at time $$t$$ is modeled by

$$B( \mathbf{s} ,t)^p = \mu( \mathbf{s} ) + M(\mathbf{s},d(t)) + \iota (\mathbf{s},t) + I(\mathbf{s},t)$$

where $$\mu$$ and $$\iota$$ are deterministic trends (respectively at the multi-night and intra-night resolution), and $$M$$ and $$I$$ are random effects (respectively at the multi-night and intra-night resolution) ; $$d(t)$$ is a step function that maps the continuous time $$t$$ to the discrete day $$d$$ of the closest night. A power transformation $$p$$ is applied on bird densities to transform the highly skewed marginal distribution into a Gaussian distribution. 


<img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/10-paper/figure/Figure3.png">
<span style="font-size:0.8em;">Illustration of the proposed mathematical model decomposition of Eqn 1 with the exception that the power transformation was not applied. Note that the values of $$M$$,$$\iota$$,$$I$$ can be either positive or negative.</span>

| Spatial trend  | Multi-night | Curve trend  | Intra-night |
| ------------- | ------------- | ------------- | ------------- |
|  <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/trend.png"> | <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/Density_estimationMap_amplitude.gif">  | <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/curve.png">  | <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/Density_estimationMap_residu.gif">  |

More details on the model are available in the manuscript or through the inference script (see below).

## Inference
The inference of the model parameters is performed with this [MATLAB LiveScript](https://rafnuss-postdoc.github.io/BMM/2016/LiveScript/Inference.html).

## Validation
1. Cross-validation for each radar by ignoring the data of this radar and estimating the bird density at the same location and time [MATLAB LiveScript of cross-validation](https://rafnuss-postdoc.github.io/BMM/2016/LiveScript/Cross_validation.html)
2. Comparison with Birdscan radars [MATLAB LiveScript of validation with bird radar](https://rafnuss-postdoc.github.io/BMM/2016/LiveScript/Validation_birdRadar.html)

<img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/10-paper/figure/Figure4.png">
<span style="font-size:0.8em;">Comparison of the estimated bird densities (black line) and their uncertainty range (10-90 quantiles in grey) with the bird densities (red dots) observed using dedicated bird radars at twolocations in (a) Herzeele, France (50◦53’05.6"N 2◦32’40.9"E) and (b) Sempach, Switzerland (47◦07’41.0"N8◦11’32.5"E). Note that because of the power transformation, model uncertainties are larger when the migration intensity is high.  It is therefore critical to account for the uncertainty ranges (light grey) when comparing the interpolation results with the bird radars observations (red dots).</span>

## Estimation and simulation
Estimation on a 3D grid covering Europe is explained in the [MATLAB LiveScript of estimation](https://github.com/Rafnuss-PostDoc/BMM/2016/LiveScript/Estimation_map.m). 

<img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/10-paper/figure/Figure5.png">
<span style="font-size:0.8em;">Maps of bird density estimation every hour of a single night (3-4 October). The sunrise and sunset fronts are visible at 18:00 and 05:00 with lower densities close to the fronts and no value after thefront. The resemblance from hour-to-hour illustrates the high temporal continuity of the model. A rain cell above Poland blocked migration on the Eastern part of the domain. By contrast, a clear pathway is visible from Northern Germany through to Southwestern France.</span>

Simulation on the same grid is performed with SGS, as explained in the [MATLAB code of simulation](https://github.com/Rafnuss-PostDoc/BMM/blob/master/2016/5-Simulation/Simulation_map.m)

<img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/10-paper/figure/Figure6.png">
<span style="font-size:0.8em;">Snapshot of three different realisations showing peak migration (4 October 2016 21:00 UTC). The total number of birds in the air for these realisations was 125, 126 and 122 million respectively. Comparing the similarities and differences of bird density patterns among the realisations illustrates the variability allowed by the stochastic model used. The texture of these realisations is more coherent with the observations than the smooth estimation map in the figure above.</span>


## Results
[A dedicated website was built to interact with the interpolated map of bird migration: ](https://bmm.raphaelnussbaumer.com/)
[<img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/FigureS5-3.png">](https://bmm.raphaelnussbaumer.com/)

| 					| Density [bird/m<sup>2</sup>] | Flight |
| ------------- 	| ------------- 	 | ------------- |
|  Estimation Map 	|  <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/Density_estimationMap_reassamble.gif">  | <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/Flight_estimationMap.gif">  |
| Simulation Map    |  <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/Density_simulationMap_reassemble.gif"> | <img src="https://raw.githubusercontent.com/Rafnuss-PostDoc/BMM/master/2016/figure/Flight_simulationMap.gif"> |

The final interpolated spatio-temporal map can also be downloaded from zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243397.svg)](https://doi.org/10.5281/zenodo.3243397).

<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
<link rel="stylesheet" href="https://cdn.rawgit.com/jpswalsh/academicons/master/css/academicons.min.css">
