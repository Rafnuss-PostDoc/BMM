<head>  
  <link rel="shortcut icon" type="image/png" href="https://bmm.raphaelnussbaumer.com/assets/favicon.png">
  </head>



# Bird Migration Map

Bird Migration Modelling is a research project which aims at modeling the spatio-temporal patterns of nocturnal bird migration using weather radars. The main website of this project is [bmm.raphaelnussbaumer.com](http://bmm.raphaelnussbaumer.com/).
All publications and update of the project are available on [ResearchGate](https://www.researchgate.net/project/Bird-Migration-Modelling-BMM).

This repository is divided in the following projects

## Project 1: Interpolation
The first project (`/2016/`)focuses on developing a geostatisctial model to interpolate bird densities at high resolution. Interactive presentation of the methodology is presented at [rafnuss-postdoc.github.io/BMM/2016](https://rafnuss-postdoc.github.io/BMM/2016)

**Presentation:**
> Space-time interpolation of nocturnal bird migration. Raphaël Nussbaumer, Lionel Benoit, Grégoire Mariethoz, Felix Liechti, Baptiste  Schmid. *BOU 2019*. Warwick University. DOI: [10.13140/RG.2.2.11249.53605](https://doi.org/10.13140/RG.2.2.11249.53605).

**Paper:**
> A Geostatistical Approach to Estimate High Resolution Nocturnal Bird Migration Densities from a Weather Radar Network. Nussbaumer, R.; Benoit, L.; Mariethoz, G.; Liechti, F.; Bauer, S.; Schmid, B. *Remote Sens*. **2019**, 11, 2233. DOI:[10.3390/rs11192233](https://doi.org/10.3390/rs11192233).
  
**Demo:**
[birdmigrationmap.vogelwarte.ch/2016](https://birdmigrationmap.vogelwarte.ch/2018/)



## Project 2: Birdflow

In the second project (`/2018/`), we use a flow model (from fluid dynamic) to quantifies the flow of birds taking-off and landing accross western Europe. A short overview of what each script is for is available at [rafnuss-postdoc.github.io/BMM/2018](https://rafnuss-postdoc.github.io/BMM/2018).

**Paper:**
> Quantifying year-round nocturnal bird migration with a fluid dynamics model. Nussbaumer R, Bauer S, Benoit L, Mariethoz G, Liechti F, Schmid B. **2021**. *J. R. Soc. Interface* 18: 20210194. DOI: [10.1098/rsif.2021.0194](https://doi.org/10.1098/rsif.2021.0194)


**Demo:**
[birdmigrationmap.vogelwarte.ch/2018](https://birdmigrationmap.vogelwarte.ch/2018/)


## Project 3: Separation of birds, insects and weathers

This smaller third project is presenting a method to differentiate bird, insect and weather in weather radar data based on doppler product (airpseed and standard deviation of radar velocity. The MATLAB code an be found at [rafnuss-postdoc.github.io/BMM/2018/LiveScript/Insect_removal](https://rafnuss-postdoc.github.io/BMM/2018/LiveScript/Insect_removal.html).

**Paper:**
>  A Gaussian Mixture Model to Separate Birds and Insects in Single-Polarization Weather Radar Data. Raphäel Nussbaumer, Baptiste Schmid, Silke Bauer, Felix Liechti. *Remote Sens*. **2021**, 13(10), 1989. DOI:[10.3390/rs13101989](https://doi.org/10.3390/rs13101989).


  
## Project 4: Windsupport

Study the influence of windspeed on groundspeed and airspeed in time and space. The MATLAB code an be found at [rafnuss-postdoc.github.io/BMM/WindSupport/HTML/script](https://rafnuss-postdoc.github.io/BMM/WindSupport/HTML/script.html).


## Project 5: Particle Tracking

Using the interpolated departure, flight and landing information, we can generate particles in the flow to estimate individual-based information (length of journey, duration,etc...).
This project is still at an early stage and I am not working a lot on it. 
You can find the current code at [rafnuss-postdoc.github.io/BMM/Particle/HTML/Particle_tracking.html](https://rafnuss-postdoc.github.io/BMM/Particle/HTML/Particle_tracking.html) and the webdemo of the particle trajectory at [flowmap.blue](https://flowmap.blue/1de5uGWfZKLLIUqmodfHjps240PC9sRwwp1IqcbVXZRY?v=48.875000,5.375000,4.96,0,0&a=1&as=1&b=1&bo=75&c=1&ca=0&cz=3&d=1&fe=1&lt=0&lfm=ALL&t=20180307T000000,20180315T000000&col=Oranges&f=50)


## Project 6: Windfarm

We combine birdflow map with the windturbine map to estime the number of bird at risk of collision (i.e., birds flying through wind turbine swept area) and simulatanously, the potential power production. The methodology is seperated in three steps:

 1. [Windfarm processing](https://rafnuss-postdoc.github.io/BMM/WindFarm/HTML/1_windfarm_processing)
 2. [Interpolation of bird density at windturbine altitude](https://rafnuss-postdoc.github.io/BMM/WindFarm/HTML/2_interpolate_height_ratio.html)
 3. [Estimate the number of birds at risk](https://rafnuss-postdoc.github.io/BMM/WindFarm/HTML/3_bird_at_risk.html)
