# digidolph
Spatially explicit null association model for dolphin data

<a href="https://zenodo.org/badge/latestdoi/95614472"><img src="https://zenodo.org/badge/95614472.svg" alt="DOI"></a>

This code can be used to re-create the analysis in Strickland et al. 2017. The model is inspired by and built upon the R package Digiroo2 with expansions for data collected under various sampling protocols. The product is simulated association matrices for indiviuals based on a set of re-locations for each individual.

The code utilizes the following R packages. See sessionInfo for most recent versions and platform tested. 

adehabitatHR
maptools
rgdal
spatstat
Digiroo2
coda
spdep
raster
PBSmapping
gdata
rgeos

References

Strickland, K., Levengood, A. Foroughirad, V. Mann, J. Krzyszczyk, E. & Frere, C. 2017. A framework for the identification of long-term avoidance in longitudinal datasets. Royal Society Open Science. 4(8):170641.

Dwyer R, Best E, Goldizen A. Digiroo2: An application programming interface for generating null models of social contact based on individuals' space use. R package version 520 0.5. 2013.
