# R code for the paper Cretois et al. 'A suggested protocol for using Citizen Science data in habitat selection studies'

---

This repository hosts functions . 

**Note that to run the "real analysis"** the raster used to extract the covariates, the telemetry and Citizen Science dataset were too heavy to be shared on GitHub. You can find all this data at XXX. **Place the downloaded folders in data to smoothly run the scripts**

This work was supported by a PhD grant funded by the Norwegian University of Science and Technology and the Research Council of Norway.

---

## Simulation study






## Analysis with real data


To run the whole analysis without going through the individual functions you can use the [run_all](./real_study/run_all.R) functions.

### Data wrangling

Functions below aim to create a simplified and standardized dataset between species.

[make_moose_dataset](./real_study/custom_functions/make_moose_dataset.R): Open the moose GPS telemetry data, the citizen science dataset and create a simplified dataset that contains only summer observations. Buffer around moose minimum convex polygon can be specified to sample more or less citizen science observations. Thinning rate can also be specified to include more or less GPS telemetry locations.

[make_roedeer_dataset](./real_study/custom_functions/make_roedeer_dataset.R): Open the roe deer GPS telemetry data and VHF data, the citizen science dataset and create a simplified dataset that contains only summer observations. Buffer around moose minimum convex polygon can be specified to sample more or less citizen science observations. Thinning rate can also be specified to include more or less GPS telemetry locations.

[make_reindeer_dataset](./real_study/custom_functions/make_moose_dataset.R): Open the wild reindeer GPS telemetry data, the citizen science dataset and create a simplified dataset that contains only summer observations. Buffer around moose minimum convex polygon can be specified to sample more or less citizen science observations. Thinning rate can also be specified to include more or less GPS telemetry locations.

### 
































