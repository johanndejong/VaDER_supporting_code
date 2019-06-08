Supporting code for the manuscript "Deep learning for clustering of multivariate clinical patient trajectories with missing values".

The code depends on data from ADNI and PPMI (http://adni.loni.usc.edu/ and https://www.ppmi-info.org/ ), which the license agreement does not allow me to make public here. Hence, I have supplied artificial patient data as input for the following scripts:

- ADNI_hyperparameter_optimization.r
- ADNI_optimal_model.r
- PPMI_hyperparameter_optimization.r
- PPMI_optimal_model.r

The artifical data has been randomly sampled from the latent Gaussian mixture distribution that we learn as part of training VaDER (https://github.com/johanndejong/VaDER) on the ADNI and PPMI data, and therefore represents the original data very well, also in terms of missing values. VaDER is needed to run these scripts.

