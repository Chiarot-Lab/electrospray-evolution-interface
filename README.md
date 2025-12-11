# Data for "Evolution of Electrosprayed Particles at a Static Air-Water Interface on Multiple Time Scales"

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jcis.2025.139611-blue)](https://doi.org/10.1016/j.jcis.2025.139611)

### Authors:
**Joseph M. Prisaznuk<sup>1</sup>**, **Xin Yong<sup>2</sup>**, **Paul R. Chiarot<sup>1\*</sup>**

<sup>1</sup>Department of Mechanical Engineering, Binghamton University, Binghamton, NY 13902, USA  
<sup>2</sup>Department of Mechanical and Aerospace Engineering, University at Buffalo, Buffalo, NY 14260, USA  

## Directory Structure

| Folder             | Description                                                                                     |
|--------------------|-------------------------------------------------------------------------------------------------|
| `data/`            | Processed data saved in `.mat` format, partial raw data in `.tif` image format                  |
| `image_overlays/`  | MATLAB scripts to load `.tif` images and overlay the elapsed time, using saved metadata         |
| `msd/`             | Mean Squared Displacement (MSD) analysis                                                        |
| `rdf_psi6/`        | Radial Distribution Function (RDF) and sixfold orientation order parameter $\psi_6$ analysis    |
| `UV-vis_spectra/`  | Results of ultraviolet to visible (UV-vis) spectrometry                                         |

## Visualizing the Data

1. Navigate to the relevant folder.
2. Open the code file `*.m`.
3. Run the script in MATLAB.

## Notes
- Shared data used by multiple scripts is in the top-level `data/` directory. Relative paths are used in the code to handle this as needed. 
- Full datasets are available upon reasonable request; please contact the corresponding author.

## Citation

If you use this data or code, please cite:

Prisaznuk, J. M.; Yong, X.; Chiarot, P. R. Evolution of Electrosprayed Particles at a Static Air-Water Interface on Multiple Time Scales. *Journal of Colloid and Interface Science* **2026**, *706*, 139611. https://doi.org/10.1016/j.jcis.2025.139611

## License

<a href="https://github.com/Chiarot-Lab/electrospray-evolution-interface">electrospray-evolution-interface</a> © 2025 is licensed under <a href="https://creativecommons.org/licenses/by-nc/4.0/">CC BY-NC 4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">
