# Nanoparticle Tracking Analysis

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/c4664a19-8f4c-4774-a13a-0958ea282b6a

# Contents

- [Overview](#overview)
- [Time lapse data](#time-lapse-data)
- [Diffusing beads example](#diffusing-beads-example)

# Overview

Here we include MATLAB codes and data for the three-dimensional NTA performed on [diffusing beads](#diffusing-beads-example). The analysis was conducted both for sparse and dense samples. Below we provide drive links to measured [time lapse data](#time-lapse-data) with the Tetrapod PSF used for our analysis. This data was first processed with [DeepSTORM3D](https://github.com/EliasNehme/DeepSTORM3D) to recover the 3D localizations at each time step. Afterwards, the localizations were linked to disjoint tracks on a frame-by-frame basis using the [Hungarian algorithm](<NTA/Hungarian algorithm>) as illustrated in the [diffusing beads example](#diffusing-beads-example) below.

# Time lapse data

The folder [`NTA`](NTA) should contain two data subfolders: `Incucyte_TP_low_density` and `Incucyte_TP_high_density`, corresponding to a time lapse of 2D images (saved as individual tif images) of a sparse/dense FOV acquired with the Tetrapod PSF over time. The time span of the experiment was 100 frames with an exposure time of 0.4 seconds, and a delay of 2 seconds in between consecutive frames. These images can be downloaded from google drive using the following url:
```
https://drive.google.com/drive/folders/1uybrj42V_9C1rwZZ8TFTtqeFbRh2ujKO?usp=sharing
```

# Diffusing beads example

As previously mentioned, the first step in the analysis was to apply a trained [DeepSTORM3D](https://github.com/EliasNehme/DeepSTORM3D) to recover the 3D localizations at each time step. The resulting localizations are given in the `localizations_deepstorm3d_*_density.csv` where `*` is either `high` or `low`. The script `NTA_main.m` does the following:
* First, it links the localizations into tracks using the [Hungarian algorithm](<NTA/Hungarian algorithm>). The resulting tracks are thresholded by length, such that only tracks with localizations in at least 50% of the frames are saved as `tracks_TP_*_density.mat`.
* Second, the script computes the ensemble mean square displacement (MSD) in 3D and for each axis separately. The result of this calculation is saved as `ensemble_msds_TP_*_density.mat`.
* Third, the ensemble MSDs are fit using a linear line to determine the diffusion coefficient and the localization precision. The result of this last step is saved as `ensemble_msd_fits_TP_*_density.mat`.

Note that to run this script you do NOT need to download the time lapse data. Everything needed is provided already within this folder.