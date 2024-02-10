# EDOF imaging and deconvolution


# Contents

- [Overview](#overview)
- [Beads in a gel](#beads-in-a-gel)
- [Spheroid measurement](#spheroid-measurement)
- [Deconvolution example](#deconvolution-example)

# Overview

Here we include MATLAB codes and data for the EDOF quantitative analysis performed on [beads suspeneded in a gel](#beads-in-a-gel). In addition, we provide Drive links to raw measurements of a [spheroid](#spheroid-measurement) with the standard PSF (unmodified objective lens) and an EDOF PSF. Finally, we provide a [demo example](#deconvolution-example) of applying a few iterations of Lucy-Richardson Deconvolution to improve the resolution of the EDOF measurements.

# Beads in a gel

The folder `Beads/Low Density/` should contain two data subfolders: `EDOF` and `ST` corresponding to z-stacks (saved as individual tif images) of the same FOV one acquired with the standard PSF (ST) and the other with an EDOF PSF. These zstacks can be downloaded from google drive using the following url:
```
https://drive.google.com/drive/folders/1tSzk91CfnRzaedVWYDP_33X9TeWfexG5?usp=sharing
```
[Thunderstorm](https://zitmen.github.io/thunderstorm/) was run on each zstack to quickly find the beads positions. These are saved in excel files `ThunStormMAX*.csv`. In addition Thunderstorm proximity filter was applied to avoid beads that are too close. Any group of proximal beads is registered in excel files `*close_locs.csv`, specifying beads that will be excluded from analysis. To reproduce the analysis from the paper, you need to run the scripts `Analyze_multiple_beads.m` and `Display_DOF_Multiple_beads.m`:
* `Analyze_multiple_beads.m` - This code crops the useful beads, and performs 2D gaussian fit per frame around the focus of each bead (takes a long time to run). The generated results are `Crops_*.mat` (fast) and `fitted*.mat` (slow), where `*` is either `ST` or `EDOF`.
* `Display_DOF_Multiple_beads.m` - This code displays the results from `fitted*.mat` and performs 1D gaussian fitting to the recovered PSF widths (also quite slow). The results are saved as `multi_bead_fits_results.mat`.

Expected results of running these scripts are included inside `Beads/Low Density/`.

# Spheroid measurement

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/9b4106a8-d516-45ed-8afa-074db0bce2b7

The raw data of a spheroid measured with both the standard PSF (unmodified objective lens) and the EDOF PSF can be downloaded from google drive using the following url:
```
https://drive.google.com/drive/folders/1eUIc_yBiKpVcyHe7Bf3BLBDfWKuwF-my?usp=sharing
```

# Deconvolution example

`NoamXX.m` applies Lucy-Richardson deconvolution to the measured spheroid with the EDOF PSF to improve its lateral resolution. Note that to run this demo you are required to download the EDOF data from the [url](https://drive.google.com/drive/folders/1eUIc_yBiKpVcyHe7Bf3BLBDfWKuwF-my?usp=sharing) in the previous section.
