# Compact PSF Engineering

This code accompanies the paper: "Depth-enhanced high throughput microscopy by compact PSF engineering"

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/a6422ffc-1f4f-41f1-8d2a-89b1da81cd53

# Contents

- [Overview](#overview)
- [EDOF imaging and deconvolution](#edof-imaging-and-deconvolution)
- [Snapshot 3D imaging and CellSnap](#snapshot-3d-imaging-and-cellsnap)
- [System requirements](#system-requirements-and-installation-instructions)
- [Code structure](#code-structure)
- [Experimental dataset](#experimental-dataset)
- [Demo example](#demo-example)
- [Nanoparticle tracking analysis](#nanoparticle-tracking-analysis)

# Overview

Our proposed hardware prototype for compact PSF engineering is compatible with any phasemask/PSF. Specifically, in the paper we demonstrated two common applications: (i) Extended-Depth-Of-Field (EDOF) imaging using a [depth-insensitive PSF](https://ieeexplore.ieee.org/document/9439955), and (ii) Snapshot 3D imaging using the [Tetrapod](https://pubs.acs.org/doi/10.1021/acs.nanolett.5b01396) PSF. In [EDOF imaging](#edof-imaging-and-deconvolution), the resulting measurements can be used either directly or after a few iterations of Lucy-Richardson deconvolution. As for [snapshot 3d imaging](#snapshot-3d-imaging-and-cellsnap), the resulting 2D measurements need to be post-processed to extract the underlying 3D information. For spheroid imaging, this is achieved with [CellSnap](#snapshot-3d-imaging-and-cellsnap); A tailored deep neural network architecture trained on an [experimental dataset](#experimental-dataset) of matched 2D inputs and 3D outputs. Alternatively, for [nanoparticle tracking analysis](#nanoparticle-tracking-analysis), the measurements are processed with [DeepSTORM3D](https://github.com/EliasNehme/DeepSTORM3D), and the resulting tracks are linked sequentially using the [Hungarian algorithm](<NTA/Hungarian algorithm>).

# EDOF imaging and deconvolution

The folder [`EDOF`](EDOF) includes MATLAB codes and links to experimental measurements with two different realizations of EDOF PSFs. The first set of measurements is of beads embedded in a 3D gel. This sample was used to quantitatively assess the gain in performance brought by the EDOF PSF compared to a standard unmodified objective lens. The second set of measurements includes the application of an EDOF PSF to spheroid imaging, with the possibility of further improving the result via Lucy-Richardson deconvolution.

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/9b4106a8-d516-45ed-8afa-074db0bce2b7

# Snapshot 3D imaging and CellSnap

The training of CellSnap is comprised of two phases: training a focus finder and afterward training a conditional 3D segmentation model. The dataset for training/testing has been curated from multiple scans of four 96-well plates. After discarding non-spherical and low snr spheroids, the resulting dataset consisted of 592 spheroid "views", out of which we used 532 for training and 60 for validation. `Conditional_3D_segmentation_testing.ipynb` demonstrates the application of a pre-trained CellSnap model on another 20 test spheroids not seen during training/validation.

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/91268082-5d8c-40db-b526-355306f0c5e2

# System requirements and installation instructions
* CellSnap was tested on a *Linux* system with Ubuntu version 18.0, equipped with an Nvidia Titan RTX GPU with 24 GB of memory.
* The [conda](https://docs.conda.io/en/latest/) environment for this project is given in `environment.yml`. To replicate the environment on a Linux system, use the command: `conda env create -f environment.yml` from within this directory. This should take a couple of minutes.
* After activating the environment using `conda activate cellsnap`, you're set to go.

# Code structure

* Training
    * `Focus_finder_training.ipynb` implements the training of the focus finder using pre-processed z-stacks of the [Tetrapod](https://pubs.acs.org/doi/10.1021/acs.nanolett.5b01396) PSF.
    * `Conditional_3D_segmentor_training.ipynb` implements the training of the conditional 3D segmentation model given the pre-trained focus finder and matched pairs of Tetrapod z-stacks and the corresponding 3D segmentation outputs obtained by post-processing the results of [Cellpose](https://www.nature.com/articles/s41592-020-01018-x) applied to the accompanying standard PSF z-stacks.
* Testing
    * `Focus_finder_testing.ipynb` implements the testing of the focus finder.
    * `Conditional_3D_segmentor_testing.ipynb` implements the testing of CellSnap after training both components.
 
 # Experimental dataset
* The [`Exp Dataset`](<Exp Dataset>) folder includes the following:
    * `Training` folder - contains 592 z-stacks of spheroid views with the Tetrapod PSF under the `TP` folder, the standard PSF under the `ST` folder, and the corresponding segmentations under the `CM` folder.
    * `Testing` folder - contains 20 z-stacks of spheroid views with the Tetrapod PSF under the `TP` folder, the standard PSF under the `ST` folder, and the corresponding segmentations under the `CM` folder.

Note that the [`Exp Dataset`](<Exp Dataset>) folder should be in the same working directory of the notebooks for running the code.

 # Demo example

`Conditional_3D_segmentor_testing.ipynb` demonstrates the applicability of CellSnap on 20 test spheroids with demo pre-trained weights hardcoded in the notebook. Running the notebook on a linux system with Titan RTX GPU takes ~2 mins. The expected results are embedded in the saved notebook.

* The `Models` folder includes the following:
    * `focus_finder_model_best_ckpt.pth` - pre-trained model weights of the focus finder.
    * `conditional_3D_segmentor_model_best_ckpt.pth` - pre-trained model weights of the conditional 3D segmentation model.

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/cbfbdc75-b444-491b-adc4-abcf90a61798

# Nanoparticle tracking analysis

In the paper, we also tested the Tetrapod PSF for the application of three-dimensional nanoparticle tracking analysis (NTA). To achieve this, we measured a sample of diffusing beads with the Tetrapod PSF and applied [DeepSTORM3D](https://github.com/EliasNehme/DeepSTORM3D) to the resulting time-lapse of 2D images. This resulted in 3D localizations over time. These were then linked to 3D tracks using a custom MATLAB code. The folder [`NTA`](NTA/) includes the relevant code and links to a measured experimental time-lapse with the Tetrapod PSF. See inside for more details.

https://github.com/EliasNehme/HTPmicroscopy/assets/32178070/c4664a19-8f4c-4774-a13a-0958ea282b6a
