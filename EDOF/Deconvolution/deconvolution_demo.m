%% Lucy-Richardson deconvolution

% start with a clean slate
close all; clear all; clc;

%% Standard deconvolution with Gaussian PSF

% reading image
ST_image_path = 'input/Standard_frame.tif';
ST_image = double(imread(ST_image_path));

% optional - cropping image
ST_image_crop_idx_x = 475;
ST_image_crop_idx_y = 685;
crop_size = 400;
ST_image = ST_image(ST_image_crop_idx_x:ST_image_crop_idx_x+crop_size, ST_image_crop_idx_y:ST_image_crop_idx_y+crop_size);

% Gaussian parameters for Gaussian PSF model
gaussian_parameters.amplitude = 800;
gaussian_parameters.stddev = 0.6;
gaussian_parameters.mean = 12;
gaussian_parameters.size = 23;

% other parameters
bg_removal_flag = 1;
save_image_flag = 1;
output_path = 'output/ST_deconvolved_image.tif';

% deconvolve and save output
ST_deconvolved_image = deconvolve_ST(ST_image, gaussian_parameters, bg_removal_flag, save_image_flag, output_path);

%% EDOF deconvolution with Lorentzian PSF

% reading image
EDOF_image_path = 'input/EDOF_frame.tif';
EDOF_image = double(imread(EDOF_image_path));

% optional - cropping image
EDOF_image_crop_idx_x = 384;
EDOF_image_crop_idx_y = 523;
crop_size = 400; % pixels
EDOF_image = EDOF_image(EDOF_image_crop_idx_x:EDOF_image_crop_idx_x+crop_size, EDOF_image_crop_idx_y:EDOF_image_crop_idx_y+crop_size);

% Lorentzian parameters for Lorentzian PSF model
lorentzian_parameters.amplitude = 515;
lorentzian_parameters.gamma = 0.8;
lorentzian_parameters.mean = 12;
lorentzian_parameters.size = 23;

% other parameters
bg_removal_flag = 1;
save_image_flag = 1;
output_path = 'output/EDOF_deconvolved_image.tif';

% deconvolve and save output
EDOF_deconvolved_image = deconvolve_EDOF(EDOF_image, lorentzian_parameters, bg_removal_flag, save_image_flag, output_path);
