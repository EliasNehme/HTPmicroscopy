function EDOF_deconvolved_image = deconvolve_EDOF(EDOF_image, lorentzian_parameters, bg_removal_flag, save_image_flag, output_path)

% performs Lucy-Richardson deconvolution to EDOF image using a Lorentzian
% PSF model.

% input:
% EDOF_image - EDOF image, double precision
% lorentzian_parameters - struct array of the Lorentzian parameters: 
%   amplitude - default: 515
%   gamma - default: 0.8
%   mean - mean of x and y. default: 12
%   size - size of the cropped PSF image. default: 23
% bg_removal_flag - whether to remove background from EDOF image.
% save_image_flag - whether to save EDOF deconvolved image.
% output_path - output path in which to save the EDOF deconvolved image. 

% output:
% EDOF_deconvolved_image - deconvolved EDOF image

% Lorentzian parameters
amplitude = lorentzian_parameters.amplitude;
gamma = lorentzian_parameters.gamma;
mean = lorentzian_parameters.mean;
size = lorentzian_parameters.size;

% Lorentzian PSF model
[x, y] = meshgrid(1:size, 1:size);
lorentzian_PSF= amplitude ./ ((x - mean).^2 + (y - mean).^2 + gamma^2);
lorentzian_PSF = round(lorentzian_PSF);

% EDOF image processing
if bg_removal_flag
    % bg removal
    bg = EDOF_image(1,1) + EDOF_image(1,end) + EDOF_image(end,1) + EDOF_image(end,end) + EDOF_image(2,1) + EDOF_image(2,end) + EDOF_image(end-1,1) + EDOF_image(end-1,end) + EDOF_image(1,2) + EDOF_image(1,end-1) + EDOF_image(end,2) + EDOF_image(end,end-1) + EDOF_image(2,2) + EDOF_image(2,end-1) + EDOF_image(end-1,2) + EDOF_image(end-1,end-1);
    bg = round(bg/16);
    EDOF_image = max(0, EDOF_image-bg);
end

% deconvolution
EDOF_deconvolved_image = deconvlucy(EDOF_image./norm(EDOF_image(:)),lorentzian_PSF./norm(lorentzian_PSF(:)), 5);

% saving image
if save_image_flag
    % image scaling
    min_value = min(EDOF_deconvolved_image(:));
    max_value = max(EDOF_deconvolved_image(:));
    scaled_EDOF_deconvolved_image = uint16((EDOF_deconvolved_image - min_value) / (max_value - min_value) * 4096);
    
    % saving TIF image
    imwrite(scaled_EDOF_deconvolved_image, output_path, 'tif');
end

end