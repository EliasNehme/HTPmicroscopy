function ST_deconvolved_image = deconvolve_ST(ST_image, gaussian_parameters, bg_removal_flag, save_image_flag, output_path)

% performs Lucy-Richardson deconvolution to ST image using a Gaussian
% PSF model.

% input:
% ST_image - ST image, double precision
% gaussian_parameters - struct array of the Gaussian parameters: 
%   amplitude - default: 800
%   stddev - default: 0.6
%   mean - mean of x and y. default: 12
%   size - size of the cropped PSF image. default: 23
% bg_removal_flag - whether to remove background from ST image.
% save_image_flag - whether to save ST deconvolved image.
% output_path - output path in which to save the ST deconvolved image. 

% output:
% ST_deconvolved_image - deconvolved ST image

% Gaussian parameters
amplitude = gaussian_parameters.amplitude;
stddev = gaussian_parameters.stddev;
mean = gaussian_parameters.mean;
size = gaussian_parameters.size;

% Gaussian PSF model
[x, y] = meshgrid(1:size, 1:size);
gaussian_PSF = amplitude * exp(-((x - mean).^2 / (2 * stddev^2) + (y - mean).^2 / (2 * stddev^2)));
gaussian_PSF = round(gaussian_PSF);

% ST image processing
if bg_removal_flag
    % bg removal
    bg = ST_image(1,1) + ST_image(1,end) + ST_image(end,1) + ST_image(end,end) + ST_image(2,1) + ST_image(2,end) + ST_image(end-1,1) + ST_image(end-1,end) + ST_image(1,2) + ST_image(1,end-1) + ST_image(end,2) + ST_image(end,end-1) + ST_image(2,2) + ST_image(2,end-1) + ST_image(end-1,2) + ST_image(end-1,end-1);
    bg = round(bg/16);
    ST_image = max(0, ST_image-bg);
end

% deconvolution
ST_deconvolved_image = deconvlucy(ST_image./norm(ST_image(:)),gaussian_PSF./norm(gaussian_PSF(:)), 5);

% saving image
if save_image_flag
    % image scaling
    min_value = min(ST_deconvolved_image(:));
    max_value = max(ST_deconvolved_image(:));
    scaled_ST_deconvolved_image = uint16((ST_deconvolved_image - min_value) / (max_value - min_value) * 4096);
    
    % saving TIF image
    imwrite(scaled_ST_deconvolved_image, output_path, 'tif');
end

end