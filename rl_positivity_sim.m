%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple demonstration of the effect of positivity in RL deconvolution %
% James Manton, 2019 - Founder and License Holder
% Brian Northan 2019 - Contributors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Simulation parameters
max_photons = 100;
num_iter = 100;
pixel_size = 20;
spacing_px = 4;
n = 512;
lambda = 510;
numerical_aperture = 1.4;
background_level = 0;
left_bg = 0;
mid_bg = 0;%0.05;
right_bg = 0;%0.25;

SAVE_DISPLAY = 1;
USE_GPU = 0;


%% Create OTF
otf = paraxial_otf(n, lambda, numerical_aperture, pixel_size);


%% Create line pairs and add background levels
field = zeros(n);
field((n/4):(3*n/4), (n/2) - spacing_px) = 1;
field((n/4):(3*n/4), (n/2) + spacing_px) = 1;
field = field + circshift(field, [0, round(n/3)]) + circshift(field, [0, -round(n/3)]);
 
field(7*n/8 - spacing_px, (n/4):(3*n/4)) = 1;
field(7*n/8 + spacing_px, (n/4):(3*n/4)) = 1;
 
field(n/2 - spacing_px, (n/16):(15*n/16)) = 1;
field(n/2 + spacing_px, (n/16):(15*n/16)) = 1;
 
field(:, round(n/3):round(2*n/3)) = field(:, round(n/3):round(2*n/3)) + mid_bg;
field(:, round(2*n/3):end) = field(:, round(2*n/3):end) + right_bg;


%% Simulate captured data
field_imaged = real(ifft2(fft2(field) .* otf));
field_imaged = field_imaged ./ max(field_imaged(:));
field_imaged = poissrnd(field_imaged * max_photons + background_level);


%% Deconvolve data
if USE_GPU
    field_rl = gather(richardson_lucy(gpuArray(field_imaged), gpuArray(otf), num_iter, 1));
else
    field_rl = richardson_lucy(field_imaged, otf, num_iter, 1);
end
field_rl = field_rl ./ max(field_rl(:));

ssimval_rl=ssim(field, field_rl);

%% Calculate spectra
spectrum_field = log(1 + abs(fftshift(fft2(field))));
spectrum_field = spectrum_field ./ max(spectrum_field(:));
spectrum_rl = log(1 + abs(fftshift(fft2(field_rl))));
spectrum_rl = spectrum_rl ./ max(spectrum_rl(:));


%% Display (and save) results
figure(1)
display_array = [field ./ max(field(:)), spectrum_field, fftshift(otf); ...
    field_imaged ./ max(field_imaged(:)), field_rl, spectrum_rl];
imshow(display_array, [])

if SAVE_DISPLAY
   imwrite(display_array, 'rl_positivity_sim.png')
end

for i=1:6
    rf = 1/(10^(i-1))
    field_wnr{i} = deconvwnr(field_imaged, fftshift(ifftn(otf)), rf);
    field_wnr{i}=field_wnr{i}./max(field_wnr{i}(:));
    spectrum_wnr{i} = log(1 + abs(fftshift(fft2(field_wnr{i}))));
    spectrum_wnr{i} = spectrum_wnr{i} ./ max(spectrum_wnr{i}(:));
    ssimval_wnr(i) = ssim(field,field_wnr{i})
end

figure(2)
display_array_wnr = [ field_wnr{1}, field_wnr{2}, field_wnr{3}; field_wnr{4}, field_wnr{5}, field_wnr{6}];
imshow(display_array_wnr);
title('deconv wnr NSR=.1 to .000001');

figure(3)
display_array_spectrum_wnr = [ spectrum_wnr{1}, spectrum_wnr{2}, spectrum_wnr{3}; spectrum_wnr{4}, spectrum_wnr{5}, spectrum_wnr{6}];
imshow(display_array_spectrum_wnr);
title('specturm wnr NSR=.1 to .000001');

ssimval_rl
ssimval_wnr
