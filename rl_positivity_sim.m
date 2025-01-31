%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple demonstration of the effect of positivity in RL deconvolution %
% James Manton, 2019 - Founder and License Holder
% Brian Northan 2019 - Contributors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%% Simulation parameters
max_photons = 1000;

% percentage of photons detected 0-1
photons_detected = 1;
pixel_size = 20;
spacing_px = 4;
n = 512;
lambda = 510;
numerical_aperture = 1.4;
background_level = 0;

method=3;

if method==1
    num_iter=10000
    experiment_description = 'classic_rl';
elseif method ==2
    % matlab deconvlucy is accelerated so use less iterations
    num_iter=100
    experiment_description = 'deconvlucy';
elseif method ==3
    % matlab deconvblind is accelerated so use less iterations
    num_iter=100
    experiment_description = 'deconvblind';
end

% 0 for lines, 1 for points, 2 for circles of varying size, 3 for circles
% of varying intensity
sim_type=1;

if sim_type==0
    experiment_description = [experiment_description '_lines'];
elseif sim_type==1
    experiment_description = [experiment_description '_points'];
elseif sim_type==2
    experiment_description = [experiment_description '_circles_varyingsize'];
elseif sim_type==3
    experiment_description = [experiment_description '_circles_varyingintensity'];
end

experiment_description = [experiment_description '_' num2str(num_iter)];
experiment_dir = ['reports/' experiment_description '/'];

mkdir([experiment_dir]);
fileID = fopen([experiment_dir 'report.md'],'w');
fprintf(fileID,'## Parameters  \n');
fprintf(fileID,'description: %s  \n',experiment_description);
fprintf(fileID,'max photons %d  \n',max_photons);
fprintf(fileID,'num iter %d  \n  ',num_iter);
fprintf(fileID,'pixel size %d  \n  ',pixel_size);
fprintf(fileID,'spacing px %d nm  \n  ',spacing_px);
fprintf(fileID,'n %d  \n',n);
fprintf(fileID,'lambda %d  \n',lambda);
fprintf(fileID,'numerical aperture %f  \n',numerical_aperture);
fprintf(fileID,'background level %d  \n',background_level);

%fprintf(fileID,'%6.2f %12.8f\n',A);

% even if bg is 0 noise will be added to signal
% so set to 0 for completely noiseless simulation
add_noise = true;

% add option to turn off normalization of convolved image
normalize_convolved = false;

background = true;
if background
    left_bg = 0;
    mid_bg = 0.05;
    right_bg = 0.25;
else
    left_bg = 0;
    mid_bg = 0;
    right_bg = 0;
end

fprintf(fileID,'left background %f  \n',left_bg);
fprintf(fileID,'mid background %f  \n',mid_bg);
fprintf(fileID,'right background %f  \n',right_bg);

SAVE_DISPLAY = 1;
USE_GPU = 0;

%% Create OTF
otf = paraxial_otf(n, lambda, numerical_aperture, pixel_size);

field = zeros(n);
    
experiment_info='';
if (sim_type==0)
    %% Create line pairs and add background levels
    field((n/4):(3*n/4), (n/2) - spacing_px) = 1;
    field((n/4):(3*n/4), (n/2) + spacing_px) = 1;
    field = field + circshift(field, [0, round(n/3)]) + circshift(field, [0, -round(n/3)]);

    field(7*n/8 - spacing_px, (n/4):(3*n/4)) = 1;
    field(7*n/8 + spacing_px, (n/4):(3*n/4)) = 1;

    field(n/2 - spacing_px, (n/16):(15*n/16)) = 1;
    field(n/2 + spacing_px, (n/16):(15*n/16)) = 1;

elseif (sim_type==1)
    numPairs=10;
    A=1;
    
    experiment_info=[experiment_info 'points '];

    for i=1:numPairs
        field((i*floor(n/numPairs))-floor(n/numPairs/2), (n/2) - i) = A;
        field((i*floor(n/numPairs))-floor(n/numPairs/2), (n/2) + i) = A;
        distance = (2*i+1)*pixel_size;
        if i==numPairs
            experiment_info=[experiment_info num2str(distance) ' nm apart.  \n'];
        else
            experiment_info=[experiment_info num2str(distance) ', '];
        end
        
    end

    field = field + circshift(field, [0, round(n/3)]) + circshift(field, [0, -round(n/3)]);
elseif (sim_type==2)
    x = -n/2:n/2-1;
    y = -n/2:n/2-1;
    [xx, yy] = meshgrid(x+0.5,y+0.5);
    
    numCircles=10;
    A=1;
    for i=1:numCircles
        yc=i*floor(n/numCircles)-floor(n/numCircles/2)-n/2;
        field((xx.^2+(yy-yc).^2)<(i)^2)=A;   % 
    end
    field = field + circshift(field, [0, round(n/3)]) + circshift(field, [0, -round(n/3)]);
elseif (sim_type==3)
    x = -n/2:n/2-1;
    y = -n/2:n/2-1;
    [xx, yy] = meshgrid(x+0.5,y+0.5);
    
    r=8;
    numCircles=10;
   
    for i=1:numCircles
        yc=i*floor(n/numCircles)-floor(n/numCircles/2)-n/2;
        field((xx.^2+(yy-yc).^2)<(r)^2)=2^i/max_photons;   % radius 100, center at the origin
    end
    field = field + circshift(field, [0, round(n/3)]) + circshift(field, [0, -round(n/3)]);
end

field(:, round(n/3):round(2*n/3)) = field(:, round(n/3):round(2*n/3)) + mid_bg;
field(:, round(2*n/3):end) = field(:, round(2*n/3):end) + right_bg;

%% Simulate captured data
field=field*max_photons;
field_imaged = real(ifft2(fft2(field) .* otf));

if normalize_convolved == true
    field_imaged = field_imaged ./ max(field_imaged(:));
end

% add noise
if (add_noise==true)
    field_imaged = poissrnd(field_imaged * photons_detected + background_level);
end

if (method==1)
    %% Deconvolve data
    if USE_GPU
        field_rl = gather(richardson_lucy(gpuArray(field_imaged), gpuArray(otf), num_iter, 1));
    else
        field_rl = richardson_lucy(field_imaged, otf, num_iter, 1);
    end
elseif method==2
    field_rl = deconvlucy(field_imaged, fftshift(ifftn(otf)), num_iter);
elseif method==3
    psf_start = fftshift(ifftn(otf));
    [field_rl, psf_blind] = deconvblind(field_imaged, psf_start, num_iter);
    psf_array = [psf_start(n/4:3*n/4,n/4:3*n/4)./max(psf_start(:)) psf_blind(n/4:3*n/4,n/4:3*n/4)./max(psf_blind(:))];
end

line_plot_fig=figure
hold on
for i=1:9
        yc=(i*floor(n/10))-floor(n/10/2)
        n/2
        xshift=n/3;
        subplot(3,3,i);
        hold on
        plot(field(yc, n/2-10-xshift:(n/2)+10-xshift), 'r');
        plot(field_imaged(yc, n/2-10-xshift:(n/2)+10-xshift), 'g');
        plot(field_rl(yc, n/2-10-xshift:(n/2)+10-xshift), 'b');
end

%legend('field', 'imaged', 'deconvolved');
line_plot_fix.WindowState='maximuzed';
saveas(line_plot_fig, [experiment_dir 'line_plots.png']);
fprintf(fileID,'## Line Plots (medium noise)  \n');
fprintf(fileID,experiment_info);
fprintf(fileID,'red: ground truth  \n');
fprintf(fileID,'green: imaged   \n');
fprintf(fileID,'blue: restored   \n');

fprintf(fileID,'![](line_plots.png)  \n');

field_rl = field_rl ./ max(field_rl(:));
ssimval_rl=ssim(field, field_rl);

%% Calculate spectra
spectrum_field = log(1 + abs(fftshift(fft2(field))));
spectrum_field = spectrum_field ./ max(spectrum_field(:));
spectrum_rl = log(1 + abs(fftshift(fft2(field_rl))));
spectrum_rl = spectrum_rl ./ max(spectrum_rl(:));

%% Display (and save) results
figure
display_array = [field ./ max(field(:)), spectrum_field, fftshift(otf); ...
    field_imaged ./ max(field_imaged(:)), field_rl, spectrum_rl];
imshow(display_array, [])

if SAVE_DISPLAY
   imwrite(display_array, 'rl_positivity_sim.png')
end

imwrite(display_array, [experiment_dir 'images.png'])
fprintf(fileID,'## Images  \n');
fprintf(fileID,'![](images.png)  \n');

if method==3
    imwrite(psf_start(n/8:7*n/8,n/8:7*n/8)./max(psf_start(:)), [experiment_dir 'psfstart.png'])
    imwrite(psf_blind(n/8:7*n/8,n/8:7*n/8)./max(psf_blind(:)), [experiment_dir 'psfblind.png'])
    fprintf(fileID,'## PSFs  \n');
    fprintf(fileID,'left starting PSF, right blind PSF  \n');
    fprintf(fileID,'![](psfstart.png) ![](psfblind.png)  \n');
end


fclose(fileID);


for i=1:6
    rf = 1/(10^(i-1))
    field_wnr{i} = deconvwnr(field_imaged, fftshift(ifftn(otf)), rf);
    field_wnr{i}=field_wnr{i}./max(field_wnr{i}(:));
    spectrum_wnr{i} = log(1 + abs(fftshift(fft2(field_wnr{i}))));
    spectrum_wnr{i} = spectrum_wnr{i} ./ max(spectrum_wnr{i}(:));
    ssimval_wnr(i) = ssim(field,field_wnr{i})
end

figure
display_array_wnr = [ field_wnr{1}, field_wnr{2}, field_wnr{3}; field_wnr{4}, field_wnr{5}, field_wnr{6}];
imshow(display_array_wnr);
title('deconv wnr NSR=1 to .00001');

figure
display_array_spectrum_wnr = [ spectrum_wnr{1}, spectrum_wnr{2}, spectrum_wnr{3}; spectrum_wnr{4}, spectrum_wnr{5}, spectrum_wnr{6}];
imshow(display_array_spectrum_wnr);
title('specturm wnr NSR=1 to .00001');

ssimval_rl
ssimval_wnr
