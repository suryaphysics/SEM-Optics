close all;
clc;
clear all;
%% CODE AUTHOR - Surya Kamal , surya.kamal@rit.edu

%% Phase Retrieval of the electron probe wavefunction 
% INPUT --> Focal series of the probe intensities  
% GOAL --> Visualization of the Point-Spread Function of SEM Optics

%% Section 1 - Read images (probe intensities at different z-planes) 

N = 256; % Define the image size (image dimension is already 256x256 pixels)

% Define the folder containing the images
imageFolder = 'probe_intensity_focal_series';

% total number of plane is always odd to get a central plane, and is
% either 17,19,21 in our experiments. 
n = 9; % Specify the range of image indices (e.g., from -n to +n)
num_planes = 2 * n + 1; % For example, this will give images from -n to +n (total 2n+1 images)

% Initialize the focal_series matrix (size: N x N x (2n+1))
intensity_focal_series = zeros(N, N, num_planes);

% Loop through the probe intensity indices and read the images
for i = -n:n
    % Correct the index for MATLAB indexing (shift by +n to map -n to 1, ..., +n to (2n+1))
    idx = i + n + 1;  % This maps -n to 1, ..., +n to (2n+1)

    % Handle the special case for i = 0
    if i == 0
        filename = '0_PSF.png';
    else
        filename = sprintf('%+d_PSF.png', i);
    end

    % Construct the full file path
    filePath = fullfile(imageFolder, filename);

    % Check if the file exists before attempting to read it
    if exist(filePath, 'file')
        % Read the image, and store it in the appropriate slice of the matrix
        intensity_focal_series(:,:,idx) = im2double(imread(filePath));
    else
        warning('File %s does not exist in the folder "%s". Skipping...', filename, imageFolder);
    end

end

% ~~~~~~~~~

% Visualization of probe intensities in z-plane
figure(1);
for i=1:num_planes
    imshow(intensity_focal_series(:,:,i));
    title("Plane "+i)
end

%% Section 2 - Optical and Sampling Parameters (x,y,z)

load('z_sampling.mat');                 % Load z-sampling data (longitudnal)
V = 20e+03;                             % Gun voltage
lambda = volt2wavelen(V);               % Relativistic wavelength
k = 2*pi/lambda;                        % Wavenumber

L = 256e-9;                             % Observation plane length
dx = L/N;                               % image plane sample interval
x = -L/2 : dx : L/2 - dx;               % array of sampling locations 
y = x;
[X,Y] = meshgrid(x,y);                  % Meshgrid of x-y location (transverse)
fx = -1/(2*dx) : 1/L : 1/(2*dx)-(1/L);  % Frequency coordinates
fy = fx;
[Fx,Fy]= meshgrid(fx,fy);               % Meshgrid in frequency domain

%% Section 3 - Iterative Phase Retrieval 

iterations = 23400;   % number of iterations
[focal_series, error] = defocus_phase_retrieval_hard(Fx,Fy,z_mm,intensity_focal_series,lambda,iterations); 

in_focus_probe = focal_series(:,:,n+1); % complete probe wavefunction in the infocus plane

% ~~~~~~~~~

% Visualization of recovered probe wavefunction phase
figure(4);
imagesc(x,y,angle(in_focus_probe)); axis image;colormap turbo;colorbar;
title("Recovered phase of infocus electron probe wavefunction"); 

%% Section 4 - Visualizing the PSF Optics using recovered phase and intensity

psf_visualization(X,Y,x,y,in_focus_probe);
