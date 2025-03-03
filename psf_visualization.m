function [] = psf_visualization(X,Y,x,y,probe)

% Function to visualize the PSF of SEM Optics based on recovered probe wavefunction

d_g = 2.1e-09;  % Demagnified ideal probe width parameter for the infocus plane probe
psi_g = exp(-((X.^2 + Y.^2)/(d_g)^2));    % Ideal Gaussian probe at the sample, demagnified geometrical image of the virtual source.

PSI_g = fftshift(fft2(ifftshift(psi_g)));  % Geometrical probe (Gaussian) Fourier spectrum
PROBE = fftshift(fft2(ifftshift((probe)))); % Experimentally recovered probe Fourier spectrum 

H = PROBE./PSI_g; % Inverse function
psf = fftshift(ifft2(ifftshift(H))); % PSF- the blurring function

% ~~~~~~~~~

% Visualization of PSF Optics

psf_optics = imresize(abs(psf),[1024,1024]);
figure(5);imagesc(x,y,psf_optics); axis image; colormap turbo
title("Point Spread function of SEM Optics - |PSF_{optics}|");
xlabel('m')
ylabel('m')
set(gca,'Fontsize',18);

