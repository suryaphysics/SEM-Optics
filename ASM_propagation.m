function [psi_out] = ASM_propagation(psi_in,z,Fx,Fy,lambda)

% Propagation function based on the angular spectrum method

k = 2*pi/lambda;
% Transfer Function for Propagation
    H_ASM = exp(1i* k * z * sqrt(1 - (lambda*Fx).^2 - (lambda*Fy).^2 ) )...;
        .*(sqrt((lambda*Fx).^2 + (lambda*Fy).^2)<(1/lambda));

    angular_spectrum_phi_in = fftshift(fft2(ifftshift(psi_in)));

    angular_spectrum_phi_out = angular_spectrum_phi_in.*H_ASM;

    psi_out = (fftshift(ifft2(ifftshift(angular_spectrum_phi_out))));

end