function [focal_series, error] = defocus_phase_retrieval(X, Y, z, focal_series_intensity, lambda, iter)

% Function to perform iterative phase retrieval based on defocus variation

% Inputs:
% X, Y: grid coordinates
% z: array of the positions of the planes
% focal_series_intensity: the intensity data for the focal series (3D matrix)
% lambda: wavelength
% iter: number of iterations

% Outputs:
% focal_series: the complete data (recovered phase and intensity) for the focal series (3D matrix)
% error: SSE over the iterations

% Number of planes (total number of focal planes - ALWAYS odd)
num_planes = size(focal_series_intensity, 3);

% Central plane for odd number of planes
center_plane = (num_planes + 1)/2;

% Initialize the center intensity for the focal series
I_center = sqrt(focal_series_intensity(:,:,center_plane));  % Start from the center plane
focal_series = zeros(size(focal_series_intensity));

error = zeros(1, iter);

% Phase initializations
focal_series(:,:,center_plane) = I_center .* exp(1i * 2 * pi * ones(size(I_center))); % constant initialization

for i = 1:iter
    %% Forward Propagation: from center plane to last plane
    for j = center_plane:num_planes-1
        focal_series(:,:,j+1) = ASM_propagation(focal_series(:,:,j), (z(j+1) - z(j)), X, Y, lambda);
        focal_series(:,:,j+1) = sqrt(focal_series_intensity(:,:,j+1)) .* exp(1i * angle(focal_series(:,:,j+1)));
    end

    %% Backward Propagation: from last plane back to plane 1
    for j = num_planes:-1:2
        focal_series(:,:,j-1) = ASM_propagation(focal_series(:,:,j), -(z(j) - z(j-1)), X, Y, lambda);
        focal_series(:,:,j-1) = sqrt(focal_series_intensity(:,:,j-1)) .* exp(1i * angle(focal_series(:,:,j-1)));
    end

    %% Final Forward Propagation: from plane 1 back to center plane
    focal_series(:,:,2) = ASM_propagation(focal_series(:,:,1), (z(2) - z(1)), X, Y, lambda);
    focal_series(:,:,2) = sqrt(focal_series_intensity(:,:,2)) .* exp(1i * angle(focal_series(:,:,2)));

    for j = 3:center_plane
        focal_series(:,:,j) = ASM_propagation(focal_series(:,:,j-1), (z(j) - z(j-1)), X, Y, lambda);
        if j==center_plane
            temp = abs(focal_series(:,:,j)); % amplitude at the center plane after the iteration before it is replaced by known probe intensity
        end
        focal_series(:,:,j) = sqrt(focal_series_intensity(:,:,j)) .* exp(1i * angle(focal_series(:,:,j)));
    end

    %% Calculate error (SSE) after the full cycle (one interation)
    error_temp = (abs(I_center) - abs(temp)).^2;
    error(1, i) = sum(error_temp(:)) / sum(abs(I_center(:)).^2);
  
  % Optionally, display the phase angle of the central image after the cycle
    figure(2);
    imagesc(angle(focal_series(:,:,center_plane)));axis image; colormap turbo;colorbar;
    title("Iteration "+i);
     
    % Optionally, you can also plot the error over iterations
    % figure(3)
    % plot(1:i, error(1:i));
end
end
