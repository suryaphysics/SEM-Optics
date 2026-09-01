function [focal_series,error,U_center] = defocus_phase_retrieval( ...
    X,Y,z,focal_series_intensity,lambda,iter,in_focus_probe)

% Stable defocus-diversity phase retrieval
%
% Inputs:
%   X,Y                   : spatial coordinate grids
%   z                     : defocus positions
%   focal_series_intensity: measured intensity stack
%   lambda                : wavelength
%   iter                  : number of iterations
%
% Outputs:
%   focal_series          : recovered complex fields
%   error                 : normalized amplitude error
%   U_center              : recovered center-plane field


%% Initialization

num_planes = size(focal_series_intensity,3);
center_idx = ceil(num_planes/2);

A_center = sqrt(max(focal_series_intensity(:,:,center_idx),0));

focal_series = zeros(size(focal_series_intensity),'like',1i);
error = zeros(1,iter);

% Initial center-plane phase estimate
% U_center = A_center .* exp(1i);
U_center = in_focus_probe;


%% Projection parameters

gamma = 1;
amin = 1e-8;


%% Iterative reconstruction

for it = 80001:iter

    focal_series(:,:,center_idx) = U_center;


    % Center -> last plane
    for k = center_idx+1:num_planes

        dz = z(k)-z(k-1);

        U = ASM_propagation(focal_series(:,:,k-1),dz,X,Y,lambda);
        U = sanitize_field(U);

        focal_series(:,:,k) = relaxed_amplitude_projection( ...
            U,focal_series_intensity(:,:,k),gamma,amin);

    end


    % Last plane -> first plane
    for k = num_planes-1:-1:1

        dz = z(k)-z(k+1);

        U = ASM_propagation(focal_series(:,:,k+1),dz,X,Y,lambda);
        U = sanitize_field(U);

        focal_series(:,:,k) = relaxed_amplitude_projection( ...
            U,focal_series_intensity(:,:,k),gamma,amin);

    end


    % First plane -> center
    for k = 2:center_idx

        dz = z(k)-z(k-1);

        U = ASM_propagation(focal_series(:,:,k-1),dz,X,Y,lambda);
        U = sanitize_field(U);

        focal_series(:,:,k) = relaxed_amplitude_projection( ...
            U,focal_series_intensity(:,:,k),gamma,amin);

    end


    % Update center estimate
    U_center = focal_series(:,:,center_idx);


    % Convergence error
    A_est = abs(U_center);

    error(it) = sum((A_center(:)-A_est(:)).^2) / ...
        (sum(A_center(:).^2)+eps);


    % fprintf('Iteration %d/%d | Error %.6e | Max amplitude %.3e\n',...
    %     it,iter,error(it),max(A_est(:)));
    % 
    % % Display progress
    % fprintf('Iter %d / %d | Error = %.6e | Max amplitude = %.3e\n', ...
    %     it, iter, error(it), max(abs(U_center(:))));

    % Update phase display every 10 iterations
    if mod(it,10) == 0 || it == 1 || it == iter

        figure(4);
        imagesc(angle(U_center));
        axis image;
        colormap parula;
        colorbar;
        title(sprintf('Recovered phase (Iteration %d)', it));

        drawnow;

    end
end

end


%% Relaxed amplitude projection

function Uout = relaxed_amplitude_projection(U,Imeas,gamma,amin)

Ameas = sqrt(max(Imeas,0));
Ameas = max(Ameas,amin);

Aold = abs(U);

phaseU = angle(U + eps);
phaseU(~isfinite(phaseU)) = 0;

Anew = (1-gamma).*Aold + gamma.*Ameas;

Uout = sanitize_field(Anew .* exp(1i*phaseU));

end


%% Remove invalid numerical values

function U = sanitize_field(U)

U(~isfinite(real(U)) | ~isfinite(imag(U))) = 0;

max_allowed = 1e12;

amp = abs(U);
bad = amp > max_allowed;

U(bad) = max_allowed .* exp(1i*angle(U(bad)));

end