function [focal_series,error,temp] = defocus_phase_retrieval_hard(X,Y,z,focal_series_intensity,lambda,iter)

I_center = sqrt(focal_series_intensity(:,:,10));
focal_series = zeros(size(focal_series_intensity));
error = zeros(1,iter);

%% Phase initializations
focal_series(:,:,10) = I_center.*exp(1i*2*pi*ones(size(I_center))); % constant initialization

for i = 1:iter

    focal_series(:,:,11) = ASM_propagation(  focal_series(:,:,10),   (z(11)-z(10))   ,X,Y,lambda);
    focal_series(:,:,11) =  sqrt(focal_series_intensity(:,:,10)).*exp(1i*angle(focal_series(:,:,11)));

    focal_series(:,:,12) = ASM_propagation(  focal_series(:,:,11),   (z(12)-z(11))   ,X,Y,lambda);
    focal_series(:,:,12) =  sqrt(focal_series_intensity(:,:,12)).*exp(1i*angle(focal_series(:,:,12)));

    focal_series(:,:,13) = ASM_propagation(  focal_series(:,:,12),   (z(13)-z(12))   ,X,Y,lambda);
    focal_series(:,:,13) =  sqrt(focal_series_intensity(:,:,13)).*exp(1i*angle(focal_series(:,:,13)));

    focal_series(:,:,14) = ASM_propagation(  focal_series(:,:,13),  (z(14)-z(13))  ,X,Y,lambda);
    focal_series(:,:,14) =  sqrt(focal_series_intensity(:,:,14)).*exp(1i*angle(focal_series(:,:,14)));

    focal_series(:,:,15) = ASM_propagation(  focal_series(:,:,14),   (z(15)-z(14))   ,X,Y,lambda);
    focal_series(:,:,15) =  sqrt(focal_series_intensity(:,:,15)).*exp(1i*angle(focal_series(:,:,15)));

    focal_series(:,:,16) = ASM_propagation(  focal_series(:,:,15),   (z(16)-z(15))   ,X,Y,lambda);
    focal_series(:,:,16) =  sqrt(focal_series_intensity(:,:,16)).*exp(1i*angle(focal_series(:,:,16)));

    focal_series(:,:,17) = ASM_propagation(  focal_series(:,:,16),   (z(17)-z(16))   ,X,Y,lambda);
    focal_series(:,:,17) =  sqrt(focal_series_intensity(:,:,17)).*exp(1i*angle(focal_series(:,:,17)));

    focal_series(:,:,18) = ASM_propagation(  focal_series(:,:,17),   (z(18)-z(17))   ,X,Y,lambda);
    focal_series(:,:,18) =  sqrt(focal_series_intensity(:,:,18)).*exp(1i*angle(focal_series(:,:,18)));

    focal_series(:,:,19) = ASM_propagation(  focal_series(:,:,18),   (z(19)-z(18))   ,X,Y,lambda);
    focal_series(:,:,19) =  sqrt(focal_series_intensity(:,:,19)).*exp(1i*angle(focal_series(:,:,19)));
    %%
    focal_series(:,:,18) = ASM_propagation(  focal_series(:,:,19),   -(z(19)-z(18))   ,X,Y,lambda);
    focal_series(:,:,18) =  sqrt(focal_series_intensity(:,:,18)).*exp(1i*angle(focal_series(:,:,18)));

    focal_series(:,:,17) = ASM_propagation(  focal_series(:,:,18),   -(z(18)-z(17))   ,X,Y,lambda);
    focal_series(:,:,17) =  sqrt(focal_series_intensity(:,:,17)).*exp(1i*angle(focal_series(:,:,17)));

    focal_series(:,:,16) = ASM_propagation(  focal_series(:,:,17),   -(z(17)-z(16))   ,X,Y,lambda);
    focal_series(:,:,16) =  sqrt(focal_series_intensity(:,:,16)).*exp(1i*angle(focal_series(:,:,16)));
    %
    focal_series(:,:,15) = ASM_propagation(  focal_series(:,:,16),   -(z(16)-z(15))   ,X,Y,lambda);
    focal_series(:,:,15) =  sqrt(focal_series_intensity(:,:,15)).*exp(1i*angle(focal_series(:,:,15)));

    focal_series(:,:,14) = ASM_propagation(  focal_series(:,:,15),   -(z(15)-z(14))   ,X,Y,lambda);
    focal_series(:,:,14) =  sqrt(focal_series_intensity(:,:,14)).*exp(1i*angle(focal_series(:,:,14)));

    focal_series(:,:,13) = ASM_propagation(  focal_series(:,:,14),   -(z(14)-z(13))   ,X,Y,lambda);
    focal_series(:,:,13) =  sqrt(focal_series_intensity(:,:,13)).*exp(1i*angle(focal_series(:,:,13)));

    focal_series(:,:,12) = ASM_propagation(  focal_series(:,:,13),   -(z(13)-z(12))   ,X,Y,lambda);
    focal_series(:,:,12) =  sqrt(focal_series_intensity(:,:,12)).*exp(1i*angle(focal_series(:,:,12)));

    focal_series(:,:,11) = ASM_propagation(  focal_series(:,:,12),   -(z(12)-z(11))   ,X,Y,lambda);
    focal_series(:,:,11) =  sqrt(focal_series_intensity(:,:,11)).*exp(1i*angle(focal_series(:,:,11)));

    focal_series(:,:,10) = ASM_propagation(  focal_series(:,:,11),   -(z(11)-z(10))   ,X,Y,lambda);
    focal_series(:,:,10) =  sqrt(focal_series_intensity(:,:,10)).*exp(1i*angle(focal_series(:,:,10)));


    focal_series(:,:,9) = ASM_propagation(  focal_series(:,:,10),   -(z(10)-z(9))   ,X,Y,lambda);
    focal_series(:,:,9) =  sqrt(focal_series_intensity(:,:,9)).*exp(1i*angle(focal_series(:,:,9)));


    focal_series(:,:,8) = ASM_propagation(  focal_series(:,:,9),   -(z(9)-z(8))   ,X,Y,lambda);
    focal_series(:,:,8) =  sqrt(focal_series_intensity(:,:,8)).*exp(1i*angle(focal_series(:,:,8)));


    focal_series(:,:,7) = ASM_propagation(  focal_series(:,:,8),   -(z(8)-z(7))   ,X,Y,lambda);
    focal_series(:,:,7) =  sqrt(focal_series_intensity(:,:,7)).*exp(1i*angle(focal_series(:,:,7)));


    focal_series(:,:,6) = ASM_propagation(  focal_series(:,:,7),   -(z(7)-z(6))   ,X,Y,lambda);
    focal_series(:,:,6) =  sqrt(focal_series_intensity(:,:,6)).*exp(1i*angle(focal_series(:,:,6)));

    focal_series(:,:,5) = ASM_propagation(  focal_series(:,:,6),   -(z(6)-z(5))   ,X,Y,lambda);
    focal_series(:,:,5) =  sqrt(focal_series_intensity(:,:,5)).*exp(1i*angle(focal_series(:,:,5)));


    focal_series(:,:,4) = ASM_propagation(  focal_series(:,:,5),   -(z(5)-z(4))   ,X,Y,lambda);
    focal_series(:,:,4) =  sqrt(focal_series_intensity(:,:,4)).*exp(1i*angle(focal_series(:,:,4)));


    focal_series(:,:,3) = ASM_propagation(  focal_series(:,:,4),   -(z(4)-z(3))   ,X,Y,lambda);
    focal_series(:,:,3) =  sqrt(focal_series_intensity(:,:,3)).*exp(1i*angle(focal_series(:,:,3)));


    focal_series(:,:,2) = ASM_propagation(  focal_series(:,:,3),   -(z(3)-z(2))   ,X,Y,lambda);
    focal_series(:,:,2) =  sqrt(focal_series_intensity(:,:,2)).*exp(1i*angle(focal_series(:,:,2)));


    focal_series(:,:,1) = ASM_propagation(  focal_series(:,:,2),   -(z(2)-z(1))   ,X,Y,lambda);
    focal_series(:,:,1) =  sqrt(focal_series_intensity(:,:,1)).*exp(1i*angle(focal_series(:,:,1)));
    %%
    focal_series(:,:,2) = ASM_propagation(  focal_series(:,:,1),   (z(2)-z(1))   ,X,Y,lambda);
    focal_series(:,:,2) =  sqrt(focal_series_intensity(:,:,2)).*exp(1i*angle(focal_series(:,:,2)));

    focal_series(:,:,3) = ASM_propagation(  focal_series(:,:,2),   (z(3)-z(2))   ,X,Y,lambda);
    focal_series(:,:,3) =  sqrt(focal_series_intensity(:,:,3)).*exp(1i*angle(focal_series(:,:,3)));


    focal_series(:,:,4) = ASM_propagation(  focal_series(:,:,3),   (z(4)-z(3))   ,X,Y,lambda);
    focal_series(:,:,4) =  sqrt(focal_series_intensity(:,:,4)).*exp(1i*angle(focal_series(:,:,4)));

    focal_series(:,:,5) = ASM_propagation(  focal_series(:,:,4),   (z(5)-z(4))   ,X,Y,lambda);
    focal_series(:,:,5) =  sqrt(focal_series_intensity(:,:,5)).*exp(1i*angle(focal_series(:,:,5)));

    focal_series(:,:,6) = ASM_propagation(  focal_series(:,:,5),   (z(6)-z(5))   ,X,Y,lambda);
    focal_series(:,:,6) =  sqrt(focal_series_intensity(:,:,6)).*exp(1i*angle(focal_series(:,:,6)));

    focal_series(:,:,7) = ASM_propagation(  focal_series(:,:,6),   (z(7)-z(6))   ,X,Y,lambda);
    focal_series(:,:,7) =  sqrt(focal_series_intensity(:,:,7)).*exp(1i*angle(focal_series(:,:,7)));

    focal_series(:,:,8) = ASM_propagation(  focal_series(:,:,7),   (z(8)-z(7))   ,X,Y,lambda);
    focal_series(:,:,8) =  sqrt(focal_series_intensity(:,:,8)).*exp(1i*angle(focal_series(:,:,8)));

    focal_series(:,:,9) = ASM_propagation(  focal_series(:,:,8),   (z(9)-z(8))   ,X,Y,lambda);
    focal_series(:,:,9) =  sqrt(focal_series_intensity(:,:,9)).*exp(1i*angle(focal_series(:,:,9)));

    focal_series(:,:,10) = ASM_propagation(  focal_series(:,:,9),   (z(10)-z(9))   ,X,Y,lambda);
    temp = abs(focal_series(:,:,10));
    focal_series(:,:,10) =  sqrt(focal_series_intensity(:,:,10)).*exp(1i*angle(focal_series(:,:,10)));
    %
    error_temp =  (abs(I_center) - abs(temp)).^2;
    error(1,i) = sum(error_temp(:))/sum(abs(I_center(:).^2));
    %
    i
    % 
    % figure(4);
    % imagesc(angle(focal_series(:,:,10))); axis image;colormap parula
    % figure(5)
    % plot(1:iter,(error));
end
end