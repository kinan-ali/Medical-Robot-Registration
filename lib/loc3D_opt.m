function [Rsensor_obj, tsensor_obj] = loc3D_opt(pts_in_sensor, pts_in_obj)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of eigen frame from object information
Og = mean(pts_in_obj, 2);

pts_obj_in_eig = pts_in_obj - Og;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of eigen frame from sensor information
Sg = mean(pts_in_sensor, 2);

pts_sensor_in_eig = pts_in_sensor - Sg;

% Covariance matrix
M = pts_sensor_in_eig * pts_obj_in_eig';

[U, S, V] = svd(M);

if det(U*V') > 0
    Rsensor_obj = U*V';
else
    Rsensor_obj = U*[1 0 0; 0 1 0; 0 0 -1]*V';
end

tsensor_obj = Sg - Rsensor_obj*Og;

%Tsensor_obj = [Rsensor_obj, tsensor_obj; 0 0 0 1];

end

