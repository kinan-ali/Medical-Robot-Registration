function [T] = loc_3Dhorn(pts_frame1, pts_frame2)
% Function to solve the 3D localization using the Horn method
% Input: points in the first frame and correspondent points in the second
% frame
% Output: the transform frame1_T_frame2

% Centres of mass
G1 = mean(pts_frame1,2);
G2 = mean(pts_frame2,2);

% Recenter points in the centre of mass
P_g1 = pts_frame1 - G1;
P_g2 = pts_frame2 - G2;

% Compute M matrix
M = P_g1 * (P_g2');

% Find optimal rotation
[U,~,V] = svd(M);
R = U * (V');

% Find optimal translation
t = G1 - (R*G2);

% Compose T
T = [R t; [0 0 0 1]];

end