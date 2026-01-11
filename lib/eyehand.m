function [X] = eyehand(listA, listB)
% Function to solve the eye-hand calibration.
% Inputs: ordered list of A and B transforms as three dimensional matrices
% Output: fixed X transform

% Find the axis of rotation for each A and B transformation
nb_configurations = size(listA,3); 
list_axisA = [];
list_axisB = [];

for ii = 1:nb_configurations
    [~, axisA] = r2thetau(listA(1:3,1:3,ii));
    list_axisA = [list_axisA, axisA];
end

for ii = 1:nb_configurations
    [~, axisB] = r2thetau(listB(1:3,1:3,ii));
    list_axisB = [list_axisB, axisB];
end

% Use 3D localization to find Rx
T = loc_3Dhorn(list_axisA, list_axisB);
Rx = T(1:3, 1:3);

% Find t through minimization of a cost function (linear solving)
for ii = 1:nb_configurations
    k = (ii-1)*3 +1;
    A(k:k+2,1:3) = listA(1:3,1:3,ii) - eye(3);
    b(k:k+2,1) = Rx * listB(1:3,4,ii) - listA(1:3,4,ii);
end

tx = pinv(A)*b;
X = [Rx tx; 0 0 0 1];

end


% NOTE on solving AX = XB for translation
% - 1) take AX = XB and split in two equations, one for rotation
% (solved before) and one for translation: Ra*tx + ta = Rx*tb + tx
% - 2) Group for tx: (Ra - I) tx = Rx*tb - ta which is in the form Ax = b
% The structure of the total system becomes: 
% [Ra1 - I;         [Rx*tb1 - ta1;
%  Ra2 - I;  * tx =  Rx*tb2 - ta2;
%  ...               ...
%  RaN - I]          Rx*tbN - taN]
% Therefore each movement contributes with a 3 equations block
% - 3) Solving with Moore-Penrose pseudoinverse as it finds tx by
% minimizing the norm of the residual: min||A*tx - b||^2 (see cost
% function)