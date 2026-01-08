function R = theta2r( theta )
%
% Christophe Doignon (LSiiT-GRAViR)
%
% Le 21 Fevrier 2001
%
% Creation de la matrice de rotation a partir de trois angles fixes
%
R      = zeros(3,3);
R(1,1) =  cos(theta(2))*cos(theta(3));
R(1,2) = -cos(theta(1))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3));
R(1,3) =  sin(theta(1))*sin(theta(3)) + cos(theta(1))*sin(theta(2))*cos(theta(3));
R(2,1) =  cos(theta(2))*sin(theta(3));
R(2,2) =  cos(theta(1))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3));
R(2,3) = -sin(theta(1))*cos(theta(3)) + cos(theta(1))*sin(theta(2))*sin(theta(3));
R(3,1) = -sin(theta(2));
R(3,2) =  sin(theta(1))*cos(theta(2));
R(3,3) =  cos(theta(1))*cos(theta(2));