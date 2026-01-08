function plot_ellipsoid(center, rot_mat, axe_lengthes, scale)

% plot_ellipsoid(center, rot_mat, axe_lengthes, scale)
%
% Trace un ellipsoid sur la figure active, centré en center,
% d'axes principaaux de longueurs axe_lengthes et orienté 
% par la matrice de rotation rot_mat par rapport au repère de 
% la figure. scale est un paramètre optionnel permettant de tracer
% l'ellipsoid plus grand (scale > 1) ou plus petit (scale < 1) 
% selon les besoins d'affichage.

if nargin < 4
	scale = 1;
end


	[X,Y,Z] = ellipsoid(0,0,0, scale*axe_lengthes(1), scale*axe_lengthes(2), scale*axe_lengthes(3));
	
	Xline = reshape(X, 1, size(X,1)*size(X,2));
	Yline = reshape(Y, 1, size(X,1)*size(X,2));
	Zline = reshape(Z, 1, size(X,1)*size(X,2));

	ell_line = [Xline; Yline; Zline];

	ell_mod = rot_mat*ell_line + repmat(center, 1, size(X,1)*size(X,2));
	
	Xmod = reshape(ell_mod(1,:), size(X,1), size(X,2));
	Ymod = reshape(ell_mod(2,:), size(X,1), size(X,2));
	Zmod = reshape(ell_mod(3,:), size(X,1), size(X,2));
	
	
	hold on;
	surf(Xmod, Ymod, Zmod);