function Tcam_inst_est = GetInstrumentPosition(flag_noise)

% Retourne la position de l'instrument dans le repère de
% la caméra endoscopique dans le cas de l'utilisation
% de marqueurs visuels sur l'instrument


global Tworld_loc Tworld_cam Tworld_base Tbase_eff Teff_inst Tcam_mark_cam Tinst_mark_inst camera_real mark_inst im_pts_mark_inst

if nargin < 1
	flag_noise = 0;
end

fprintf(1, 'obtaining position of instrument in camera frame\n');

Tcam_inst = inv(Tworld_cam)*Tworld_base*Tbase_eff*Teff_inst;

% On verifie si le marqueur de l'instrument (bout et un point 
%intermediaire) sont dans l'image
[im_pts_mark_inst, liste_visible_points] = camera_proj(Tcam_inst, camera_real, mark_inst, 0);


if((liste_visible_points(1) == 0) || (liste_visible_points(2) == 0))
	Tcam_inst_est = [];
	return;

else
	Tcam_inst_est = Tcam_inst;
end


%% On enlève la rotation autour de l'axe de l'instrument
%Rcam_inst = Tcam_inst(1:3,1:3);
%
%theta = r2eulerxyz(Rcam_inst); 
%
%Rcam_inst_est = theta2r([theta(1);0;0]) * theta2r([0;theta(2);0]);
%
%Tcam_inst_est = [Rcam_inst_est Tcam_inst(1:3,4); 0 0 0 1];

if(flag_noise ~= 0)
    %Tcam_inst_est
	Tcam_inst_est = add_noise(Tcam_inst_est);
    %Tcam_inst_est
    %pause;
end
	
fprintf(1, 'instrument positions obtained\n');

return;


function Tloc_mark_noisy = add_noise(Tloc_mark)
	dist = norm(Tloc_mark(1:3,4));
	noise_trans = 2*noise_trans_func(dist)*(rand(3,1) - 0.5*ones(3,1));
	noise_rot = 2*noise_rot_func(dist)*(rand(1,1) - 0.5);
	Rloc_mark = Tloc_mark(1:3,1:3);
	tloc_mark = Tloc_mark(1:3,4)
	[theta,u] = r2thetau(Rloc_mark);
	ang_noisy = theta + noise_rot;
	Rloc_mark_noisy = thetau2r(ang_noisy*u); 
	tloc_mark_noisy = tloc_mark + noise_trans;
	
	Tloc_mark_noisy = [Rloc_mark_noisy, tloc_mark_noisy; 0 0 0 1]; 

	noisy_minus_ex = Rloc_mark_noisy - Rloc_mark;
return;



function noise_trans = noise_trans_func(dist)
	
	%noise_trans = 0;
	%noise_trans = dist/10000;
	noise_trans = 1; %0.02;

return

function noise_rot = noise_rot_func(dist)
	
	%noise_rot = 0;
	%noise_rot = dist/10000 * pi/180;
	noise_rot = 0.5*pi/180; %0.02*pi/180;

return


