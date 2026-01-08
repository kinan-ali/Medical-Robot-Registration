function targ_pos = GetTargetPosition(flag_noise)

global Tworld_cam Tworld_org torg_target camera_real im_pts

if nargin < 1
	flag_noise = 0;
end

if flag_noise == 0
	noise_ampl = 0;
else
	noise_ampl = 0.5; %1;
end

%Tworld_cam
%Tworld_org
Tcam_org = inv(Tworld_cam)*Tworld_org;

[im_pts, liste_visible_points] = camera_proj(Tcam_org, camera_real, torg_target, noise_ampl);

if liste_visible_points(1) == 0
	fprintf(1, 'point non visible\n');	
end

targ_pos = im_pts(1:2);

