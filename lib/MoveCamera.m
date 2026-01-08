function MoveCamera(trans, thetaz)

global Tworld_cam tworld_troc_cam

fprintf(1, 'Moving camera : %f %f %f %f\n', trans(1), trans(2), trans(3), thetaz);

% Les arguments d'entrée sont des indications, mais
% le mouvement n'est pas diretement contrôle, donc incertain
trans_app = trans + 2*(rand(3,1)-0.5*ones(3,1));

thetaz_app = thetaz*pi/180 + 20*pi/180*(rand(1,1) - 0.5);

Rworld_cam = Tworld_cam(1:3,1:3);
tworld_cam = Tworld_cam(1:3,4);

tworld_cam_new = tworld_cam + Rworld_cam*trans_app;

% axe de l'endoscope dans repere monde
axis_in_world = tworld_cam - tworld_troc_cam;
axis_in_world = axis_in_world/(norm(axis_in_world));

axis_in_world_new = tworld_cam_new - tworld_troc_cam;
axis_in_world_new = axis_in_world_new/(norm(axis_in_world_new));

axis_rot = cross(axis_in_world, axis_in_world_new);
sang = norm(axis_rot);
ang = asin(sang);

if abs(sang) > 1e-5
	axis_rot = axis_rot/sang; 
else
	axis_rot = [0;0;1];
end

axis_rot_in_cam = Rworld_cam' * axis_rot;

Rcam_cam_new = thetau2r(ang*axis_rot_in_cam) * theta2r([0;0;thetaz_app]);

Rworld_cam_new = Rworld_cam * Rcam_cam_new;


%tworld_ccd = tworld_cam_new - length_endoscope*axis_in_world;
%
%tworld_eff_rob_cam = tworld_ccd 
%
%Tbase_rob_cam_world = inv(Tworld_base_rob_cam);
%tbase_rob_cam_world = Tbase_rob_cam_world(1:3,4);
%
%
%
%tbase_eff_rob_cam = tbase_rob_cam_world + Rbase_rob_cam_world * (tworld_cam_new - length_cam* Rbse_rob_caminv(Tworld_base_rob_cam) * 



Tworld_cam = [Rworld_cam_new tworld_cam_new; 0 0 0 1];

fprintf(1, 'camera moved\n')