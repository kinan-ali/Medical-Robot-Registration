function mesure = GetLocalizerInformation(flag_noise)

% Retourne la position mesuree des marqueurs


global Tworld_loc Tworld_cam Tworld_base Tbase_eff Teff_inst Tcam_mark_cam TIM %Tinst_mark_inst

if nargin < 1
	flag_noise = 0;
end

fprintf(1, 'obtaining localizer information\n');


Tloc_mark_cam = inv(Tworld_loc)*Tworld_cam * Tcam_mark_cam;

Tloc_mark_inst = inv(Tworld_loc)*Tworld_base * Tbase_eff * Teff_inst * TIM; %Tinst_mark_inst;

if(flag_noise ~= 0)
	mesure.mark(1).T = add_noise(Tloc_mark_cam);
	mesure.mark(2).T = add_noise(Tloc_mark_inst);
else
	mesure.mark(1).T = Tloc_mark_cam;
	mesure.mark(2).T = Tloc_mark_inst;
end
	
fprintf(1, 'localizer information obtained\n');

return;


function Tloc_mark_noisy = add_noise(Tloc_mark)
	dist = norm(Tloc_mark(1:3,4));
	noise_trans = 2*noise_trans_func(dist)*(rand(3,1) - 0.5*ones(3,1));
	noise_rot = 2*noise_rot_func(dist)*(rand(1,1) - 0.5);
	Rloc_mark = Tloc_mark(1:3,1:3);
	tloc_mark = Tloc_mark(1:3,4);
	[theta,u] = r2thetau(Rloc_mark);
	ang_noisy = theta + noise_rot;
	Rloc_mark_noisy = thetau2r(ang_noisy*u); 
	tloc_mark_noisy = tloc_mark + noise_trans;
	
	Tloc_mark_noisy = [Rloc_mark_noisy, tloc_mark_noisy; 0 0 0 1]; 

	noisy_minus_ex = Rloc_mark_noisy - Rloc_mark;
return;


return;

function noise_trans = noise_trans_func(dist)
	
	%noise_trans = 0;
	%noise_trans = dist/10000;
	noise_trans = 0.1;

return

function noise_rot = noise_rot_func(dist)
	
	%noise_rot = 0;
	%noise_rot = dist/10000 * pi/180;
	noise_rot = 0.2*pi/180;

return


