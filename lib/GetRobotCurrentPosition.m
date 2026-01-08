function transf_base_eff = GetRobotCurrentPosition(flag_noise)

global Tbase_eff;

if nargin < 1
	flag_noise = 0;
end

if flag_noise == 0
	transf_base_eff = Tbase_eff;
else
	transf_base_eff = add_noise(Tbase_eff);
	
	
	
end
return;


function Tbase_eff_noisy = add_noise(Tbase_eff)
	dist = norm(Tbase_eff(1:3,4));
	noise_trans = 2*noise_trans_func(dist)*(rand(3,1) - 0.5*ones(3,1));
	noise_rot = 2*noise_rot_func(dist)*(rand(1,1) - 0.5)
	Rbase_eff = Tbase_eff(1:3,1:3);
	tbase_eff = Tbase_eff(1:3,4);
	[theta,u] = r2thetau(Rbase_eff)
	ang_noisy = theta + noise_rot
	Rbase_eff_noisy = thetau2r(ang_noisy*u) 
	tbase_eff_noisy = tbase_eff + noise_trans;
	
	Tbase_eff_noisy = [Rbase_eff_noisy, tbase_eff_noisy; 0 0 0 1]; 

return;



function noise_trans = noise_trans_func(dist)
	
	%noise_trans = 0;
	%noise_trans = dist/10000;
	noise_trans = 0.1; %0.02;

return

function noise_rot = noise_rot_func(dist)
	
	%noise_rot = 0;
	%noise_rot = dist/10000 * pi/180;
	noise_rot = 0.1*pi/180; %0.02*pi/180;

return


