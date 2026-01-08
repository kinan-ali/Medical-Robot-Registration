function MoveEffVelocity(velocity, period)

% Commande en vitesse opérationnelle dans le repère de l'effecteur
% velocity ne contient que les vitesses cartesiennes de l'effecteur
% puisque le poignet est passif

global Tbase_eff tbase_troc_inst joint_position;

if ((size(velocity,1) ~= 3) || (size(velocity,2)~=1))
	fprintf('dimension incorrect du vecteur de vitesse : "velocity" doit être un vecteur colonne de dimension 3\n');
	return;
end

if((size(period,1) ~=1) || (size(period,2) ~= 1))
	fprintf(1, 'period doit être une durée en secondes \n');
	return;
end


fprintf(1, 'Moving effector with velocity : %f %f %f\n', velocity(1), velocity(2), velocity(3));


tbase_eff = Tbase_eff(1:3,4);
Rbase_eff = Tbase_eff(1:3,1:3);

tbase_eff_new = tbase_eff + Rbase_eff*velocity*period;

joint_position = MGI_pos(tbase_eff_new);

%tbase_troc_inst
Tbase_eff = MGD(joint_position, tbase_troc_inst);

%test = (tbase_troc_inst - Tbase_eff(1:3,4))./Tbase_eff(1:3,3)

fprintf(1, 'effector moved\n');