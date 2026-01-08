function MoveEffPosition(position)


% la position est donnee par rapport à la base du robot
% On ne donne pas l'orientation (passive)

global Tbase_eff tbase_troc_inst joint_position;

fprintf(1, 'Moving effector to : %f %f %f\n', position(1), position(2), position(3));

if ((size(position,1) ~= 3) | (size(position,2)~=1))
	fprintf(2, 'error : position must be given by a 3D vector\n');
	return;
end

%position

joint_position = MGI_pos(position);

%tbase_troc_inst
Tbase_eff = MGD(joint_position, tbase_troc_inst);

%test = (tbase_troc_inst - Tbase_eff(1:3,4))./Tbase_eff(1:3,3)

fprintf(1, 'effector moved\n');