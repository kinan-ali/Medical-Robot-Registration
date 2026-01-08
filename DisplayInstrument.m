function DisplayInstrument()

global Tworld_base Tbase_eff Teff_inst instrument_radius instrument_length

Tworld_inst = Tworld_base * Tbase_eff * Teff_inst;

scale = instrument_length;
[X,Y,Z] = cylinder(instrument_radius/scale*[1;1]);

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

Rworld_inst = Tworld_inst(1:3,1:3);
tworld_inst = Tworld_inst(1:3,4);
inst_moved = Rworld_inst*[X_line;Y_line;Z_line] + repmat(tworld_inst - Rworld_inst*[0;0;instrument_length], 1, size(X,1)*size(X,2));

X_moved = reshape(inst_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(inst_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(inst_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);