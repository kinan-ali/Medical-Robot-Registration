function DisplayPatient

global Tworld_patient patient_radius tworld_troc_cam tworld_troc_inst

patient_length = 500;

scale = patient_length;
[X,Y,Z] = cylinder(patient_radius/scale*[1;1]);

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

Rworld_patient = Tworld_patient(1:3,1:3);
tworld_patient = Tworld_patient(1:3,4);
patient_moved = Rworld_patient*[X_line;Y_line;Z_line] + repmat(tworld_patient - Rworld_patient*[0;0;patient_length/2], 1, size(X,1)*size(X,2));

X_moved = reshape(patient_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(patient_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(patient_moved(3,:), size(X,1), size(X,2)); 

mesh(X_moved, Y_moved, Z_moved);


scale = 10;
[X,Y,Z] = sphere;

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_troc_cam, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_troc_inst, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);
