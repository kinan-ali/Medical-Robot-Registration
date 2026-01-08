function DisplayEndoscope()

global Tworld_cam endoscope_length endoscope_radius

plot_camera(0, Tworld_cam, 'r', 1, 0)

scale = endoscope_length;
[X,Y,Z] = cylinder(endoscope_radius/scale*[1;1]);

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

Rworld_cam = Tworld_cam(1:3,1:3);
tworld_cam = Tworld_cam(1:3,4);
endo_moved = Rworld_cam*[X_line;Y_line;Z_line] + repmat(tworld_cam - Rworld_cam*[0;0;endoscope_length], 1, size(X,1)*size(X,2));

X_moved = reshape(endo_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(endo_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(endo_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);