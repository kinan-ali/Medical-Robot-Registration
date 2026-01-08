function DisplayOrgan

global Tworld_org torg_target 

x = [-50:10:50];
y = [-50:10:50];

[X,Y] = meshgrid(x,y);

Z = zeros(length(y), length(x));

X_line = reshape(X, 1, size(X,1)*size(X,2));
Y_line = reshape(Y, 1, size(X,1)*size(X,2));
Z_line = reshape(Z, 1, size(X,1)*size(X,2));

Rworld_org = Tworld_org(1:3,1:3);
tworld_org = Tworld_org(1:3,4);
org_moved = Rworld_org*[X_line;Y_line;Z_line] + repmat(tworld_org, 1, size(X,1)*size(X,2));

X_moved = reshape(org_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(org_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(org_moved(3,:), size(X,1), size(X,2)); 

mesh(X_moved, Y_moved, Z_moved);



scale = 10;
[X,Y,Z] = sphere;

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

tworld_target = Rworld_org * torg_target + tworld_org;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_target, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);
