function DisplayLocalizer

global Tworld_loc Tworld_cam Tcam_mark_cam Tworld_base Tbase_eff Teff_inst TIM %Tinst_mark_inst

Pt(:,1) = [-250; -50; 0];
Pt(:,2) = [250; -50; 0];
Pt(:,3) = [250; 50; 0];
Pt(:,4) = [-250; 50; 0];

Pt(:,5) = [-250; -50; -100];
Pt(:,6) = [250; -50; -100];
Pt(:,7) = [250; 50; -100];
Pt(:,8) = [-250; 50; -100];

Rworld_loc = Tworld_loc(1:3,1:3);
tworld_loc = Tworld_loc(1:3,4); 
Pt_moved = Rworld_loc*Pt + repmat(tworld_loc, 1, 8);

plot3([Pt_moved(1,1) Pt_moved(1,2)], [Pt_moved(2,1) Pt_moved(2,2)], [Pt_moved(3,1) Pt_moved(3,2)]);
plot3([Pt_moved(1,2) Pt_moved(1,3)], [Pt_moved(2,2) Pt_moved(2,3)], [Pt_moved(3,2) Pt_moved(3,3)]);
plot3([Pt_moved(1,3) Pt_moved(1,4)], [Pt_moved(2,3) Pt_moved(2,4)], [Pt_moved(3,3) Pt_moved(3,4)]);
plot3([Pt_moved(1,4) Pt_moved(1,1)], [Pt_moved(2,4) Pt_moved(2,1)], [Pt_moved(3,4) Pt_moved(3,1)]);

plot3([Pt_moved(1,5) Pt_moved(1,6)], [Pt_moved(2,5) Pt_moved(2,6)], [Pt_moved(3,5) Pt_moved(3,6)]);
plot3([Pt_moved(1,6) Pt_moved(1,7)], [Pt_moved(2,6) Pt_moved(2,7)], [Pt_moved(3,6) Pt_moved(3,7)]);
plot3([Pt_moved(1,7) Pt_moved(1,8)], [Pt_moved(2,7) Pt_moved(2,8)], [Pt_moved(3,7) Pt_moved(3,8)]);
plot3([Pt_moved(1,8) Pt_moved(1,5)], [Pt_moved(2,8) Pt_moved(2,5)], [Pt_moved(3,8) Pt_moved(3,5)]);

plot3([Pt_moved(1,1) Pt_moved(1,5)], [Pt_moved(2,1) Pt_moved(2,5)], [Pt_moved(3,1) Pt_moved(3,5)]);
plot3([Pt_moved(1,2) Pt_moved(1,6)], [Pt_moved(2,2) Pt_moved(2,6)], [Pt_moved(3,2) Pt_moved(3,6)]);
plot3([Pt_moved(1,3) Pt_moved(1,7)], [Pt_moved(2,3) Pt_moved(2,7)], [Pt_moved(3,3) Pt_moved(3,7)]);
plot3([Pt_moved(1,4) Pt_moved(1,8)], [Pt_moved(2,4) Pt_moved(2,8)], [Pt_moved(3,4) Pt_moved(3,8)]);
 

scale = 20;
[X,Y,Z] = sphere;

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

tworld_sphere = Rworld_loc * [-200;0;0] + tworld_loc;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);



X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

tworld_sphere = Rworld_loc * [200;0;0] + tworld_loc;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);

% Marqueurs
Pt(:,1) = [-25; -25; 0];
Pt(:,2) = [25; -25; 0];
Pt(:,3) = [25; 25; 0];
Pt(:,4) = [-25; 25; 0];

Tworld_mark = Tworld_cam * Tcam_mark_cam; 
Rworld_mark = Tworld_mark(1:3,1:3);
tworld_mark = Tworld_mark(1:3,4); 
Pt_moved = Rworld_mark*Pt + repmat(tworld_mark, 1, 8);

plot3([Pt_moved(1,1) Pt_moved(1,2)], [Pt_moved(2,1) Pt_moved(2,2)], [Pt_moved(3,1) Pt_moved(3,2)]);
plot3([Pt_moved(1,2) Pt_moved(1,3)], [Pt_moved(2,2) Pt_moved(2,3)], [Pt_moved(3,2) Pt_moved(3,3)]);
plot3([Pt_moved(1,3) Pt_moved(1,4)], [Pt_moved(2,3) Pt_moved(2,4)], [Pt_moved(3,3) Pt_moved(3,4)]);
plot3([Pt_moved(1,4) Pt_moved(1,1)], [Pt_moved(2,4) Pt_moved(2,1)], [Pt_moved(3,4) Pt_moved(3,1)]);

scale = 5;
[X,Y,Z] = sphere;

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

tworld_sphere = Rworld_mark * [-20;-20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


tworld_sphere = Rworld_mark * [20;-20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


tworld_sphere = Rworld_mark * [20;20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


tworld_sphere = Rworld_mark * [-20;20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);



Tworld_mark = Tworld_base * Tbase_eff * Teff_inst * TIM; %Tinst_mark_inst; 
Rworld_mark = Tworld_mark(1:3,1:3);
tworld_mark = Tworld_mark(1:3,4); 
Pt_moved = Rworld_mark*Pt + repmat(tworld_mark, 1, 8);

plot3([Pt_moved(1,1) Pt_moved(1,2)], [Pt_moved(2,1) Pt_moved(2,2)], [Pt_moved(3,1) Pt_moved(3,2)]);
plot3([Pt_moved(1,2) Pt_moved(1,3)], [Pt_moved(2,2) Pt_moved(2,3)], [Pt_moved(3,2) Pt_moved(3,3)]);
plot3([Pt_moved(1,3) Pt_moved(1,4)], [Pt_moved(2,3) Pt_moved(2,4)], [Pt_moved(3,3) Pt_moved(3,4)]);
plot3([Pt_moved(1,4) Pt_moved(1,1)], [Pt_moved(2,4) Pt_moved(2,1)], [Pt_moved(3,4) Pt_moved(3,1)]);

scale = 5;
[X,Y,Z] = sphere;

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

tworld_sphere = Rworld_mark * [-20;-20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


tworld_sphere = Rworld_mark * [20;-20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


tworld_sphere = Rworld_mark * [20;20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);


tworld_sphere = Rworld_mark * [-20;20;0] + tworld_mark;

sphere_moved = [X_line;Y_line;Z_line] + repmat(tworld_sphere, 1, size(X,1)*size(X,2));

X_moved = reshape(sphere_moved(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(sphere_moved(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(sphere_moved(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);



