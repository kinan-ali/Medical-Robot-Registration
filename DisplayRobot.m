function DisplayRobot

global Tworld_base joint_position robot_axes_radius axe1_d axe2_d axe3_d axe4_d robot_axe2_length robot_axe3_length;

H01 = [eye(3) [0; 0; joint_position(1) + axe1_d];
	0 0 0 1];

H12 = [theta2r([0; pi/2;0]) [joint_position(2); axe2_d; 0];
	0 0 0 1];

	
H23 = [theta2r([-pi/2;0;0]) [axe3_d; joint_position(3); 0 ];
	0 0 0 1];	 
	
H34 = [theta2r([0;pi/2; 0;0]) [0; axe4_d; 0];
	0 0 0 1];	 
		
H04 = H01 * H12*H23*H34;	

tbase_eff = H04(1:3,4);


scale = joint_position(1) + robot_axes_radius; %+axe1_d;
[X,Y,Z] = cylinder(robot_axes_radius/scale*[1;1]);

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2));

Rworld_base = Tworld_base(1:3,1:3);
tworld_base = Tworld_base(1:3,4);
body1 = Rworld_base*[X_line;Y_line;Z_line] + repmat(tworld_base, 1, size(X,1)*size(X,2));

X_moved = reshape(body1(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(body1(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(body1(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);

%pause;

scale = robot_axe2_length; %joint_position(2); % +axe2_d;
[X,Y,Z] = cylinder(robot_axes_radius/scale*[1;1]);

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2)) - robot_axes_radius*ones(1,size(X,1)*size(X,2));

Tworld_base2 = Tworld_base*H01*H12;
Rworld_base2 = Tworld_base2(1:3,1:3);
tworld_base2 = Tworld_base2(1:3,4);
body2 = Rworld_base2*[X_line;Y_line;Z_line] + repmat(tworld_base2, 1, size(X,1)*size(X,2));

X_moved = reshape(body2(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(body2(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(body2(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);

%pause;

scale = robot_axe3_length; %joint_position(3); %+axe3_d;
[X,Y,Z] = cylinder(robot_axes_radius/scale*[1;1]);

X_line = scale*reshape(X, 1, size(X,1)*size(X,2));
Y_line = scale*reshape(Y, 1, size(X,1)*size(X,2));
Z_line = scale*reshape(Z, 1, size(X,1)*size(X,2)); % - robot_axe3_length*ones(1,size(X,1)*size(X,2));
;

Tworld_base3 = Tworld_base2*H23;
Rworld_base3 = Tworld_base3(1:3,1:3);
tworld_base3 = Tworld_base3(1:3,4);
body3 = Rworld_base3*[X_line;Y_line;Z_line] + repmat(tworld_base3, 1, size(X,1)*size(X,2));

X_moved = reshape(body3(1,:), size(X,1), size(X,2)); 
Y_moved = reshape(body3(2,:), size(X,1), size(X,2)); 
Z_moved = reshape(body3(3,:), size(X,1), size(X,2)); 

surf(X_moved, Y_moved, Z_moved);

%pause;