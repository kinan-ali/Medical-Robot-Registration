%Robot Registral tion Lab 
%Initialization: 
clear all; 
clc; 
InitConfig();

%% Question 3: 

%First of all we'll move the camera landoscop a few times to take images of
%the target, and then we can use these images to do 3D localization and 
% calculate the transformation between NAVIGATION SYS and the target. However 
% we will use a lot of points because the movement of the camera is done by
% "the user" so it's not accurate. (but still not too much so the process is tiring)


K_camera = [400, 0, 380; 
    0, 400, 285;
    0, 0, 1]; %the intrinsic parameters of the indoscopic camera

Tmark_cam_cam = inv(Tcam_mark_cam);   %transformation from the camera to the marker camera

n_pts = 4; % number of movements of the camera to take the imgs for 3D Loc

noise = 0;   % noise = 0 turn off noise, 1 turn on noise

OCs = []; %optical centers in the navigation system frame
RDs = []; %ray directions in the navigation system frame
for i=1:n_pts
    targ_pos = GetTargetPosition(noise);
    
    p = GetTargetPosition(noise); % the u,v of the camera in pixels
    p = [p;1];
    p_normalized = inv(K_camera)*p;  %the position of the target in the normalized img frame of the camera 

    measure = GetLocalizerInformation(noise);

    Tnav_mark_cam = measure.mark(1).T; %transformation between the camera marker and the navigation system

    Tnav_cam = Tnav_mark_cam * Tmark_cam_cam; %transformation between the camera and navigation system
    

    %Ray direction of the position of the point 
    ray_direction = p_normalized / norm(p_normalized); % the direction of the target pnt 

    %Converting the ray direction from camera to navigation (rotating the vector)
    R_nav_cam = Tnav_cam((1:3),(1:3));
    current_ray_direction_nav_frame = R_nav_cam*ray_direction;

    current_optical_center_nav_frame = Tnav_cam((1:3), 4);

    OCs = [OCs, current_optical_center_nav_frame];
    RDs = [RDs, current_ray_direction_nav_frame];

    %Moving the camera in a random movement
    MoveCamera(round(rand(3,1)*10), round(rand*5)); %note: round coz the function only takes int values
end

%TRIANGULATION: Least Squares Intersection   
%now that we have the ray directions RDs and the optical centers OCs in the
%navigation system frame, we perform triangulation (Linear System of
%equations to calculate the distances) (NOTE: we are doing the
%triangulation directly in the navigation system frame) 

% We want to find the point 'P_target' closest to all rays defined by OCs and RDs.

A = zeros(3,3);
b = zeros(3,1);

for k = 1:n_pts
    %get origin (O) and direction (u) for the current ray
    u = RDs(:, k);
    O = OCs(:, k);
    
    % Projection matrix onto the plane perpendicular to the ray
    % Formula: (I - u*u')
    proj_matrix = eye(3) - (u * u');
    
    % Accumulate matrices for Least Squares (Ax = b)
    A = A + proj_matrix;
    b = b + (proj_matrix * O);
end

% Solve for the Target Position in Navigation Frame
P_target_nav = A \ b; % Left division solves Ax = b

fprintf('Target Position in Navigation Frame:\n');
disp(P_target_nav);



%Going forward to the robot, doing the hand-eye calibration 

%%

% 1) get the transformation between the nav & needle marker
noise = 0 ;
% position A:
measure = GetLocalizerInformation(noise);
Tnav_mark_inst_A = measure.mark(2).T; %transformation between the camera marker and the navigation system

T_base_eff = GetRobotCurrentPosition();
T_base_inst_A = T_base_eff*Teff_inst;


% position B
current_position = T(1:3, 4);
new_position = current_position + [40,0,0]';
MoveEffPosition(new_position);

measure = GetLocalizerInformation(noise);
Tnav_mark_inst_B = measure.mark(2).T; %transformation between the camera marker and the navigation system

T_base_eff = GetRobotCurrentPosition();
T_base_inst_B = T_base_eff*Teff_inst;


% position C 
current_position = T(1:3, 4);
new_position = current_position + [0,40,0]';
MoveEffPosition(new_position);

measure = GetLocalizerInformation(noise);
Tnav_mark_inst_C = measure.mark(2).T; %transformation between the camera marker and the navigation system

T_base_eff = GetRobotCurrentPosition();
T_base_inst_C = T_base_eff*Teff_inst;


% Constr












