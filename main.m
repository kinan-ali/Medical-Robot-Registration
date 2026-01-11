% Robot Registration Lab 
% Descrioption: Mediacal Robot Registration Lab, HealthTech program master,
% University of Strasbourg 
% Authors: Kinan ALI & Filippo MORINI
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
P_nav_target = A \ b; % Left division solves Ax = b

fprintf('Target Position in Navigation Frame:\n');
disp(P_nav_target);



%Going forward to the robot, doing the hand-eye calibration 


%% Hand-Eye Calibration

% 1. Data Collection
% We need multiple movements to solve AX=XB reliably 
% We will store the absolute poses first
n_robot_moves = 8; 
Tbase_inst_list = {}; % Stores absolute instrument poses in Base frame
Tnav_mark_list = {};  % Stores absolute marker poses in Nav frame


% Move robot to random positions and collect data
for i = 1:n_robot_moves
    % move to random positions, we will move X and Y for safty
    % we want to stay in the safe workspace of the robot
    % rand*100-50 gives an interval of 10 cm around 0
    curr_robot_position = GetRobotCurrentPosition(noise);
    curr_robot_z = curr_robot_position(3,4);
    rand_pos = [round(rand*100-50), round(rand*100-50), curr_robot_z]';
    MoveEffPosition(rand_pos);
    pause(0.1); % for updating the simulation

    % getting the navigation measurment: 
    measure = GetLocalizerInformation(noise);
    Tnav_mark = measure.mark(2).T; 

    % getting robot kinematics: 
    Tbase_eff = GetRobotCurrentPosition(noise);
    % moving from the effector frame to the INSTRUMENT frame: 
    Tbase_inst = Tbase_eff * Teff_inst;

    % store transformations: 
    Tbase_inst_list{end+1} = Tbase_inst; 
    Tnav_mark_list{end+1} = Tnav_mark;

end

% 2. Construct Relative Transformations (Lists A and B)
% A_i * X = X * B_i
% A corresponds to instrument motion, B to marker motion
listeA = {};
listeB = {};
for i = 1:(n_robot_moves-1)

    % A relative motions for the instrument:
    T1 = Tbase_inst_list{i};
    T2 = Tbase_inst_list{i+1};
    A = inv(T1) * T2;
    listeA{end+1} = A;


    % B: Relative motions of the marker (same as before)
    U1 = Tnav_mark_list{i};
    U2 = Tnav_mark_list{i+1};
    B = inv(U1) * U2; 
    listeB{end+1} = B;

end

% solving the hand-eye calibration problem AX = XB : 
Tinst_mark = AXeqXB(listeA, listeB)

% Calculating the static world transform base -> nav :
last_T_base_inst = Tbase_inst_list{end};
last_T_nav_mark  = Tnav_mark_list{end};
Tbase_nav = last_T_base_inst * Tinst_mark * inv(last_T_nav_mark)

%% Following the chaine: 
% now that we have the elements of the chaine, we can follow it: 
% tareget to base frame:
P_nav_target_homo = [P_nav_target; 1];
P_base_target_homo = Tbase_nav * P_nav_target_homo;
P_base_target = P_base_target_homo(1:3);

% we got the target position in the base, HOWEVER this is where the NEEDLE
% should be. we know the offset between the needle and eff: 350 mm: 
% FURTHERMORE: this 350 offest is in the eff frame not in the base! 
% to be able to apply this offset and caculate the desired position of the
% EFFECTOR in the BASE frame, we need to know how the needle is oriented,
% So we need to find where's the trocar! 

% 1. Calculate Trocar Position (The Pivot Point)
% We use the stored positions from your calibration loop (Tbase_inst_list).
% All these needle positions pass through the single fixed trocar point.
% that's because the problem assumes passive wrist and the instrument
% always passes through the trocar:
A_trocar = zeros(3,3);
b_trocar = zeros(3,1);

for i = 1:length(Tbase_inst_list)
    % needle position pose in the base frame: 
    T = Tbase_inst_list{i};
    % The needle axis is the Z-axis of this transformation (3rd column)
    u = T(1:3, 3); 
    O = T(1:3, 4); % The tip position

    % Least Squares Line Intersection (Same math as triangulation)
    proj = eye(3) - u*u';
    A_trocar = A_trocar + proj;
    b_trocar = b_trocar + proj * O;
end

P_base_trocarInst = A_trocar \ b_trocar;

% 2. Calculate the Command using the Trocar
% We know the target P_base_target (where the tip goes).
% We know the needle must pivot through P_trocar.

% A. Find the direction vector (The specific "Tilt" for this target)
insertion_vector = P_base_target - P_base_trocarInst;
needle_axis_z = insertion_vector / norm(insertion_vector);



% B. Apply the 350mm offset ALONG that specific axis
% This is the "subtraction in the effector frame" projected correctly into World space.
target_effector_pos = P_base_target - (350 * needle_axis_z);


% 3. Execute
MoveEffPosition(target_effector_pos);
final_error = ComputeTRE();
fprintf('Final TRE (error of positioning in the target): %.4f mm\n', final_error);


