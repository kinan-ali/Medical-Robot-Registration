%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ROBOT REGISTRATION %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: Filippo Morini & Kinan Ali

clear
clc
close all


%% INITIAL CONDITIONS

% Display the initial configuration
InitConfig();
%DisplayConfig()

% Set noise
noise_robot = 0;
noise_nav = 0;
noise_cam = 0;
noise_inst = 0;


%% 3) REGISTRATION USING THE NAVIGATION SYSTEM

% Known intrinsic parameters
K_camera = [400 0 380;   % in pixels
            0 400 285;
            0 0 1];
K_camera_inv = K_camera \ eye(size(K_camera));

%%%%% POSITION OF THE TARGET IN THE NAVIGATION SYSTEM FRAME 
% The estimation of the position of the target in the navigation system
% frame is carried out by using triangulation from several perspective images.
% Triangulation is performed in the navigation system frame.
% Implemented function: geometric_triangulation

% 1) Ray directions and optical centers for triangulation
% -------------------------------------------------------------------------
n_positions = 4;    % Number of different positions of the endoscopic camera

% Transform between marker and endoscopic camera (inverse of Tcam_mark_cam)
R_markCam_cam = Rcam_mark_cam';
t_markCam_cam = -R_markCam_cam * tcam_mark_cam;
T_markCam_cam = [R_markCam_cam, t_markCam_cam;  
                   0 0 0 1];                   
% Determining RDs and OCs to then perform triangulation
RDs = [];   % ray directions in the navigation frame
OCs= [];    % optical centers in the navigation system frame
move = [5 0 0 -8;
         -2 5 0 6;
         0 -6 3 3;
         1 1 1 1];
for i = 1:n_positions
    
    % Transform between navigation system and marker
    measure = GetLocalizerInformation(noise_nav);
    T_nav_markerCam = measure.mark(1).T;
    R_nav_markerCam = T_nav_markerCam(1:3,1:3);
    t_nav_markerCam = T_nav_markerCam(1:3,4);

    % Position of the target in the endoscopic camera image
    targ_pos = GetTargetPosition(noise_cam);

    % Ray directions in the navigation system frame
    m = K_camera_inv * [targ_pos; 1];               % normalized image plane vector
    m_nav = R_nav_markerCam * R_markCam_cam * m;    % m in the navigation system frame
    ray_direction = m_nav/norm(m_nav);              % unit vector of m
    
    % Position of the optical center in the navigation system frame
    opt_center_nav = T_nav_markerCam * [t_markCam_cam;1];
    
    RDs = [RDs, ray_direction];
    OCs = [OCs, opt_center_nav];

    % Inputs of MoveCamera are translation [mm] and angle of rotation [°]
    MoveCameraInitialPosition()
    MoveCamera([move(1,i),move(2,i), move(3,i)]', move(4,i))
end

% 2) Triangulation in navigation systema frame by using the 
% implemented function "geometric_triangulation"
% -------------------------------------------------------------------------
RDs = RDs(1:3,:);   % ray directions
OCs = OCs(1:3,:);   % optical centers
target_nav = geometric_triangulation(RDs, OCs);


%%%%% CALIBRATION OF THE FIXED TRANSFORM BETWEEN 
%%%%% THE INSTRUMENT AND THE ATTACHED RIGID BODY 
% Sensors: navigation system (for the marker) + robot (for the instrument)
% Implemented function: eyehand, loc_3Dhorn
% A(kl) * X = X * B(kl) --> Need at least 2 movements (3 configurations: A, B, C) to solve for X
%                           Here 4 configurations were used for robustness

% 1) Collect the transforms for the A(kl) and B(kl) list of matrices
% -------------------------------------------------------------------------
n_config = 4;
% Take the reference configuration as starting point for robustness to noise
Tbaseeff = GetRobotCurrentPosition(noise_robot);
current_position = Tbaseeff(1:3,4);
% Movements (of end effector)
move = [30 20 5;
        -20 30 -10;
        -10 -40 15;
        0 0 0]';                   % last movement won't be measured
% Initialize the transformations
T_nav_markerInst_list = zeros(4,4,n_config);
T_base_inst_list = zeros(4,4,n_config);
T_base_eff_list = zeros(4,4,n_config);

for n = 1:n_config
    measure = GetLocalizerInformation(noise_nav);
    T_nav_markerInst_list(:,:,n) = measure.mark(2).T;              % Pose of the marker in the navigation system frame
    T_base_eff_list(:,:,n) = GetRobotCurrentPosition(noise_robot);
    T_base_inst_list(:,:,n) = T_base_eff_list(:,:,n)*Teff_inst;    % Pose of the instrument in the robot base frame
    new_position = current_position + move(:,n);                   % Movement to be performed 
    MoveEffPosition(new_position);
end

% 2) Computation the Akl and Bkl list of matrices
% -------------------------------------------------------------------------
list_Akl = zeros(4,4,6);    % contains T_instK_instL
list_Bkl = zeros(4,4,6);    % contains T_markerK_markerL

list_Akl(:,:,1) = inv(T_base_inst_list(:,:,1))*T_base_inst_list(:,:,2);
list_Akl(:,:,2) = inv(T_base_inst_list(:,:,1))*T_base_inst_list(:,:,3);
list_Akl(:,:,3) = inv(T_base_inst_list(:,:,1))*T_base_inst_list(:,:,4);
list_Akl(:,:,4) = inv(T_base_inst_list(:,:,2))*T_base_inst_list(:,:,3);
list_Akl(:,:,5) = inv(T_base_inst_list(:,:,2))*T_base_inst_list(:,:,4);
list_Akl(:,:,6) = inv(T_base_inst_list(:,:,3))*T_base_inst_list(:,:,4);

list_Bkl(:,:,1) = inv(T_nav_markerInst_list (:,:,1))*T_nav_markerInst_list (:,:,2);
list_Bkl(:,:,2) = inv(T_nav_markerInst_list (:,:,1))*T_nav_markerInst_list (:,:,3);
list_Bkl(:,:,3) = inv(T_nav_markerInst_list (:,:,1))*T_nav_markerInst_list (:,:,4);
list_Bkl(:,:,4) = inv(T_nav_markerInst_list (:,:,2))*T_nav_markerInst_list (:,:,3);
list_Bkl(:,:,5) = inv(T_nav_markerInst_list (:,:,2))*T_nav_markerInst_list (:,:,4);
list_Bkl(:,:,6) = inv(T_nav_markerInst_list (:,:,3))*T_nav_markerInst_list (:,:,4);

% 3) Use the implemented function "eyehand" (contains "loc_3Dhorn") to solve
% -------------------------------------------------------------------------
T_inst_markerInst = eyehand(list_Akl, list_Bkl);


%%%% TROCARD POSITION
% Use geometric triangulation in the robot base frame
% 1 - Move the robot in order to get many configurations
% 2 - Compute the lines that passes from the needle and the end effector
% 3 - the trocard is at their intersection (use of the implemented function "geometric_calibration")

% 1) Use the configurations of eye-hand calibration
% -------------------------------------------------------------------------
% Positions of the end effector in the robot base frame in four configurations
OCs = [T_base_eff_list(1:3,4,1), T_base_eff_list(1:3,4,2), T_base_eff_list(1:3,4,3), T_base_eff_list(1:3,4,4)];
% Same for the instrument
inst_base = [T_base_eff_list(:,:,1)*[teff_inst;1], T_base_eff_list(:,:,2)*[teff_inst;1], ...
                    T_base_eff_list(:,:,3)*[teff_inst;1], T_base_eff_list(:,:,4)*[teff_inst;1]];
inst_base = inst_base(1:3,:);

% 2) Ray directions and optical centers in the robot base frame
% -------------------------------------------------------------------------
RDs = inst_base - OCs;                          % ray directions vectors
for i=1:size(RDs,2)
    unitvec = RDs(1:3,i) / norm(RDs(1:3,i));
    RDs(1:3,i) = unitvec;                       % ray directions unit vectors
end

% 3) Find the intersection through triangulation
% -------------------------------------------------------------------------
trocard_base = geometric_triangulation(RDs, OCs);

%%%%% COMPUTATION OF TRAJECTORIES IN THE ROBOT BASE FRAME 
% The aim is to bring the instrument in the target position
% 1 - We have the entry point E (trocard) in the robot base frame,
% therefore we need to compute the target T in the robot base frame
% 2 - Compute the trajectory (line that passes from E and T) and then compute
% the position in which to bring the end effector
% 3 - Move the end effector

InitConfig() % return to initial configuration

% 1) Compute the target in the robot base frame
% -------------------------------------------------------------------------
% Computation of target point in robot base frame
measure = GetLocalizerInformation(noise_nav);
T_nav_markerInst = measure.mark(2).T;           
T_base_eff = GetRobotCurrentPosition(noise_robot);
% Use a tranform chain to go from navigation system to base frame
target_base = T_base_eff * Teff_inst * T_inst_markerInst * inv(T_nav_markerInst) * [target_nav;1];
target_base = target_base(1:3);                 % target in the base frame

% 2) Computation of trajectory
% -------------------------------------------------------------------------
trajectory = (target_base - trocard_base) / norm(target_base - trocard_base);
eff_position = target_base - trajectory*350;    % desired end effector position

% 3) Move the effector in the wanted position
% -------------------------------------------------------------------------
MoveEffPosition(eff_position)
DisplayConfig()
ComputeTRE()


%% 4) REGISTRATION WITHOUT THE NAVIGATION SYSTEM
% To obtain the position of the target in the endoscopic camera frame:
% 1 - Move the camera in eight positions from which it can see both the
% target and the instrument marker
% 2 - Compute the relative transforms between the camera positions by using
% the pose of the marker obtained by GetInstrumentPosition()
% 3 - Geometric triangulation (using the implemented function "geometric_triangulation") 

InitConfig()                    % return to initial configuration

% 1) Move the camera and compute the optical centers and ray directions
% -------------------------------------------------------------------------
% Movements: each of these movements provide an image in which both the target
% and the marker are visible
nmove = 8;
move_t = [10 5 -15;                 % translations     
          14 6 -8;          
          8 3 -20;
          12 7 -15];
move_angle = [0,-40, 30, 50];       % angles
% Initialization
targ_images = zeros(2,nmove); 
Tcam_inst = zeros(4,4,nmove);
T_cam1_camN = zeros(4,4,nmove-1);
OCs = zeros(3,nmove);
RDs = zeros(3,nmove);

for ii = 1:nmove
    MoveCameraInitialPosition()                             % To start from a repeatable position
    if mod(ii,2) ~= 0
        MoveCamera(move_t(round(ii/2),:)', move_angle(round(ii/2)));  
    else
        MoveCamera((move_t(round(ii/2),:)' + rand(3,1)), (move_angle(round(ii/2)) + rand))
    end
    targ_images(:,ii) = GetTargetPosition(noise_cam);       % target positions in the images
    Tcam_inst(:,:,ii) = GetInstrumentPosition(noise_inst);  % Pose of the instrument marker in camera
    % Relative pose of cameras (reference: camera in position 1)
    if ii > 1
        T_cam1_camN(:,:,ii-1) = Tcam_inst(:,:,1) * inv(Tcam_inst(:,:,ii));    
    end
    % Optical centers in camera position 1 frame
    if ii < 2
        OCs(:,ii) = [0 0 0]';
    else
        OCs(:,ii) = T_cam1_camN(1:3,4,ii-1);
    end
end

% Normalized image plane points in the correspondent image plane
m = K_camera_inv * [targ_images; ones(1, size(targ_images,2))];
m = m(1:3,:);
% Ray directions in camera position 1 frame
RDs(:,1) = m(:,1)/norm(m(:,1));
for jj = 1:(nmove-1)
    dir = T_cam1_camN(1:3,1:3,jj)*m(:,jj+1);
    RDs(:,jj+1) = dir/ norm(dir);
end

% 2) Geometric triangulation in the frame of the device when we took image1
% -------------------------------------------------------------------------
target_camera1 = geometric_triangulation(RDs, OCs);

T_base_eff = GetRobotCurrentPosition(noise_robot);

% 3) Target in the robot base frame
% -------------------------------------------------------------------------
Tcam_inst1 = Tcam_inst(:,:,1);
target_base = T_base_eff * Teff_inst * inv(Tcam_inst1) * [target_camera1;1];
target_base = target_base(1:3);

% NOTE: From here the solving methods are the same as before in exercise 3
% Trocart in camera frame: is the variable 'trocard_base'

% 4) Computation of trajectory
% -------------------------------------------------------------------------
trajectory = (target_base - trocard_base) / norm(target_base - trocard_base);
eff_position = target_base - trajectory*350;    % desired end effector position

% 5) Move the effector in the wanted position
% -------------------------------------------------------------------------
MoveCameraInitialPosition()
MoveEffPosition(eff_position)
DisplayConfig()
ComputeTRE()



%% 5) ASSESSMENT OF ERRORS PROPAGATION
% ---------------------------------------------------------
% PART 1: Setup and Analytical Solution (Variance Propagation)
% ---------------------------------------------------------
clc;
InitConfig();
fprintf('--- Q5.1: Analytical Variance Propagation ---\n');

% 1. Define Geometry (Simulation Configuration)
% We explicitly define a realistic configuration in the Robot Base Frame 
% for this analysis.
% ---------------------------------------------------------
% L: Length of the instrument 
L_inst = 350; %mm

noise = 0;
T_base_eff = GetRobotCurrentPosition(noise);

P_base_eff = T_base_eff(1:3,4);
P_base_trocar = [-350; -100; 800]; % computed earlier in Q3 (triangulation)

% Calculate the Mean Tip Position (Y) based on this geometry
% Vector V goes from Effector to Trocar
V = P_base_trocar - P_base_eff;
dist_V = norm(V);
u_hat = V / dist_V; % Unit direction vector of the needle

% The Tip is located distance L along this direction starting from Effector
P_base_tip = P_base_eff + L_inst * u_hat;

% 2. Define Input Uncertainty (Sigma_X)
% ---------------------------------------------------------
% Instructions: "Standard deviation ... isotropic and of amplitude 1 mm"
sigma_trocar = 1.0; % 1 mm standard deviation
sigma_trocar_x = 1/3; 
Sigma_X = eye(3) * sigma_trocar_x; % Covariance matrix of the Trocar

% 3. Calculate Jacobian (J)
% ---------------------------------------------------------
% Formula derived: J = (L / ||V||) * (I - u*u')
% This maps small changes in Trocar (X) to changes in Tip (Y)
I3 = eye(3);
J = (L_inst / dist_V) * (I3 - (u_hat * u_hat'));

% 4. Analytical Error Propagation
% ---------------------------------------------------------
Sigma_Y_base = J * Sigma_X * J';

fprintf('Analytical Covariance Matrix (Base Frame) calculated.\n');
disp(Sigma_Y_base);

% 5. Prepare for Plotting (Convert to World Frame)
% ---------------------------------------------------------
% The Rotation matrix from Base to World.
R_world_base = Tworld_base(1:3, 1:3);
t_world_base = Tworld_base(1:3, 4);

% Rotate the Covariance Matrix to World Frame
% Sigma_world = R * Sigma_base * R'
Sigma_Y_world_analytic = R_world_base * Sigma_Y_base * R_world_base';

% Transform the Center Point (Tip) to World Frame for plotting later
P_world_tip_analytic = Tworld_base * [P_base_tip; 1];
P_world_tip_analytic = P_world_tip_analytic(1:3);

fprintf('Converted Analytical results to World Frame for plotting.\n');

% ---------------------------------------------------------
% PART 2: Numerical Evaluation (Monte Carlo Simulation)
% ---------------------------------------------------------
fprintf('\n--- Q5.2: Numerical Simulation (Monte Carlo) ---\n');

% 1. Simulation Parameters
% ---------------------------------------------------------
n_samples = 1000; % Number of iterations 
P_base_tips_noisy = zeros(3, n_samples); % Storage for results

fprintf('Simulating %d noisy samples...\n', n_samples);

% 2. Monte Carlo Loop
% ---------------------------------------------------------
for i = 1:n_samples
    % A. Generate Noise
    % isotropic noise = 1mm -> sigma_x^2 = 1/3
    % randn() generates Normal distribution with mean=0, std=1.
    noise_vector = randn(3, 1) * sqrt(sigma_trocar_x); 
    
    % B. Perturb the Input (Trocar)
    % We assume ONLY the trocar has error 
    P_base_trocar_noisy = P_base_trocar + noise_vector;
    
    % C. Recalculate Output (Tip Position) using the EXACT Geometric Model
    % Y = P_eff + L * (u / ||u||)
    
    % Recalculate vector V_noisy from Effector -> Noisy Trocar
    V_noisy = P_base_trocar_noisy - P_base_eff;
    
    % Normalize
    u_hat_noisy = V_noisy / norm(V_noisy);
    
    % Calculate new Tip Position
    P_tip_noisy = P_base_eff + L_inst * u_hat_noisy;
    
    % D. Store result
    P_base_tips_noisy(:, i) = P_tip_noisy;
end

% 3. Calculate Numerical Covariance
% ---------------------------------------------------------
% Use Matlab's built-in cov() function. 
% Note: cov expects rows = samples, so we transpose (') the matrix.
Sigma_Y_base_numeric = cov(P_base_tips_noisy');

fprintf('Numerical Covariance Matrix (Base Frame) calculated.\n');
disp(Sigma_Y_base_numeric);

% 4. Convert to World Frame for Comparison
% ---------------------------------------------------------
% We must rotate the numerical covariance to World Frame just like we did
% for the analytical one.
Sigma_Y_world_numeric = R_world_base * Sigma_Y_base_numeric * R_world_base';

% We also transform the cloud of points to World Frame for plotting
P_world_tips_noisy = zeros(3, n_samples);
for i = 1:n_samples
    pt_base_homo = [P_base_tips_noisy(:, i); 1];
    pt_world_homo = Tworld_base * pt_base_homo;
    P_world_tips_noisy(:, i) = pt_world_homo(1:3);
end

fprintf('Converted Numerical results to World Frame for plotting.\n');

% ---------------------------------------------------------
% PART 3: Visualization and Comparison
% ---------------------------------------------------------
fprintf('\n--- Q5.3: Visualization & Comparison ---\n');

% 1. Calculate Eigenvalues and Eigenvectors
% ---------------------------------------------------------
% Analytical:
[V_ana, D_ana] = eig(Sigma_Y_world_analytic);

% Use abs() to handle tiny negative noise (e.g., -1e-16) preventing complex
% numbers in sqrt()
radii_ana = sqrt(abs(diag(D_ana))); 

% Numerical:
[V_num, D_num] = eig(Sigma_Y_world_numeric);

% Same here for numerical stability
radii_num = sqrt(abs(diag(D_num)));

% Print Comparison to Console
fprintf('\nComparison of Principal Axes (Standard Deviations):\n');
fprintf('%-15s %-15s %-15s\n', 'Axis', 'Analytical', 'Numerical');
for k = 1:3
    fprintf('Axis %d:         %.4f          %.4f\n', k, radii_ana(k), radii_num(k));
end
fprintf('(Note: One axis should be close to 0, representing the needle axis direction)\n');

% 2. Plotting
% ---------------------------------------------------------
figure('Name', 'Q5: Error Propagation Comparison'); clf; hold on; grid on; axis equal;
xlabel('X World (mm)'); ylabel('Y World (mm)'); zlabel('Z World (mm)');
title('Comparison of Error Ellipsoids (World Frame)');
view(3); % Set 3D view

% A. Plot the "Cloud" of Monte Carlo points (Grey dots)
plot3(P_world_tips_noisy(1,:), P_world_tips_noisy(2,:), P_world_tips_noisy(3,:), ...
      '.', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Monte Carlo Cloud');

% B. Plot the Mean Tip Position (Black Marker)
plot3(P_world_tip_analytic(1), P_world_tip_analytic(2), P_world_tip_analytic(3), ...
      'k+', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Mean Tip Position');

% C. Plot Analytical Ellipsoid
% Note: plot_ellipsoid(center, rot_mat, axe_lengthes, scale)
% We use scale=1 (1 Standard Deviation)
% We rely on the function's default coloring. 
fprintf('\nPlotting Analytical Ellipsoid...\n');
plot_ellipsoid(P_world_tip_analytic, V_ana, radii_ana, 1);

% D. Plot Numerical Ellipsoid
fprintf('Plotting Numerical Ellipsoid...\n');
% We assume the numerical one will overlay the analytical one.
plot_ellipsoid(P_world_tip_analytic, V_num, radii_num, 1);

legend('Location', 'best');
hold off;



%% 6) POSITIONING UNDER THE FEEDBACK OF THE CAMERA
% Before the robot was moving in an open loop
% --> now closed-loop based on position based visual servoing
% REFERENCE: target position in sensor frame
% FEEDBACK MEASURE: position of the tip in the sensor frame
% --> try to minimize the ERROR between reference and feedback
% Jacobian which links the robot motion to the variation of the error
% Variation of error: dot(e) = - dot(Ptip) --> in the image device frame!
% Robot motion: J(q)*dot(q) = Vrobot (velocity of the end effector: it is a
% operational space control? YES, operational space control --> see * (end of page)
% Jacobian is the Rotation R_cam_base * (R that inverts the x and y
% rotations: due to the trocard inversion)  NOTE: leverage effect is
% ignored

% 1) Find the target in a chosen camera position (will be used as reference)
% -------------------------------------------------------------------------
% Move in a position in which is possible to see both target and instrument
% tip
InitConfig()
MoveCamera([10 5 -15]', 0); % this position is found in exercise 4

% Target in current camera position (through relative pose)
Tcam_inst = GetInstrumentPosition(noise_inst);
T_cam_cam1 = Tcam_inst*inv(Tcam_inst1);         % from camera1 (see ex4) to current camera position
target_camera = T_cam_cam1*[target_camera1;1]; 
target_camera = target_camera(1:3);             % target in the current camera position

% 2) Point Based Visual Servoing
% -------------------------------------------------------------------------
dt = 0.1;
max_iter = 300;
lambda = 0.5;               % is the gain of the controller
threshold_error = 0.1;      % max error allowed
L = teff_inst(3);           % distance from end effector to instrument tip

% Initial conditions
inst_camera = Tcam_inst(1:3,4);         % initial position of the instrument tip in camera frame

error = target_camera - inst_camera;    % initial error in camera frame
error_plot = norm(error);               % error plot will save errors through iterations

max_vel = 50;                           % [mm/s] is the max effector velocity (for safety) 
iter = 0;

while (norm(error) > threshold_error) && (iter < max_iter)
    % Quantities used to update the Jacobian
    Tbase_eff = GetRobotCurrentPosition(noise_robot);
    eff_base = Tbase_eff(1:3,4);                % effector in the robot frame
    d_vect = trocard_base - eff_base;           % vector from effector to trocard
    d = norm(d_vect);
    u = d_vect / d;                             % unit vector of d
    J = ((1 - (L/d))*eye(3) + ((L/d)*(u*u')));      % Jacobian (to be used in the robot frame)
    % Transfer error in the robot base frame
    Tbase_cam = Tbase_eff*Teff_inst*inv(Tcam_inst); % transform from camera to robot base
    error = Tbase_cam(1:3,1:3)*error;               % error in robot frame
    % Control law
    eff_vel = inv(J)*lambda*error;
    eff_vel = Tbase_eff(1:3,1:3)'*eff_vel;          % effector velocity in effector frame
    if norm(eff_vel) > max_vel                      % control to not have a too high velocity 
        eff_vel = (eff_vel / norm(eff_vel)) * max_vel;
    end
    MoveEffVelocity(eff_vel,dt)                 % use the robot manipulator (Inverse Differential Kinematics)
    % Update the error with new instrument tip position
    Tcam_inst = GetInstrumentPosition(noise_inst);   
    inst_camera = Tcam_inst(1:3,4);             % instrument tip position in camera frame
    error = target_camera - inst_camera;        % update the error
    error_plot = [error_plot, norm(error)];
    iter = iter +1;
end

% 3) Display results
% -------------------------------------------------------------------------
DisplayConfig()
final_error = ComputeTRE()
f1 = figure;
plot(error_plot, 'LineWidth',2)
title('error trend for PBVS')
xlabel = ('Iterations');
ylabel = ('Error [mm]');
