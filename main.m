%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ROBOT REGISTRATION %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Display the initial configuration
InitConfig();
% DisplayConfig()

%% 3) REGISTRATION USING THE NAVIGATION SYSTEM

% Known parameters

K_camera = [400 0 380;   % in pixels
            0 400 285;
            0 0 1];

%% Pose estimation of the target in the endoscopic camera

% We will use triangulation using several images to recover the target
% position in the navigation system frame.

% We will move the endoscopic camera, then for each position find the line
% that passes from the normalized image point. After, we will represent
% each line in the navigation system frame through a transformation
% composed by: nav_T_markerCamera * markerCamera_T_camera.

% Then, we will do a triangulation in the navigation system.

n_positions = 4;    % Number of different positions of the endoscopic camera
R_markCam_cam = Rcam_mark_cam';
t_markCam_cam = -R_markCam_cam * tcam_mark_cam;
T_markCam_cam = [R_markCam_cam, t_markCam_cam;  % Transform between marker and 
                   0 0 0 1];                    % endoscopic camera (inverse of Tcam_mark_cam)
noise = 0; % Noise

% We will use the GEOMETRIC SOLVING for triangulation (better for several
% perspective images)
RDs = []; % ray directions in the navigation frame
OCs= []; % optical centers in the navigation system frame
for i = 1:n_positions

    % Transform between navigation system and marker
    measure = GetLocalizerInformation();
    T_nav_markerCam = measure.mark(1).T;
    R_nav_markerCam = T_nav_markerCam(1:3,1:3);
    t_nav_markerCam = T_nav_markerCam(1:3,4);

    % Position of the target in the endoscopic camera image
    targ_pos = GetTargetPosition(noise);

    % Position of the target in the normalized image plane
    K_camera_inv = K_camera \ eye(size(K_camera));
    m = K_camera_inv * [targ_pos; 1];   % m = the vector from camera origin to the target
    m_nav = R_nav_markerCam * R_markCam_cam * m;    % m in the navigation system frame
    ray_direction = m_nav/norm(m_nav);              % Unit vector of m
    
    % Position of the optical center in the navigation system frame
    opt_center_nav = T_nav_markerCam * [t_markCam_cam;1];

    RDs = [RDs, ray_direction];
    OCs = [OCs, opt_center_nav];

    % Inputs of MoveCamera are translation [mm] and angle of rotation [°]    
    MoveCamera(rand(3,1)*3 +1, round(rand*5) + 1)
end

% Now triangulation by using the implemented function
% "geometric_triangulation"
RDs = RDs(1:3,:);
OCs = OCs(1:3,:);
target_nav = geometric_triangulation(RDs, OCs);

%% Calibration of the fixed transform between the instrument and the attached
% rigid body
% NOTE: the calibration is not between the end effector and the rigid body
% because the rigid body is attached to the instrument, therefore the
% transform between the rigid body and the end effector is not fixed.

% Sensors: navigation system + robot

% Position A
measureA = GetLocalizerInformation(noise);
Tnav_markerInst_A = measureA.mark(2).T;    % Transformation from the instrument marker
                                            % to the navigation system
T_base_effA = GetRobotCurrentPosition(noise);
T_base_inst_A = T_base_effA*Teff_inst;

% Position B
current_position = T_base_effA(1:3, 4);
new_position = current_position + [40,0,0]';
MoveEffPosition(new_position);

measureB = GetLocalizerInformation(noise);
Tnav_markerInst_B = measureB.mark(2).T; %transformation between the camera marker and the navigation system

T_base_effB = GetRobotCurrentPosition(noise);
T_base_inst_B = T_base_effB*Teff_inst;

% Position C 
current_position = T_base_effB(1:3, 4);
new_position = current_position + [0,40,0]';    % The three positions need to be NOT aligned
MoveEffPosition(new_position);

measureC = GetLocalizerInformation(noise);
Tnav_markerInst_C = measureC.mark(2).T; %transformation between the camera marker and the navigation system

T_base_effC = GetRobotCurrentPosition();
T_base_inst_C = T_base_effC*Teff_inst;

% Now compute the Akl and Bkl matrices
list_Akl = zeros(4,4,3);    % contains T_instK_instL
list_Bkl = zeros(4,4,3);    % contains T_markerK_instL

list_Akl(:,:,1) = inv(T_base_inst_A)*T_base_inst_B;
list_Akl(:,:,2) = inv(T_base_inst_B)*T_base_inst_C;
list_Akl(:,:,3) = inv(T_base_inst_A)*T_base_inst_C;

list_Bkl(:,:,1) = inv(Tnav_markerInst_A)*Tnav_markerInst_B;
list_Bkl(:,:,2) = inv(Tnav_markerInst_B)*Tnav_markerInst_C;
list_Bkl(:,:,3) = inv(Tnav_markerInst_A)*Tnav_markerInst_C;

% Use the implemented function "eyehand" to solve
T_inst_markerInst = eyehand(list_Akl, list_Bkl);


%% Trocard position

% Use geometric triangulation in the robot base frame
% 1 - Move the robot in order to get three configurations
% 2 - Compute the lines that passes from the needle and the end effector
% 3 - the trocard is at their intersection

% Positions of the end effector in the robot base frame in three
% configurations
OCs = [T_base_effA(1:3,4), T_base_effB(1:3,4),T_base_effC(1:3,4)];

% Same for the instrument
inst_base = [T_base_effA*[teff_inst;1], T_base_effB*[teff_inst;1], T_base_effC*[teff_inst;1]];
inst_base = inst_base(1:3,:);

RDs = inst_base - OCs;    % ray directions vectors
for i=1:size(RDs,2)
    unitvec = RDs(1:3,i) / norm(RDs(1:3,i));
    RDs(1:3,i) = unitvec;
end

trocard_base = geometric_triangulation(RDs, OCs);


%% Computation of trajectories in the robot base frame
% We have E, T in the robot base

InitConfig() % return to initial configuration

measure = GetLocalizerInformation(noise);
T_nav_markerInst = measure.mark(2).T;
T_base_eff = GetRobotCurrentPosition(noise);
% Target in the base frame
target_base = T_base_eff * Teff_inst * T_inst_markerInst * inv(T_nav_markerInst) * [target_nav;1];
target_base = target_base(1:3);

% Computation of trajectory
trajectory = (target_base - trocard_base) / norm(target_base - trocard_base);
eff_position = target_base - trajectory*350;

% Move the effector in the wanted position
MoveEffPosition(eff_position)
DisplayConfig()
ComputeTRE()

%% Transform chain
% - from target in endoscopic camera to navigation system through rigid
% body
% - from nav. system to rigid body of the instrument then to instrument, then
% to end effector and robot base
% --> this is for the target in the robot base frame, the instrument is
% more simple (instrument->end effector->base)

% NOTE: The frame of the robot effector Feff is centered on the robot 
% effector and oriented the same way as Finst, so that eff_R_inst=I3
% and eff_t_inst=(0;0;350)' mm. This means that we don't need the position
% of the trocard!

% NOTE: see the text, there are some unknown elements. Explain why one of
% those needs to be found, while the other are not necessary! E. g. the
% position of the trocard, or the position of navigation frame (slide 294:
% absolute information are of no use). Why, instead, we need other frames
% positions?


%% 4) REGISTRATION WITHOUT THE NAVIGATION SYSTEM

% ADVANTAGE: do not use the navigation system, no need of line of view, no 
% need of markers on the instrument nor the endoscopic camera (therefore no
% need to calibrate the fixed transform between them! --> less transform
% involved usually lead to a smaller error)
% DISADVANTAGE: need to have the instrument tip and target in the same
% image (?), marker on the instrument tip, might be large, how precise is
% the position of the target detected (I am using the images to first
% obtain the pose of the instrument, then use that pose to have the
% position of the target: I am using twice the images which are the main
% source of errors!)

% + see the effect of noise: how robust are the methods?

% USED MARKER? Propose a structure (+ photo of it + explain why is suited
% for that)

% Data: cam_T_instr
% FOR instrument: from instrument to end effector and to robot base (as before)
% FOR target:
% - triangulation to obtain the target in the endoscopic camera frame (the
% relative pose between images will be computed referring to the instrument
% which will be in a fixed position and will "play the role of the
% navigation system")
% - from camera to instrument to end effector to robot base

InitConfig() % return to initial configuration
MoveCameraInitialPosition()

% First: find an image in which both the target and the marker are visible
MoveCamera([10 5 -15]', 0);

% Use of three images for triangulation (two movements)
targ_image1 = GetTargetPosition(noise);
Tcam_inst1 = GetInstrumentPosition(noise);

% First movement   
MoveCamera(rand(3,1) +1, round(rand*3) + 1)   % Inputs of MoveCamera are translation [mm]
                                                % and angle of rotation [°]
targ_image2 = GetTargetPosition(noise);
Tcam_inst2 = GetInstrumentPosition(noise);

% Second movement
MoveCameraInitialPosition()
MoveCamera([5 5 -20]', 0);
MoveCamera(rand(3,1) +1, round(rand*3) + 1)
targ_image3 = GetTargetPosition(noise);
Tcam_inst3 = GetInstrumentPosition(noise);

targ_images = [targ_image1, targ_image2, targ_image3];
% Normalized image plane points in the correspondent image plane
m = K_camera_inv * [targ_images; ones(1, size(targ_images,2))];
m = m(1:3,:);

% Relative pose between optical centers
Tcam1_cam2 = Tcam_inst1 / Tcam_inst2;
Tcam1_cam3 = Tcam_inst1 / Tcam_inst3;

% Choose first optical center as reference frame
OC1 = [0 0 0]';
OC2 = Tcam1_cam2(1:3,4);
OC3 = Tcam1_cam3(1:3,4);
OCs = [OC1, OC2, OC3];

RD1 = m(:,1)/norm(m(1:3,1));
dir2 = Tcam1_cam2(1:3,1:3)*m(:,2);
RD2 = dir2 / norm(dir2);
dir3 = Tcam1_cam3(1:3,1:3)*m(:,3);
RD3 = dir3 / norm(dir3);
RDs = [RD1, RD2, RD3];

target_camera1 = geometric_triangulation(RDs, OCs);

T_base_eff = GetRobotCurrentPosition(noise);

% Target in the robot base frame
target_base = T_base_eff * Teff_inst * inv(Tcam_inst1) * [target_camera1;1];
target_base = target_base(1:3);


% FROM NOW THE SOLVING METHODS ARE THE SAME AS BEFORE

% Trocart in camera frame: is the variable 'trocard_base'

% Computation of trajectory
trajectory = (target_base - trocard_base) / norm(target_base - trocard_base);
eff_position = target_base - trajectory*350;

% Move the effector in the wanted position
MoveCameraInitialPosition()
MoveEffPosition(eff_position)
DisplayConfig()
ComputeTRE()

