%initVS_eye_in_hand

%global object configuration lim_average_im_err Te nb_samp_max display_mode im_pts_ref_ini im_pts_ref_fin lambda im_pts_ini type_jacobian Timg_obj_ref_ini Timg_obj_ref_fin Timg_obj_ini imager_real imager_est robot Teff_tool Teff_tool_real Tbase_targ type_robot_control noise_ampl imaging_device VisualServoing pose sing_val_handling traj_img_obj_ref traj_im_pts_ref traj_liste_visible_points nb_info ind_fin_traj ind_fin_traj_targ save_mode video_mode;

global Tworld_cam_ini Tworld_cam Tworld_org Tworld_loc Tworld_base Tbase_eff Teff_inst tworld_troc_cam tworld_troc_inst tbase_troc_inst Tcam_mark_cam torg_target camera_real endoscope_length endoscope_radius axe1_d axe2_d axe3_d axe4_d robot_axes_radius joint_position instrument_radius instrument_length patient_radius Tworld_patient im_pts robot_axe2_length robot_axe3_length mark_inst im_pts_mark_inst TIM %Tinst_mark_inst ;

fprintf(1,'Initialization\n');

load_patient_body;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres de dimensions
endoscope_length = 350; %500;
endoscope_radius = 5;
robot_axes_radius = 50;
robot_axe2_length = 400;
robot_axe3_length = 400;


instrument_length = 400; %700;
instrument_radius = 2.5;

axe1_d = 0;
axe2_d = robot_axes_radius; %100; %500;
axe3_d = robot_axes_radius; %100;
axe4_d = 0;


% Marqueurs sur l'instrument
mark_inst = [[0;0;0], [0;0;-10]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position des differents reperes

% Reperes fixes
thetax = 0;
thetay = 0;
thetaz = 0;
Rworld_base = theta2r([thetax,0,0])*theta2r([0,thetay,0])*theta2r([0,0,thetaz]);
tworld_base = [500; 200; -500]; %[500; 500; -500];
Tworld_base = [Rworld_base tworld_base; 0 0 0 1];

thetax = 0;
thetay = 0;
thetaz = 0;
Rworld_org = theta2r([thetax,0,0])*theta2r([0,thetay,0])*theta2r([0,0,thetaz]);
tworld_org = [0; 0; 100];
Tworld_org = [Rworld_org tworld_org; 0 0 0 1];


patient_radius = 282.84;
Tworld_patient = [theta2r([-pi/2;0;0])  [0;0;100]; 0 0 0 1];

% Les trocarts sont a� la surface du patient
tworld_troc_cam = [100; -100; 300]; %[200; 0; 300];

tworld_troc_inst = [150; 100; 300]; %[200; 200; 300];

Tbase_world = inv(Tworld_base);
tbase_troc_inst = Tbase_world(1:3,4) + Tbase_world(1:3,1:3)*tworld_troc_inst;


thetax = -pi/2 - pi/4;
thetay = 0;
thetaz = pi;
Rworld_loc = theta2r([0,0,thetaz])*theta2r([thetax,0,0])*theta2r([0,thetay,0]);
tworld_loc = [0; 1000; 1500];
Tworld_loc = [Rworld_loc tworld_loc; 0 0 0 1];
%pause;

thetax = pi/2;
thetay = 0;
thetaz = 0;
Rcam_mark_cam = theta2r([thetax,0,0])*theta2r([0,thetay,0])*theta2r([0,0,thetaz]);
tcam_mark_cam = [0; -20; -350]; %[0; -20; -500]
Tcam_mark_cam = [Rcam_mark_cam tcam_mark_cam; 0 0 0 1];


torg_target = [0;0;0];

thetax = pi/2;
thetay = 0;
thetaz = 0;
Rinst_mark_inst = theta2r([thetax,0,0])*theta2r([0,thetay,0])*theta2r([0,0,thetaz]);
tinst_mark_inst = [0; -20; -400]; %[0; -20; -700];
TIM = [Rinst_mark_inst tinst_mark_inst; 0 0 0 1];
%Tinst_mark_inst

Reff_inst = eye(3);
teff_inst = [0; 0; 350]; %[0; 0; 500];
Teff_inst = [Reff_inst teff_inst; 0 0 0 1];


% Reperes mobiles
thetax = -pi/2 - pi/3;
thetay = pi/6;
thetaz = pi/2;
Rworld_cam_ini = theta2r([0,0,thetaz])*theta2r([thetax,0,0])*theta2r([0,thetay,0]);
tworld_cam_ini = tworld_troc_cam + Rworld_cam_ini*[0;0;40];
Tworld_cam_ini = [Rworld_cam_ini tworld_cam_ini; 0 0 0 1];

Tworld_cam = Tworld_cam_ini;


q_ini = [1100; -200; -100];
joint_position = q_ini;
Tbase_eff_ini = MGD(q_ini, tbase_troc_inst);
Tbase_eff = Tbase_eff_ini;
%pause;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paramètres caméra
camera_real.K = [400 0 380; 0 400 285; 0 0 1];
camera_real.corner_im = [1 760 1 570]; 
%imager_real.corner_im = [1 1600 1 1200]; 

fprintf(1,'Initialization done\n');

%



%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition des parametres graphiques
%%display_mode = 'on '; %'off'; %'on '
%display_mode = 'off';
%
%%save_mode = 'off';
%save_mode = 'on '; %'off';
%
%video_mode = 'on '; %'off'; % pour affichage de traj 3D en mode video (image par image)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition de la scene
%
%% Definition du type d'imageur
%imaging_device = 'slice_proj';
%%imaging_device = 'persp_proj';
%
%% Definition du type d'asservissement visuel
%%VisualServoing = '3D';
%VisualServoing = '2D';
%
%
%% Definition de l'objet
%
%object = [[0;0;0], [0;0;100], [0;0;0], [40; 40; 100], [40; 0;0], [40;0;100], [40;0;0], [40;40;100], [40;40;0], [40;40;100], [0;40;0], [0; 40; 100], [40;40;0], [40; 0; 100]];
%
%%object = [[100; -100; 0], [100; 100; 0], [-100; 100; 0], [-100; -100; 0], [0; 0; 50]];
%%object = [[100; -100; 0], [100; 100; 0], [-100; 100; 0], [-100; -100; 0]];
%%object = [[50; 0; 0], [50; 50; 0], [0; 50; 0], [0; 0; 0]];
%%object = [[50; 0; 0], [50; 50; 0], [0; 50; 0]];
%%object = [[50; 0; 0], [50; 50; 0], [0; 50; 0], [0; 0; 0], [30; 30; 30]];
%%object = [[0;0;0] [25; 0; 10], [25; 25; 10], [-25; 25; -10], [-25; -25; -10]];
%%object = [[50; 0; 0], [50; 50; 0], [0; 50; 0], [0; 0; 0], [30; 30; 30], [10; 10; -30]];
%%object = [[10; 0; 0], [10; 10; 0], [0; 10; 0]];
%
%% Definition de la configuration
%%configuration = 'eyeinhand'; % 'eyetohand'
%configuration = 'eyetohand'
%
%
%
%% Definition du mouvement de l'objet
%movement = 'n';
%%movement = 'y';
%
%% Definition de la consigne
%reference_type = 'position  ';
%%reference_type = 'trajectory';
% 
%
%
%
%% Definition de la {position desiree de l'imageur par rapport à l'objet
%timg_obj_ref_ini = [0;0;-50];
%thetax_ref_ini = pi/8; %pi/8; %pi;
%thetay_ref_ini = pi/8; %pi/4;
%thetaz_ref_ini = pi/8; %pi/8; %pi/4;
%
%
%%timg_obj_ref_ini = [-200;0;1000];
%%thetax_ref_ini = 0; %pi/8; %pi;
%%thetay_ref_ini = pi/8; %pi/4;
%%thetaz_ref_ini = pi/2; %pi/8; %pi/4;
%
%%timg_obj_ref_ini = [0;0;1000];
%%thetax_ref_ini = 0; %pi/8; %pi;
%%thetay_ref_ini = 0; %pi/4;
%%thetaz_ref_ini = pi - pi/4; %pi/8; %pi/4;
%
%
%% ou definition par la matrice homogene
%%Timg_obj_ref_ini = [0.75 -0.58 0.30 3;
%%			0.36 0.75 0.54 -51;
%%			-0.54 -0.30 0.77 779;
%%			0 0 0 1];      % minimum local (cf Tahri 2010)
%%
%%timg_obj_ref_ini = Timg_obj_ref_ini(1:3,4);
%%ang_euler_xyz = r2eulerxyz(Timg_obj_ref_ini(1:3, 1:3));
%%thetax_ref_ini = ang_euler_xyz(1);
%%thetay_ref_ini = ang_euler_xyz(2);
%%thetaz_ref_ini = ang_euler_xyz(3);
%
%%timg_obj_ref_ini = [10;10;-30];
%%thetax_ref_ini = 0; %pi;
%%thetay_ref_ini = 0;
%%thetaz_ref_ini = 0;
%
%% Definition de la position initiale de l'imageur par rapport à l'objet
%timg_obj_ini = [0; 0; -30];
%thetax_ini = pi/8; %0; %pi;
%thetay_ini = -pi/8; %-pi/8; %pi/4; %pi;
%thetaz_ini = 0; %5*pi/8;
%
%%timg_obj_ini = [200; 0; 1000];
%%thetax_ini = pi/8; %0; %pi;
%%thetay_ini = -pi/8; %-pi/8; %pi/4; %pi;
%%thetaz_ini = pi/8; %pi/2; %5*pi/8;
%
%%timg_obj_ini = [0; 0; 1000];
%%thetax_ini = 0; %pi;
%%thetay_ini = 0; %-pi/8; %pi/4; %pi;
%%thetaz_ini = 0; %pi/2; %5*pi/8;
%
%
%%timg_obj_ini = [-10; -10; -50];
%%thetax_ini = 0; %pi;
%%thetay_ini = 0; %-pi/4; %pi;
%%thetaz_ini = pi/4;
%
%% ou definition par la matrice homogene
%%Timg_obj_ini = [0.66 0.68 -0.29 209;
%%				-0.73 0.66 -0.11 184;
%%					0.11 0.29 0.94 1149;
%%					0 0 0 1];     % minimum local (cf Tahri 2010)
%%
%%
%%timg_obj_ini = Timg_obj_ini(1:3,4);
%%ang_euler_xyz = r2eulerxyz(Timg_obj_ini(1:3, 1:3));
%%thetax_ini = ang_euler_xyz(1);
%%thetay_ini = ang_euler_xyz(2);
%%thetaz_ini = ang_euler_xyz(3);
%
%
%
%
%% Definition de la position finale  de l'imageur par rapport à l'objet (si mouvement)
%% Permet de definir un changement de consigne 
%% On ne considère pour l'instant que des mouvements rectilignes
%
%if reference_type == 'position  '
%	timg_obj_ref_fin = timg_obj_ref_ini;
%	thetax_ref_fin = thetax_ref_ini;
%	thetay_ref_fin = thetay_ref_ini;
%	thetaz_ref_fin = thetaz_ref_ini;
%	
%	ind_ini_traj = 1;
%	ind_fin_traj = 2;
%
%else % reference_type = 'trajectory'
%	
%	timg_obj_ref_fin = [50;50;500];
%	thetax_ref_fin = 0; %pi;
%	thetay_ref_fin = pi/4;
%	thetaz_ref_fin = pi/4;
%
%	ind_ini_traj = 30;
%	ind_fin_traj = 200;
%
%
%end
%
%%timg_obj_ref_fin = [-10; -10; -50];
%%thetax_ref_fin = 0; %pi;
%%thetay_ref_fin = 0; %-pi/4; %pi;
%%thetaz_ref_fin = pi/4;
%
%
%
%% Definition de la position de l'objet (eye in hand) ou de la camera (eye to hand) p/r à la base du robot
%%tbase_targ = [10; 100; 100];
%tbase_targ_ini = [500; 500; 0]; %[10; 0; 100];
%thetax_base_targ_ini = -pi/2; % pi/2
%thetay_base_targ_ini = pi; % pi/8
%thetaz_base_targ_ini = 0;
%
%
%% Definition de la position finale  de la cible par rapport au robot (si mouvement)
%% Permet de définir un mouvement de la cible dans une configuration eye-in-hand
%% ou une perturbation dans une configuration eye-to-hand 
%% On ne considère pour l'instant que des mouvements rectilignes
%if movement == 'n'
%	tbase_targ_fin = tbase_targ_ini;
%	thetax_base_targ_fin = thetax_base_targ_ini;
%	thetay_base_targ_fin = thetay_base_targ_ini;
%	thetaz_base_targ_fin = thetaz_base_targ_ini;
%
%	ind_ini_traj_targ = 1;
%	ind_fin_traj_targ = 2;
%
%else
%	tbase_targ_fin = [10; 0; 200];
%	thetax_base_targ_fin = -pi/2; % pi/2
%	thetay_base_targ_fin = pi; % pi/8
%	thetaz_base_targ_fin = 0;
%	
%	ind_ini_traj_targ = 10;
%	ind_fin_traj_targ = 400;
%
%end
%
%
%
%
%%tbase_eff = [0;0;100];
%%thetax_base_eff = -pi/2;
%%thetay_base_eff = 0;
%%thetaz_base_eff = 0;
%
%
%% Definition de la position reelle de la camera p/r à l'effecteur (eye in hand) ou de l'objet par rapport à l'effecteur (eye to hand)
%teff_tool_real = [0;0;0]
%thetax_eff_tool_real = 0;
%thetay_eff_tool_real = 0;
%thetaz_eff_tool_real = 0;
%
%% Definition de la position estimee de la camera p/r à l'effecteur (eye in hand) ou de l'objet par rapport à l'effecteur (eye to hand)
%teff_tool = teff_tool_real; % + [-50;0;0]
%thetax_eff_tool = thetax_eff_tool_real;% - 25*pi/180;
%thetay_eff_tool = thetay_eff_tool_real;
%thetaz_eff_tool = thetaz_eff_tool_real;% + 25*pi/180;
%
%
%
%
%% Definition de la camera
%%imager_real.K = [800 0 400; 0 800 300; 0 0 1];
%%imager_real.corner_im = [1 800 1 600]; 
%%%imager_real.corner_im = [1 1600 1 1200]; 
%
%%imager_est.K =  imager_real.K; %[750 0 400; 0 850 300; 0 0 1]; %imager_real.K; %[600 0 400; 0 600 300; 0 0 1]; 
%%imager_est.corner_im = imager_real.corner_im;
%
%% Definition de l'imageur à coupe
%imager_real.K = [2 0 400; 0 2 300; 0 0 1];
%imager_real.corner_im = [1 800 1 600]; 
%
%imager_est.K = imager_real.K; %[600 0 400; 0 600 300; 0 0 1];
%imager_est.corner_im = imager_real.corner_im;
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition de l'asservissement
%% Periode d'echantillonnage
%Te = 0.04;
%
%% gain du correcteur prop
%lambda = 1; %3; %0.5/Te; %1; %1/Te;
%
%% Duree max de l'asservissement
%nb_samp_max = 250; %30; %500;
%
%% erreurs limites pour arret
%lim_average_im_err = 0.2;
%
%% gestion des singularites (cas 2D)
%%sing_val_handling = 'y';
%sing_val_handling = 'n';
%
%
%% Jacobien image utilisé
%%type_jacobian = 'online '; % calcule a chaque pas
%%type_jacobian = 'final  '; % jacobien de la position finale
%%type_jacobian = 'initial'; % jacobien de la position initiale
%type_jacobian = 'known  '; % jacobien connu (purement virtuel)
%%type_jacobian = 'mean   '; % Methode Malis
%%type_jacobian = 'meanmod'; % Methode Malis modifiee par Tahri / Mezouar;  % Attention !! Ne fonctionne que si l'effecteur et la camera / la cible sont confondus !! il faudrait reprendre le code
%
%% Pour les asservissements 3D et 2D utilisant le calcul de la pose pour determiner le jacobien
%pose = 'known  '; 
%%pose = 'unknown';
%
%% Control du robot
%%type_robot_control = 'joint_vel';
%type_robot_control = 'op_vel   ';
%
%% Est-ce que tous les DOFs sont utilises ? A definir pour le
%% cas ou on ne contrôle que certains DOFs
%used_DOFs = [1:6]; 
%
%
%% bruit dans l'extraction des primitives
%noise_ampl = 0; %0.1; %2; %1;
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul automatique des donnees utiles
%
%
%Rimg_obj_ref_ini = theta2r([thetax_ref_ini;0;0]) * theta2r([0;thetay_ref_ini;0]) * theta2r([0;0;thetaz_ref_ini]);
%Timg_obj_ref_ini = [Rimg_obj_ref_ini timg_obj_ref_ini;
%				0 0 0 1]
%
%
%Rimg_obj_ref_fin = theta2r([thetax_ref_fin;0;0]) * theta2r([0;thetay_ref_fin;0]) * theta2r([0;0;thetaz_ref_fin]);
%Timg_obj_ref_fin = [Rimg_obj_ref_fin timg_obj_ref_fin;
%				0 0 0 1];
%
%Timg_obj_ref_fin = Timg_obj_ref_ini;
%
%Rimg_obj_ini = theta2r([thetax_ini;0;0]) * theta2r([0;thetay_ini;0]) * theta2r([0;0;thetaz_ini]);
%Timg_obj_ini = [Rimg_obj_ini timg_obj_ini;
%				0 0 0 1];
%
%
%		
%% Creation de la trajectoire de reference
%nb_pas_mouv = ind_fin_traj - ind_ini_traj +1
%ind = [1:nb_pas_mouv];
%
%delta_trans = (timg_obj_ref_fin - timg_obj_ref_ini)/(nb_pas_mouv-1);
%trans = repmat(timg_obj_ref_ini, 1, nb_pas_mouv) + delta_trans*(ind-1)
%
%
%delta_thetax = (thetax_ref_fin - thetax_ref_ini)/(nb_pas_mouv-1);
%thetax = repmat(thetax_ref_ini, 1, nb_pas_mouv) + (ind-1)*delta_thetax;
%
%delta_thetay = (thetay_ref_fin - thetay_ref_ini)/(nb_pas_mouv-1);
%thetay = repmat(thetay_ref_ini, 1, nb_pas_mouv) + (ind-1)*delta_thetay;
%
%delta_thetaz = (thetaz_ref_fin - thetaz_ref_ini)/(nb_pas_mouv-1);
%thetaz = repmat(thetaz_ref_ini, 1, nb_pas_mouv) + (ind-1)*delta_thetaz;
%
%
%traj_img_obj_ref = [repmat(Timg_obj_ref_ini, 1, ind_ini_traj-1) zeros(4,4*nb_pas_mouv) repmat(Timg_obj_ref_fin, 1, nb_samp_max - ind_fin_traj+1)];
%
%
%ind = (ind_ini_traj-1)*4+1;
%for kk=1:nb_pas_mouv
%	Rimg_obj = theta2r([thetax(kk);0;0]) * theta2r([0;thetay(kk);0])*theta2r([0;0;thetaz(kk)]);
%	traj_img_obj_ref(:, ind:ind+3) = [Rimg_obj trans(1:3,kk); 0 0 0 1];  
%	ind = ind + 4;
%end
%		
%traj_img_obj_ref
%pause;		
%		
%				
%Rbase_targ_ini = theta2r([thetax_base_targ_ini;0;0]) * theta2r([0;thetay_base_targ_ini;0]) * theta2r([0;0;thetaz_base_targ_ini]);
%Tbase_targ_ini = [Rbase_targ_ini tbase_targ_ini;
%				0 0 0 1];
%
%Rbase_targ_fin = theta2r([thetax_base_targ_fin;0;0]) * theta2r([0;thetay_base_targ_fin;0]) * theta2r([0;0;thetaz_base_targ_fin]);
%Tbase_targ_fin = [Rbase_targ_fin tbase_targ_fin;
%				0 0 0 1];
%
%
%% Creation de la trajectoire du mouvement
%nb_pas_mouv = ind_fin_traj_targ - ind_ini_traj_targ +1
%ind = [1:nb_pas_mouv];
%
%delta_trans = (tbase_targ_fin - tbase_targ_ini)/(nb_pas_mouv-1);
%trans = repmat(tbase_targ_ini, 1, nb_pas_mouv) + delta_trans*(ind-1)
%
%
%delta_thetax = (thetax_base_targ_fin - thetax_base_targ_ini)/(nb_pas_mouv-1);
%thetax = repmat(thetax_base_targ_ini, 1, nb_pas_mouv) + (ind-1)*delta_thetax;
%
%delta_thetay = (thetay_base_targ_fin - thetay_base_targ_ini)/(nb_pas_mouv-1);
%thetay = repmat(thetay_base_targ_ini, 1, nb_pas_mouv) + (ind-1)*delta_thetay;
%
%delta_thetaz = (thetaz_base_targ_fin - thetaz_base_targ_ini)/(nb_pas_mouv-1);
%thetaz = repmat(thetaz_base_targ_ini, 1, nb_pas_mouv) + (ind-1)*delta_thetaz;
%
%
%traj_base_targ = [repmat(Tbase_targ_ini, 1, ind_ini_traj_targ-1) zeros(4,4*nb_pas_mouv) repmat(Tbase_targ_fin, 1, nb_samp_max - ind_fin_traj_targ+1)];
%
%
%ind = (ind_ini_traj_targ-1)*4+1;
%for kk=1:nb_pas_mouv
%	Rbase_targ = theta2r([thetax(kk);0;0]) * theta2r([0;thetay(kk);0])*theta2r([0;0;thetaz(kk)]);
%	traj_base_targ(:, ind:ind+3) = [Rbase_targ trans(1:3,kk); 0 0 0 1];  
%	ind = ind + 4;
%end
%		
%traj_base_targ
%pause;		
%		
%
%
%Reff_tool_real = theta2r([thetax_eff_tool_real;0;0]) * theta2r([0;thetay_eff_tool_real;0]) * theta2r([0;0;thetaz_eff_tool_real]);
%Teff_tool_real = [Reff_tool_real teff_tool_real;
%				0 0 0 1];
%
%
%Reff_tool = theta2r([thetax_eff_tool;0;0]) * theta2r([0;thetay_eff_tool;0]) * theta2r([0;0;thetaz_eff_tool]);
%Teff_tool = [Reff_tool teff_tool;
%				0 0 0 1];
%
%
%if configuration == 'eyeinhand'
%	robot.Tbase_eff = Tbase_targ_ini * inv(Timg_obj_ini) * inv(Teff_tool_real);
%else % configuration == 'eyetohand'
%	robot.Tbase_eff = Tbase_targ_ini * Timg_obj_ini * inv(Teff_tool_real);
%end
%
%robot.Rbase_eff = robot.Tbase_eff(1:3,1:3);
%robot.tbase_eff = robot.Tbase_eff(1:3,4);
%
%
%% Calcul de la position articulire correspondant à la position initiale 
%robot.joint_pos = MGI(robot.Tbase_eff)
%pause;
%
%%robot.Tbase_eff = MGD(robot.joint_pos);
%
%
%%Rbase_eff = theta2r([thetax_base_eff;0;0]) * theta2r([0;thetay_base_eff;0]) * theta2r([0;0;thetaz_base_eff]);
%%Tbase_eff = [Rbase_eff tbase_eff;%
%%				0 0 0 1];
%%			
%%robot.Tbase_eff = Tbase_eff;
%%robot.Rbase_eff = Rbase_eff;
%%robot.tbase_eff = tbase_eff;			
%				
%				
%%if configuration == 'eyeinhand'				
%%	Tcam_obj_ini = inv(Teff_tool)*inv(robot.Tbase_eff)*Tbase_targ;
%%else % configuration = 'eyetohand'				
%%	Tcam_obj_ini = inv(Tbase_targ)*robot.Tbase_eff*Teff_tool
%%	pause;		
%%end				
%				
%% Calcul de l'image de reference
%if imaging_device == 'persp_proj'
%	imaging_proj = @camera_proj; 
%
%else
%	imaging_proj = @slice_proj;
%end
%
%[im_pts_ref_ini, liste_visible_points_ref_ini] = imaging_proj(Timg_obj_ref_ini, imager_real, 0); % noise_ampl);
%[im_pts_ref_fin, liste_visible_points_ref_fin] = imaging_proj(Timg_obj_ref_fin, imager_real, 0); %noise_ampl);
%
%if imaging_device == 'persp_proj'
%	nb_info = 2*size(object,2);
%	
%else
%	nb_info = size(object,2);
%	
%end
%traj_im_pts_ref = zeros(nb_info,nb_samp_max);
%		
%for kk=1:nb_samp_max
%	[im_pts, liste_visible_points] = imaging_proj(traj_img_obj_ref(:,4*(kk-1)+1:4*kk), imager_real, 0); % noise_ampl);
%	traj_im_pts_ref(:,kk) = reshape(im_pts(1:2,:), nb_info, 1);
%	traj_liste_visible_points_ref(:,kk) = liste_visible_points;
%end
%
%traj_im_pts_ref
%pause;
%
%% Calcul de l'image initial
%[im_pts_ini, liste_visible_points_ini] = imaging_proj(Timg_obj_ini, imager_real, noise_ampl);
%
%				