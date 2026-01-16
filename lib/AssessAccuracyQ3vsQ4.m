function AssessAccuracyQ3vsQ4()
% ASSESSACCURACYQ3VSQ4 Generates error charts and tables for Q3 vs Q4.
%
%   This function performs a Monte Carlo simulation (30 runs) to compare
%   the accuracy of the Navigation-based method (Q3) vs the Camera-only
%   method (Q4).
%
%   Outputs:
%       - A Boxplot comparing the final TRE of both methods.
%       - Two detailed tables in the Command Window (Intermediate errors & Final TRE).
%
%   Authors: Filippo Morini & Kinan Ali

    %% INITIALIZATION & SETUP
    fprintf('\n----------------------------------------------\n');
    fprintf('STARTING ERROR DECOMPOSITION ANALYSIS (Function Mode)\n');
    fprintf('----------------------------------------------\n');

    % Initial Configuration
    InitConfig();
    
    % Constants & Parameters (Matched to Main Script)
    n_runs = 30;         % Number of Monte Carlo iterations
    n_positions = 4;     % Number of camera positions for Q3
    n_config = 4;        % Number of robot configs for Hand-Eye
    nmove = 8;           % Number of camera movements for Q4
    
    % Known intrinsic parameters
    K_camera = [400 0 380; 0 400 285; 0 0 1];
    K_camera_inv = K_camera \ eye(size(K_camera));
    
    % Robot & Tool Parameters (Assumed global or standard for this lab)
    % We retrieve these from a clean state to ensure accuracy
    T_base_eff = GetRobotCurrentPosition(0); 
    % Note: Teff_inst is usually a global variable in this lab setup. 
    % We access it directly. If strictly local, ensure it's defined.
    global Teff_inst teff_inst Rcam_mark_cam tcam_mark_cam
    
    % Transform for Camera Marker (Q3)
    R_markCam_cam = Rcam_mark_cam';
    t_markCam_cam = -R_markCam_cam * tcam_mark_cam;
    
    % Movement definitions (Same as Main Script)
    move_Q3 = [5 0 0 -8; -2 5 0 6; 0 -6 3 3; 1 1 1 1];
    move_eff = [30 20 5; -20 30 -10; -10 -40 15; 0 0 0]';
    move_t_Q4 = [10 5 -15; 14 6 -8; 8 3 -20; 12 7 -15];
    move_angle_Q4 = [0, -40, 30, 50];


    %% PHASE 1: GROUND TRUTH GENERATION (NOISE = 0)
    fprintf('Computing Ground Truth vectors (Noise = 0)...\n');
    InitConfig();

    % --- A. Ground Truth for Q3 (Target in Nav Frame) ---
    RDs_GT = []; OCs_GT = [];
    for i = 1:n_positions
        measure = GetLocalizerInformation(0);
        T_nav_mark = measure.mark(1).T;
        targ_pos = GetTargetPosition(0);
        m = K_camera_inv * [targ_pos; 1];
        m_nav = T_nav_mark(1:3,1:3) * R_markCam_cam * m;
        RDs_GT = [RDs_GT, m_nav/norm(m_nav)];
        OCs_GT = [OCs_GT, T_nav_mark * [t_markCam_cam;1]];
        MoveCameraInitialPosition();
        MoveCamera([move_Q3(1,i),move_Q3(2,i), move_Q3(3,i)]', move_Q3(4,i))
    end
    target_nav_GT = geometric_triangulation(RDs_GT(1:3,:), OCs_GT(1:3,:));

    % --- B. Ground Truth for Hand-Eye & Trocar ---
    Tbaseeff = GetRobotCurrentPosition(0);
    curr_pos = Tbaseeff(1:3,4);
    T_base_eff_list_GT = zeros(4,4,n_config);
    T_nav_markerInst_list_GT = zeros(4,4,n_config);
    T_base_inst_list_GT = zeros(4,4,n_config);

    for n = 1:n_config
        measure = GetLocalizerInformation(0);
        T_nav_markerInst_list_GT(:,:,n) = measure.mark(2).T;
        T_base_eff_temp = GetRobotCurrentPosition(0);
        T_base_eff_list_GT(:,:,n) = T_base_eff_temp; % Store for Trocar GT
        T_base_inst_list_GT(:,:,n) = T_base_eff_temp*Teff_inst;
        MoveEffPosition(curr_pos + move_eff(:,n));
    end

    % Solve Hand-Eye GT
    l_Akl = zeros(4,4,6); l_Bkl = zeros(4,4,6);
    l_Akl(:,:,1) = inv(T_base_inst_list_GT(:,:,1))*T_base_inst_list_GT(:,:,2);
    l_Akl(:,:,2) = inv(T_base_inst_list_GT(:,:,1))*T_base_inst_list_GT(:,:,3);
    l_Akl(:,:,3) = inv(T_base_inst_list_GT(:,:,1))*T_base_inst_list_GT(:,:,4);
    l_Akl(:,:,4) = inv(T_base_inst_list_GT(:,:,2))*T_base_inst_list_GT(:,:,3);
    l_Akl(:,:,5) = inv(T_base_inst_list_GT(:,:,2))*T_base_inst_list_GT(:,:,4);
    l_Akl(:,:,6) = inv(T_base_inst_list_GT(:,:,3))*T_base_inst_list_GT(:,:,4);
    l_Bkl(:,:,1) = inv(T_nav_markerInst_list_GT(:,:,1))*T_nav_markerInst_list_GT(:,:,2);
    l_Bkl(:,:,2) = inv(T_nav_markerInst_list_GT(:,:,1))*T_nav_markerInst_list_GT(:,:,3);
    l_Bkl(:,:,3) = inv(T_nav_markerInst_list_GT(:,:,1))*T_nav_markerInst_list_GT(:,:,4);
    l_Bkl(:,:,4) = inv(T_nav_markerInst_list_GT(:,:,2))*T_nav_markerInst_list_GT(:,:,3);
    l_Bkl(:,:,5) = inv(T_nav_markerInst_list_GT(:,:,2))*T_nav_markerInst_list_GT(:,:,4);
    l_Bkl(:,:,6) = inv(T_nav_markerInst_list_GT(:,:,3))*T_nav_markerInst_list_GT(:,:,4);
    T_inst_markerInst_GT = eyehand(l_Akl, l_Bkl);

    % Solve Trocar GT
    OCs_tr = [T_base_eff_list_GT(1:3,4,1), T_base_eff_list_GT(1:3,4,2), T_base_eff_list_GT(1:3,4,3), T_base_eff_list_GT(1:3,4,4)];
    inst_base_tr = [T_base_eff_list_GT(:,:,1)*[teff_inst;1], T_base_eff_list_GT(:,:,2)*[teff_inst;1], ...
                    T_base_eff_list_GT(:,:,3)*[teff_inst;1], T_base_eff_list_GT(:,:,4)*[teff_inst;1]];
    RDs_tr = inst_base_tr(1:3,:) - OCs_tr;
    for j=1:size(RDs_tr,2), RDs_tr(1:3,j) = RDs_tr(1:3,j) / norm(RDs_tr(1:3,j)); end
    trocard_base_GT = geometric_triangulation(RDs_tr, OCs_tr);

    % --- C. Ground Truth for Q4 (Target in Camera Frame) ---
    InitConfig();
    OCs_Q4 = zeros(3,nmove); RDs_Q4 = zeros(3,nmove);
    targ_images = zeros(2,nmove); Tcam_inst = zeros(4,4,nmove); T_cam1_camN = zeros(4,4,nmove-1);
    for ii = 1:nmove
        MoveCameraInitialPosition();
        if mod(ii,2) ~= 0, MoveCamera(move_t_Q4(round(ii/2),:)', move_angle_Q4(round(ii/2)));
        else, MoveCamera((move_t_Q4(round(ii/2),:)' + rand(3,1)), (move_angle_Q4(round(ii/2)) + rand)); end
        targ_images(:,ii) = GetTargetPosition(0);
        Tcam_inst(:,:,ii) = GetInstrumentPosition(0);
        if ii > 1, T_cam1_camN(:,:,ii-1) = Tcam_inst(:,:,1) * inv(Tcam_inst(:,:,ii)); end
        if ii < 2, OCs_Q4(:,ii) = [0 0 0]'; else, OCs_Q4(:,ii) = T_cam1_camN(1:3,4,ii-1); end
    end
    m = K_camera_inv * [targ_images; ones(1, size(targ_images,2))];
    RDs_Q4(:,1) = m(1:3,1)/norm(m(1:3,1));
    for jj = 1:(nmove-1)
        dir = T_cam1_camN(1:3,1:3,jj)*m(1:3,jj+1);
        RDs_Q4(:,jj+1) = dir/ norm(dir);
    end
    target_cam1_GT = geometric_triangulation(RDs_Q4, OCs_Q4);


    %% PHASE 2: MONTE CARLO SIMULATION (NOISE = 1)
    fprintf('Running Monte Carlo Simulation (%d runs)...\n', n_runs);

    % Storage
    err_Q3_TargetTriang = zeros(n_runs,1);
    err_Q3_HandEye      = zeros(n_runs,1);
    err_Trocar          = zeros(n_runs,1);
    err_Q4_TargetTriang = zeros(n_runs,1);
    final_TRE_Q3 = zeros(n_runs,1);
    final_TRE_Q4 = zeros(n_runs,1);

    n_sys = 1; % Noise ON

    for k = 1:n_runs
        InitConfig();
        
        % --- Q3 NOISY RUN ---
        RDs = []; OCs = [];
        for i = 1:n_positions
            measure = GetLocalizerInformation(n_sys);
            T_nav_mark = measure.mark(1).T;
            targ_pos = GetTargetPosition(n_sys);
            m = K_camera_inv * [targ_pos; 1];
            m_nav = T_nav_mark(1:3,1:3) * R_markCam_cam * m;
            RDs = [RDs, m_nav/norm(m_nav)];
            OCs = [OCs, T_nav_mark * [t_markCam_cam;1]];
            MoveCameraInitialPosition();
            MoveCamera([move_Q3(1,i),move_Q3(2,i), move_Q3(3,i)]', move_Q3(4,i))
        end
        target_nav_k = geometric_triangulation(RDs(1:3,:), OCs(1:3,:));
        err_Q3_TargetTriang(k) = norm(target_nav_k - target_nav_GT);
        
        % Calibration & Trocar
        T_nav_markerInst_list = zeros(4,4,n_config);
        T_base_inst_list = zeros(4,4,n_config);
        T_base_eff_list = zeros(4,4,n_config); 
        
        Tbaseeff = GetRobotCurrentPosition(n_sys);
        curr_pos = Tbaseeff(1:3,4);
        
        for n = 1:n_config
            measure = GetLocalizerInformation(n_sys);
            T_nav_markerInst_list(:,:,n) = measure.mark(2).T;
            T_base_eff_current = GetRobotCurrentPosition(n_sys); 
            T_base_eff_list(:,:,n) = T_base_eff_current; 
            T_base_inst_list(:,:,n) = T_base_eff_current * Teff_inst;
            MoveEffPosition(curr_pos + move_eff(:,n));
        end
        
        l_Akl = zeros(4,4,6); l_Bkl = zeros(4,4,6);
        l_Akl(:,:,1) = inv(T_base_inst_list(:,:,1))*T_base_inst_list(:,:,2);
        l_Akl(:,:,2) = inv(T_base_inst_list(:,:,1))*T_base_inst_list(:,:,3);
        l_Akl(:,:,3) = inv(T_base_inst_list(:,:,1))*T_base_inst_list(:,:,4);
        l_Akl(:,:,4) = inv(T_base_inst_list(:,:,2))*T_base_inst_list(:,:,3);
        l_Akl(:,:,5) = inv(T_base_inst_list(:,:,2))*T_base_inst_list(:,:,4);
        l_Akl(:,:,6) = inv(T_base_inst_list(:,:,3))*T_base_inst_list(:,:,4);
        l_Bkl(:,:,1) = inv(T_nav_markerInst_list (:,:,1))*T_nav_markerInst_list (:,:,2);
        l_Bkl(:,:,2) = inv(T_nav_markerInst_list (:,:,1))*T_nav_markerInst_list (:,:,3);
        l_Bkl(:,:,3) = inv(T_nav_markerInst_list (:,:,1))*T_nav_markerInst_list (:,:,4);
        l_Bkl(:,:,4) = inv(T_nav_markerInst_list (:,:,2))*T_nav_markerInst_list (:,:,3);
        l_Bkl(:,:,5) = inv(T_nav_markerInst_list (:,:,2))*T_nav_markerInst_list (:,:,4);
        l_Bkl(:,:,6) = inv(T_nav_markerInst_list (:,:,3))*T_nav_markerInst_list (:,:,4);
        
        T_inst_markerInst_k = eyehand(l_Akl, l_Bkl);
        err_Q3_HandEye(k) = norm(T_inst_markerInst_k(1:3,4) - T_inst_markerInst_GT(1:3,4));
        
        OCs_tr = [T_base_eff_list(1:3,4,1), T_base_eff_list(1:3,4,2), T_base_eff_list(1:3,4,3), T_base_eff_list(1:3,4,4)];
        inst_base_tr = [T_base_eff_list(:,:,1)*[teff_inst;1], T_base_eff_list(:,:,2)*[teff_inst;1], ...
                        T_base_eff_list(:,:,3)*[teff_inst;1], T_base_eff_list(:,:,4)*[teff_inst;1]];
        RDs_tr = inst_base_tr(1:3,:) - OCs_tr;
        for j=1:size(RDs_tr,2), RDs_tr(1:3,j) = RDs_tr(1:3,j) / norm(RDs_tr(1:3,j)); end
        trocard_base_k = geometric_triangulation(RDs_tr, OCs_tr);
        err_Trocar(k) = norm(trocard_base_k - trocard_base_GT);
        
        % Final Move Q3
        InitConfig();
        measure = GetLocalizerInformation(n_sys);
        T_nav_markerInst = measure.mark(2).T;
        T_base_eff = GetRobotCurrentPosition(n_sys);
        target_base = T_base_eff * Teff_inst * T_inst_markerInst_k * inv(T_nav_markerInst) * [target_nav_k;1];
        traj = (target_base(1:3) - trocard_base_k) / norm(target_base(1:3) - trocard_base_k);
        MoveEffPosition(target_base(1:3) - traj*350);
        final_TRE_Q3(k) = ComputeTRE();
        
        % --- Q4 NOISY RUN ---
        InitConfig();
        OCs = zeros(3,nmove); RDs = zeros(3,nmove);
        targ_images = zeros(2,nmove); Tcam_inst = zeros(4,4,nmove); T_cam1_camN = zeros(4,4,nmove-1);
        for ii = 1:nmove
            MoveCameraInitialPosition();
            if mod(ii,2) ~= 0, MoveCamera(move_t_Q4(round(ii/2),:)', move_angle_Q4(round(ii/2)));
            else, MoveCamera((move_t_Q4(round(ii/2),:)' + rand(3,1)), (move_angle_Q4(round(ii/2)) + rand)); end
            targ_images(:,ii) = GetTargetPosition(n_sys);
            Tcam_inst(:,:,ii) = GetInstrumentPosition(n_sys);
            if ii > 1, T_cam1_camN(:,:,ii-1) = Tcam_inst(:,:,1) * inv(Tcam_inst(:,:,ii)); end
            if ii < 2, OCs(:,ii) = [0 0 0]'; else, OCs(:,ii) = T_cam1_camN(1:3,4,ii-1); end
        end
        m = K_camera_inv * [targ_images; ones(1, size(targ_images,2))];
        RDs(:,1) = m(1:3,1)/norm(m(1:3,1));
        for jj = 1:(nmove-1)
            dir = T_cam1_camN(1:3,1:3,jj)*m(1:3,jj+1);
            RDs(:,jj+1) = dir/ norm(dir);
        end
        target_cam1_k = geometric_triangulation(RDs, OCs);
        err_Q4_TargetTriang(k) = norm(target_cam1_k - target_cam1_GT);
        
        T_base_eff = GetRobotCurrentPosition(n_sys);
        Tcam_inst1 = Tcam_inst(:,:,1);
        target_base = T_base_eff * Teff_inst * inv(Tcam_inst1) * [target_cam1_k;1];
        traj = (target_base(1:3) - trocard_base_k) / norm(target_base(1:3) - trocard_base_k);
        MoveCameraInitialPosition();
        MoveEffPosition(target_base(1:3) - traj*350);
        final_TRE_Q4(k) = ComputeTRE();
    end


    %% PHASE 3: DISPLAY RESULTS
    fprintf('\n===========================================================\n');
    fprintf('     INTERMEDIATE ERROR DECOMPOSITION (Mean +- Std) [mm]   \n');
    fprintf('===========================================================\n');
    fprintf('| Error Source (Intermediate)      | Mean Error | Std Dev |\n');
    fprintf('|----------------------------------|------------|---------|\n');
    fprintf('| Q3: Target Triangulation (Nav)   | %8.4f   | %6.4f  |\n', mean(err_Q3_TargetTriang), std(err_Q3_TargetTriang));
    fprintf('| Q3: Hand-Eye Calib (Translation) | %8.4f   | %6.4f  |\n', mean(err_Q3_HandEye), std(err_Q3_HandEye));
    fprintf('| Shared: Trocar Estimation        | %8.4f   | %6.4f  |\n', mean(err_Trocar), std(err_Trocar));
    fprintf('| Q4: Target Triangulation (Cam)   | %8.4f   | %6.4f  |\n', mean(err_Q4_TargetTriang), std(err_Q4_TargetTriang));
    fprintf('===========================================================\n');

    fprintf('\n===========================================================\n');
    fprintf('           FINAL TRE COMPARISON (End Effector) [mm]        \n');
    fprintf('===========================================================\n');
    fprintf('| Method                           | Mean TRE   | Max TRE |\n');
    fprintf('|----------------------------------|------------|---------|\n');
    fprintf('| Q3 (Nav System)                  | %8.4f   | %6.4f  |\n', mean(final_TRE_Q3), max(final_TRE_Q3));
    fprintf('| Q4 (Camera Only)                 | %8.4f   | %6.4f  |\n', mean(final_TRE_Q4), max(final_TRE_Q4));
    fprintf('===========================================================\n');

    % Boxplot Figure
    figure('Name', 'TRE Comparison Q3 vs Q4', 'Color', 'w');
    group = [ones(size(final_TRE_Q3)); 2*ones(size(final_TRE_Q4))];
    boxplot([final_TRE_Q3; final_TRE_Q4], group, 'Labels', {'Q3 (Nav System)', 'Q4 (Camera Only)'});
    title('Target Registration Error Distribution (Noisy)');
    ylabel('TRE [mm]');
    grid on;

end