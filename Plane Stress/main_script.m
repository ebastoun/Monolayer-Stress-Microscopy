%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Plane stress problem adapted to Monolayer Stress Microscopy (MSM)
% Intercellular monolayer stress calculation in 2D
% By Raul Aparicio Yuste
% PhD student, University of Zaragoza & University of Tubingen
% E-mail : raparicio@unizar.es
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Plane stress problem adapted from plane stress analysis of plates
% Original code written by : Siva Srinivas Kolukula                       |
%                            Senior Research Fellow                       |
%                            Structural Mechanics Laboratory              |
%                            Indira Gandhi Center for Atomic Research     |
%                            India   
% E-mail : allwayzitzme@gmail.com    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%--------------------------------------------------------------------------


% Run block by block
clear 
clc
%%
%--------------------------------------------------------------------------
% PARAMETERS: DEFINED BY USER
%--------------------------------------------------------------------------

% UNITS:
%                       Displacement    u[um]
%                       Young's Modulus E[Pa = pN/um^2]
%                       Stress          s[Pa]


% LOAD INPUTS: 
% 2D traction field from Traction Force Microscopy (TFM) (on substrate)

load traction_x_pos12_inf   %[Pa]  Dataset of images, x component
load traction_y_pos12_inf   %[Pa]  Dataset of images, y component


% DEFINE parameters:
fcal              = 0.228;             % Calibration factor      [um/pixel] 
n_images          = 410;               % Total number of images 
time_step         = 10;                % Images taken every x min[min] 
E_cell            = 500;               % Cell Young's modulus    [Pa]
nu_cell           = 0.48;              % Poisson's ratio         [-]

%--------------------------------------------------------------------------

% Other parameters
num_pix        = 1024;                 % Original number of pixels along X or
                                       % Y dimension. This parameter should
                                       % not change                                                     
a              = num_pix*fcal;         % Length of the plate (along X-axes)[um]
b              = num_pix*fcal;         % Length of the plate (along Y-axes)[um]
TFM_resolution = size(traction_x{1},1);% Number of pixels of traction field
                                       % matrix from TFM
elementsalongX = TFM_resolution ;      % Number of elements along X-axes
elementsalongY = TFM_resolution ;      % Number of elements along Y-axes

% WARNING:
% TFM_resolution has to be divisor of 1024, 64 or 32 is recommended. When 
% this variable is increased, plotting will take more time and also the
% computational time will be increased. This variable is chosen before the
% computation of TFM

%--------------------------------------------------------------------------

% Generate coordinates (used by some plots)
pix_res       = num_pix/TFM_resolution;             % pixel_resolution for plots
x_pix         = (1:pix_res:num_pix+1)*fcal;         % Generate X coordinates
y_pix         = (1:pix_res:num_pix+1)*fcal;         % Generate Y coordinates
y_pix         = y_pix';
time          = (0:time_step:n_images*time_step-1); % Time [min]


% Initialize plane stress problem
tic
[B_store,D,bcdof,bcval] = initialize_plane_stress_problem(fcal,TFM_resolution,E_cell,nu_cell,num_pix);
toc

% Plot mesh
load coordinates.dat ;
load nodes.dat ;
PlotMesh(coordinates,nodes)           % Generated mesh for Plane Stress Problem

%% 
%--------------------------------------------------------------------------
% COMPUTE PLANE STRESS PROBLEM: CELL MONOLAYER
%--------------------------------------------------------------------------

% Load global stiffness matrix
load stiffness.mat

% Initialize variables which will be stored after running the code
[s_I_max ,s_II_max ,invariant_I_max,invariant_II_max,J_II_max,LANS_max ,MLSS_max ,tx_max,ty_max,t_mag_max,...
 s_I_min ,s_II_min ,invariant_I_min,invariant_II_min,J_II_min,LANS_min ,MLSS_min ,tx_min,ty_min,t_mag_min,...
 s_I     ,s_II     ,invariant_I    ,invariant_II    ,J_II    ,LANS     ,MLSS     ,tx    ,ty    ,t_mag    ,...
 mean_s_I,mean_s_II,mean_inv_I     ,                          mean_LANS,mean_MLSS,                        ...  
 s_Ix,s_Iy,s_IIx,s_IIy] = initialize_variables(n_images);

% Loop for the calculation of the plane stress problem
tic
parfor k = 1:n_images
    
    % Load traction field
    tx    = traction_x{k};
    ty    = traction_y{k};
    t_mag = (tx.^2+ty.^2).^(0.5);
    
    % Retrieve displacements from the cell monolayer: Solve plane stress problem 
    [displacement,UX,UY] = cell_monolayer_plane_stress(fcal,a,a,elementsalongX,elementsalongY,tx,ty,pix_res,coordinates,nodes,stiffness,bcdof,bcval);

    % Compute secondary variables: Strain and stress
    nnode   = length(coordinates) ;  % total number of nodes in system
    nel     = length(nodes) ;        % number of elements
    [e_x,e_y,gamma_xy,s_x,s_y,s_xy] =  secondary_variables(D,B_store,displacement,nnode,nel);

    %U_MAGN = (UX.^2+UY.^2).^(0.5);
    
    % Compute principal components and orientation (eigenvalues and vectors)
    [s_I{k},s_II{k},s_Ix{k},s_Iy{k},s_IIx{k},s_IIy{k},invariant_I{k},invariant_II{k},J_II{k}] = principal_components_directions(s_x,s_y,s_xy,elementsalongX);

    % Compute stress variables:
    % Local average normal stress (LANS) = (SigmaI + SigmaII) /2
    LANS{k} = (s_I{k} + s_II{k})/2;
    % Maximum local shear stress (MLSS) = (SigmaI - SigmaII) /2
    MLSS{k} = (s_I{k} - s_II{k})/2;
    
    
    % Store max min values (needed to plot later)
    [s_I_max(k),s_II_max(k),invariant_I_max(k),invariant_II_max(k),J_II_max(k),LANS_max(k),MLSS_max(k),tx_max(k),ty_max(k),t_mag_max(k),...
     s_I_min(k),s_II_min(k),invariant_I_min(k),invariant_II_min(k),J_II_min(k),LANS_min(k),MLSS_min(k),tx_min(k),ty_min(k),t_mag_min(k)]...
     = store_min_max_values(s_I{k},s_II{k},invariant_I{k},invariant_II{k},J_II{k},LANS{k},MLSS{k},tx,ty,t_mag);

    % To normalize, compute the mean value of each frame
    [mean_s_I(k),mean_s_II(k),mean_inv_I(k),mean_LANS(k),mean_MLSS(k)] = mean_field(s_I{k},s_II{k},invariant_I{k},LANS{k},MLSS{k});
    
    
    
end

toc

% Save results and max and min values, so there is no need to run this part
% of the code again when you close Matlab
tic
min_max_calculation(s_I_max,s_II_max,invariant_I_max,invariant_II_max,J_II_max,LANS_max,MLSS_max,tx_max,ty_max,t_mag_max,...
                    s_I_min,s_II_min,invariant_I_min,invariant_II_min,J_II_min,LANS_min,MLSS_min,tx_min,ty_min,t_mag_min)
                         
save stress_values.mat  s_I s_II s_Ix s_Iy s_IIx s_IIy invariant_I invariant_II J_II LANS MLSS mean_s_I mean_s_II mean_inv_I mean_LANS mean_MLSS ...
                        s_I_max s_II_max invariant_I_max LANS_max MLSS_max s_I_min s_II_min invariant_I_min LANS_min MLSS_min

toc


%%

%--------------------------------------------------------------------------
% PLOTS
%--------------------------------------------------------------------------
% Warning:write paths for phase contrast and bead images from your computer

load min_max_values.mat;
load stress_values.mat;

%%%%%%%%%%%%%%% Choose variable to plot %%%%%%%%%%%%%%%%%%%%%
% 1: Principal stress (I and II) and invariant I
% 2: Stress variables (LANS and MLSS)
% 3: Eigen vector of principal stresses (I and II)
% 4: Eigen vectors + phase contrast image
% 5: Traction field (input, tractions that the monolayer receives from the
%    susbtrate)
plot_variable = 1;



% LOAD IMAGES: PHASE CONTRAST AND BEADS
%beads            = cell(1,n_images);
phase_cont       = cell(1,n_images);

for k = 1:n_images

%     if k<100      % depending on the number of zeros on the file name
%        name = sprintf('/Users/raparicio/Desktop/June_Borrelia/TIFFS-20210607T115656Z-001/TIFFS/Pos22/beads%03d.tif', k);
%     else
%        name = sprintf('/Users/raparicio/Desktop/June_Borrelia/TIFFS-20210607T115656Z-001/TIFFS/Pos22/beads%01d.tif', k);
%     end
% 
%     
%     beads{k}      = double (   imread(name) );  % fluorescent tracer beads  
%     
    if k<100      % depending on the number of zeros on the file name
       name = sprintf('/Users/raparicio/Desktop/July_borrelia/Pos12/pc%03d.tif', k);
    else
       name = sprintf('/Users/raparicio/Desktop/July_borrelia/Pos12/pc%01d.tif', k);  
    end
    phase_cont{k} = double (   imread(name) );  % phase contrast 
    
end


tic
parfor k = 1:n_images
    
    
    if plot_variable == 1
    
        path = './outputs/s_principal/s_principal_%d.tif';
        title1 = '\sigma_I [Pa]'; title2 = '\sigma_{II} [Pa]'; title3 = 'Invariant_{I} [Pa]'; 
        Plot_3_variables(s_I{k},s_II{k},invariant_I{k},path,k,0,s_I_maxx,s_II_minn,0,invariant_I_minn,invariant_I_maxx,n_images,time/60,title1,title2,title3,plot_variable)
        G(k) = getframe(gcf) ;
    
    elseif plot_variable == 2
    
        path = './outputs/stress_variables/stress_%d.tif';
        title1 = 'Average normal stress [Pa]'; title2 = 'Maximum shear stress [Pa]'; title3 = 'Phase contrast'; 
        Plot_3_variables_1(LANS{k},MLSS{k},phase_cont{k}.^(0.5),path,k,LANS_minn,LANS_maxx,MLSS_minn,MLSS_maxx,J_II_minn,J_II_maxx,n_images,time/60,title1,title2,title3,plot_variable)
        H(k) = getframe(gcf) ;
    
    elseif plot_variable == 3
    
        path = './outputs/eigen_vectors/eigen_vectors%d.tif';
        title1 = '\sigma_I [Pa]'; title2 = '\sigma_{II} [Pa]';
        Plot_quiver_variables(x_pix,y_pix,s_I{k},0,s_I_maxx-200,s_Ix{k},s_Iy{k},s_II{k},s_II_minn+200,0,s_IIx{k},s_IIy{k},path,k,n_images,time/60,title1,title2,plot_variable)
        I(k) = getframe(gcf) ;
    
    elseif plot_variable == 4
    
        path = './outputs/eigen_vectors_pc/eigen_vectors_pc%d.tif';
        title1 = '\sigma_I [Pa]'; title2 = '\sigma_{II} [Pa]';
        Plot_quiver_variables(x_pix,y_pix,phase_cont{k}.^(0.5),0,0,s_Ix{k},s_Iy{k},phase_cont{k}.^(0.5),0,0,s_IIx{k},s_IIy{k},path,k,n_images,time/60,title1,title2,plot_variable)
        J(k) = getframe(gcf) ;
        
    elseif plot_variable == 5
        
        tx    = -traction_x{k};
        ty    = -traction_y{k};
        t_mag = (tx.^2+ty.^2).^(0.5);
        
        path = './outputs/traction_field_input/traction_field%d.tif';
        title1 = 't_x [Pa]'; title2 = 't_y[Pa]'; title3 = 't_{magnitude} [Pa]';
        Plot_3_variables(tx,ty,t_mag,path,k,tx_minn,tx_maxx,ty_minn,ty_maxx,0,t_mag_maxx,n_images,time/60,title1,title2,title3,plot_variable)
        K(k) = getframe(gcf) ;
    
    end
    
    %path = './outputs/U/displacement_%d.tif';
    %Plot_3_variables_displacement(UX,UY,U_MAGN,path,k,elementsalongX)
    %path = './outputs/e/strain_%d.tif';
    %Plot_3_variables_strain(e_x,e_y,gamma_xy,path,k,elementsalongX)
    
    %path = './outputs/s/stress_%d.tif';
    %Plot_3_variables_stress(s_x,s_y,s_xy,path,k,elementsalongX,s_x_minn,s_x_maxx,s_y_minn,s_y_maxx,s_xy_minn,s_xy_maxx,n_images,time)
    % Create a structure for the video. 
    %F(k) = getframe(gcf) ;
    
end


toc


% Make videos

%writerObj = VideoWriter('stress.avi');

if plot_variable == 1
    writerObj       = VideoWriter('principal_stress.avi');
    video_recording = G;
elseif plot_variable == 2
    writerObj       = VideoWriter('stress_variables.avi');
    video_recording = H;
elseif plot_variable == 3
    writerObj       = VideoWriter('eigen_vectors.avi');
    video_recording = I;
elseif plot_variable == 4
    writerObj       = VideoWriter('eigen_vectors_pc.avi');
    video_recording = J;
elseif plot_variable == 5
    writerObj       = VideoWriter('traction_field_input.avi');
    video_recording = K;
end

writerObj.FrameRate = 3; % set the seconds per image

% open the video writer
open(writerObj);

% write the frames to the video
for i=1:length(video_recording)
    % convert the image to a frame
    frame = video_recording(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

% Save results for further analysis
save pos12_inf.mat  s_I s_II s_Ix s_Iy s_IIx s_IIy invariant_I invariant_II J_II LANS MLSS mean_s_I mean_s_II mean_inv_I mean_LANS mean_MLSS time ...
                    s_I_max s_II_max invariant_I_max LANS_max MLSS_max s_I_min s_II_min invariant_I_min LANS_min MLSS_min

%% Crop images

crop_sI   = cell(1,n_images);
crop_sII  = cell(1,n_images);
crop_invI = cell(1,n_images);
crop_LANS = cell(1,n_images);
crop_MLSS = cell(1,n_images);

tic
for k = 1:n_images
    
    % Crop images
    crop_sI         {k} = imcrop(s_I{k}        ,[10 10 46 46]);
    crop_sII        {k} = imcrop(s_II{k}       ,[10 10 46 46]);
    crop_invI       {k} = imcrop(invariant_I{k},[10 10 46 46]);
    crop_LANS       {k} = imcrop(LANS{k}       ,[10 10 46 46]);
    crop_MLSS       {k} = imcrop(MLSS{k}       ,[10 10 46 46]);
   
end
toc

% Compute min max mean again with crop image
new_min_sI   = zeros(1,n_images);
new_min_sII  = zeros(1,n_images);
new_min_invI = zeros(1,n_images);
new_min_LANS = zeros(1,n_images);
new_min_MLSS = zeros(1,n_images);

new_max_sI   = zeros(1,n_images);
new_max_sII  = zeros(1,n_images);
new_max_invI = zeros(1,n_images);
new_max_LANS = zeros(1,n_images);
new_max_MLSS = zeros(1,n_images);

new_mean_sI   = zeros(1,n_images);
new_mean_sII  = zeros(1,n_images);
new_mean_invI = zeros(1,n_images);
new_mean_LANS = zeros(1,n_images);
new_mean_MLSS = zeros(1,n_images);

tic
for k = 1:n_images
    
    new_min_sI     (k) = min(min(crop_sI{k}     ));
    new_min_sII    (k) = min(min(crop_sII{k}    ));
    new_min_invI   (k) = min(min(crop_invI{k}   ));
    new_min_LANS   (k) = min(min(crop_LANS{k}   ));
    new_min_MLSS   (k) = min(min(crop_MLSS{k}   ));
    
    new_max_sI     (k) = max(max(crop_sI{k}     ));
    new_max_sII    (k) = max(max(crop_sII{k}    ));
    new_max_invI   (k) = max(max(crop_invI{k}   ));
    new_max_LANS   (k) = max(max(crop_LANS{k}   ));
    new_max_MLSS   (k) = max(max(crop_MLSS{k}   ));
    
    new_mean_sI    (k) = mean(mean(crop_sI{k}   ));
    new_mean_sII   (k) = mean(mean(crop_sII{k}  ));
    new_mean_invI  (k) = mean(mean(crop_invI{k} ));
    new_mean_LANS  (k) = mean(mean(crop_LANS{k} ));
    new_mean_MLSS  (k) = mean(mean(crop_MLSS{k} ));
   
end
toc
save pos12_inf_crop.mat    crop_sI      crop_sII      crop_invI       crop_LANS       crop_MLSS    ...
                           new_min_sI   new_min_sII   new_min_invI    new_min_LANS    new_min_MLSS ...
                           new_max_sI   new_max_sII   new_max_invI    new_max_LANS    new_max_MLSS ...
                           new_mean_sI  new_mean_sII  new_mean_invI   new_mean_LANS   new_mean_MLSS
                                  
