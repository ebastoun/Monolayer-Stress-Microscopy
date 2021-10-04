function [s_I_max ,s_II_max ,invariant_I_max,invariant_II_max,J_II_max,LANS_max ,MLSS_max ,tx_max,ty_max,t_mag_max,...
          s_I_min ,s_II_min ,invariant_I_min,invariant_II_min,J_II_min,LANS_min ,MLSS_min ,tx_min,ty_min,t_mag_min,...
          s_I     ,s_II     ,invariant_I    ,invariant_II    ,J_II    ,LANS     ,MLSS     ,tx    ,ty    ,t_mag    ,...
          mean_s_I,mean_s_II,mean_inv_I     ,                          mean_LANS,mean_MLSS,                        ...       
          s_Ix,s_Iy,s_IIx,s_IIy] = initialize_variables(n_images)


    % s_x_max              = zeros(n_images,1);
    % s_y_max              = zeros(n_images,1);
    % s_xy_max             = zeros(n_images,1);
    s_I_max              = zeros(n_images,1);
    s_II_max             = zeros(n_images,1);
    invariant_I_max      = zeros(n_images,1);
    invariant_II_max     = zeros(n_images,1);
    J_II_max             = zeros(n_images,1); 
    LANS_max             = zeros(n_images,1);
    MLSS_max             = zeros(n_images,1);
    tx_max               = zeros(n_images,1);
    ty_max               = zeros(n_images,1);
    t_mag_max            = zeros(n_images,1);
    mean_s_I             = zeros(n_images,1);
    mean_s_II            = zeros(n_images,1);
    mean_inv_I           = zeros(n_images,1);
    mean_LANS            = zeros(n_images,1);
    mean_MLSS            = zeros(n_images,1);


    % s_x_min              = zeros(n_images,1);
    % s_y_min              = zeros(n_images,1);
    % s_xy_min             = zeros(n_images,1);
    s_I_min              = zeros(n_images,1);
    s_II_min             = zeros(n_images,1);
    invariant_I_min      = zeros(n_images,1);
    invariant_II_min     = zeros(n_images,1);
    J_II_min             = zeros(n_images,1); 
    LANS_min             = zeros(n_images,1);
    MLSS_min             = zeros(n_images,1);
    tx_min               = zeros(n_images,1);
    ty_min               = zeros(n_images,1);
    t_mag_min            = zeros(n_images,1);

    s_I           = cell(n_images,1);
    s_II          = cell(n_images,1);
    s_Ix          = cell(n_images,1);
    s_Iy          = cell(n_images,1);
    s_IIx         = cell(n_images,1);
    s_IIy         = cell(n_images,1);
    LANS          = cell(n_images,1);
    MLSS          = cell(n_images,1);
    invariant_I   = cell(n_images,1);
    invariant_II  = cell(n_images,1);
    J_II          = cell(n_images,1);
    tx            = cell(n_images,1);
    ty            = cell(n_images,1);
    t_mag         = cell(n_images,1);


end