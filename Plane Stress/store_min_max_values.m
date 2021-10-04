function [s_I_max,s_II_max,invariant_I_max,invariant_II_max,J_II_max,LANS_max,MLSS_max,tx_max,ty_max,t_mag_max,...
          s_I_min,s_II_min,invariant_I_min,invariant_II_min,J_II_min,LANS_min,MLSS_min,tx_min,ty_min,t_mag_min]...
          = store_min_max_values(s_I,s_II,invariant_I,invariant_II,J_II,LANS,MLSS,tx,ty,t_mag)

 
    
%     s_x_max              = max(s_x);
%     s_y_max              = max(s_y);
%     s_xy_max             = max(s_xy);
    s_I_max              = max(max(s_I));
    s_II_max             = max(max(s_II));
    invariant_I_max      = max(max(invariant_I));
    invariant_II_max     = max(max(invariant_II));
    J_II_max             = max(max(J_II));
    LANS_max             = max(max(LANS));
    MLSS_max             = max(max(MLSS));
    tx_max               = max(max(tx));
    ty_max               = max(max(ty));
    t_mag_max            = max(max(t_mag));
    
%     s_x_min              = min(s_x);
%     s_y_min              = min(s_y);
%     s_xy_min             = min(s_xy);
    s_I_min              = min(min(s_I));
    s_II_min             = min(min(s_II));
    invariant_I_min      = min(min(invariant_I));
    invariant_II_min     = min(min(invariant_II));
    J_II_min             = min(min(J_II));
    LANS_min             = min(min(LANS));
    MLSS_min             = min(min(MLSS));
    tx_min               = min(min(tx));
    ty_min               = min(min(ty));
    t_mag_min            = 0 ; %min(min(t_mag));

end