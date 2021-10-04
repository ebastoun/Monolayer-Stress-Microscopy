function min_max_calculation(s_I_max,s_II_max,invariant_I_max,invariant_II_max,J_II_max,LANS_max,MLSS_max,tx_max,ty_max,t_mag_max,...
                             s_I_min,s_II_min,invariant_I_min,invariant_II_min,J_II_min,LANS_min,MLSS_min,tx_min,ty_min,t_mag_min)

%     s_x_maxx              = max(s_x_max);
%     s_y_maxx              = max(s_y_max);
%     s_xy_maxx             = max(s_xy_max);
    s_I_maxx              = max(s_I_max);
    s_II_maxx             = max(s_II_max);
    invariant_I_maxx      = max(invariant_I_max);
    invariant_II_maxx     = max(invariant_II_max);
    J_II_maxx             = max(J_II_max);
    LANS_maxx             = max(LANS_max);
    MLSS_maxx             = max(MLSS_max);
    tx_maxx               = max(tx_max);
    ty_maxx               = max(ty_max);
    t_mag_maxx            = max(t_mag_max);

%     s_x_minn              = min(s_x_min);
%     s_y_minn              = min(s_y_min);
%     s_xy_minn             = min(s_xy_min);
    s_I_minn              = min(s_I_min);
    s_II_minn             = min(s_II_min);
    invariant_I_minn      = min(invariant_I_min);
    invariant_II_minn     = min(invariant_II_min);
    J_II_minn             = min(J_II_min);
    LANS_minn             = min(LANS_min);
    MLSS_minn             = min(MLSS_min);
    tx_minn               = min(tx_min);
    ty_minn               = min(ty_min);
    t_mag_minn            = min(t_mag_min);

    save min_max_values.mat  s_I_maxx s_II_maxx invariant_I_maxx invariant_II_maxx J_II_maxx LANS_maxx MLSS_maxx tx_maxx ty_maxx t_mag_maxx...
                             s_I_minn s_II_minn invariant_I_minn invariant_II_minn J_II_minn LANS_minn MLSS_minn tx_minn ty_minn t_mag_minn
        % s_x_maxx s_y_maxx s_xy_maxx
        % s_x_minn s_y_minn s_xy_minn

end