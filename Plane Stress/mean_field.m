function [mean_s_I,mean_s_II,mean_inv_I,mean_LANS,mean_MLSS] = mean_field(s_I,s_II,invariant_I,LANS,MLSS)

    mean_s_I   = mean(mean(s_I)        );
    mean_s_II  = mean(mean(s_II)       );
    mean_inv_I = mean(mean(invariant_I));
    mean_LANS  = mean(mean(LANS)       );
    mean_MLSS  = mean(mean(MLSS)       );

end