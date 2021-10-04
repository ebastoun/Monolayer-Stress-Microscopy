function [e_x, e_y,gamma_xy] = strain_calculation(u,v, fcal, res)


    e_x      = zeros(size(u,1));
    e_y      = zeros(size(v,1));
    gamma_xy = zeros(size(u,1));
    
    
    % e_x   e_x   e_x   e_x   e_x   e_x   e_x   e_x   e_x   e_x   e_x
    
    for i = 1:size(u,1)              % rows
        
            for j = 1:(size(u,2)-1)  % columns
                
                e_x(i,j) =  ( u(i,j+1) - u(i,j) ) / (fcal*res)   ;
                    
            end
    end
    
    
    % e_y   e_y   e_y   e_y   e_y   e_y   e_y   e_y   e_y   e_y   e_y
    
    for j = 1:size(v,2)              % columns
        
            for i = 1:(size(u,1)-1)  % rows
                
                e_y(i,j) =  ( v(i,j)-v(i+1,j) ) / (fcal*res)  ;
                    
            end
    end
    
    
    % gamma_xy    gamma_xy    gamma_xy    gamma_xy    gamma_xy    gamma_xy
    
    % two contributions
    
    % (1) delta_u/deltay
    
    for j = 1:(size(v,2)-1)          % columns
        
            for i = 1:(size(u,1)-1)  % rows
                
                gamma_xy(i,j) = ( u(i,j)-u(i+1,j) ) / (fcal*res)  ;
                    
            end
    end
    
    % (2) delta_v/deltax
    
    for i = 1:(size(u,1)-1)          % rows
        
            for j = 1:(size(u,2)-1)  % columns
                
                auxiliar       =  (  v(i,j)-v(i,j+1) ) / (fcal*res)  ;
                gamma_xy(i,j)  = gamma_xy(i,j) + auxiliar;
                    
            end
    end





end