function [e_x,e_y,gamma_xy,s_x,s_y,s_xy] =  secondary_variables(D,B_store,displacement,nnode,nel)

    e_x      = zeros(nnode,1);
    e_y      = zeros(nnode,1);
    gamma_xy = zeros(nnode,1);
    s_x      = zeros(nnode,1);
    s_y      = zeros(nnode,1);
    s_xy     = zeros(nnode,1);
    
    nnel  = 4;                          % number of nodes per element
    ndof  = 2;                             % number of dofs per node (UX,UY)
    load nodes.dat ;
    
    for iel=1:nel                       % loop for the total number of elements
        
        nd      = nodes(iel,:);             % extract connected node for (iel)-th element
        index   = elementdof(nd,nnel,ndof); % extract system dofs associated with element
    
        for innel = 1:nnel  
            strain_node     = B_store{innel} * displacement(index);
            e_x (nd(innel))     = strain_node(1);
            e_y (nd(innel))     = strain_node(2);
            gamma_xy(nd(innel)) = strain_node(3);
            stress_node     = D*strain_node;
            s_x (nd(innel)) = stress_node(1);
            s_y (nd(innel)) = stress_node(2);
            s_xy(nd(innel)) = stress_node(3);
        end
    
    end

end