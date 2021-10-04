function [e_I,e_II,e_I_vector_u_iter,e_I_vector_v_iter,e_II_vector_u_iter,e_II_vector_v_iter,invariant_I,invariant_II,J_II] = ...
    principal_components_directions(component1,component2,component3,elementsalongX)

     n_el      = elementsalongX;

     new_image1 = component1(1:n_el*n_el);
     new_image2 = component2(1:n_el*n_el);
     new_image3 = component3(1:n_el*n_el);
     
     new_image1 = reshape(new_image1,n_el,n_el);
     new_image2 = reshape(new_image2,n_el,n_el);
     new_image3 = reshape(new_image3,n_el,n_el);
     
     new_image1 = flip(new_image1);
     new_image2 = flip(new_image2);
     new_image3 = flip(new_image3);
     
     auxa       = component1(n_el*n_el+1:n_el*n_el+n_el);
     auxb       = component2(n_el*n_el+1:n_el*n_el+n_el);
     auxc       = component3(n_el*n_el+1:n_el*n_el+n_el);
     
     new_image1 = [auxa' ; new_image1];
     new_image2 = [auxb' ; new_image2];
     new_image3 = [auxc' ; new_image3];
     
     aux1       = n_el*n_el+n_el+1;
     aux2       = n_el*n_el+n_el+n_el;
     
     add1       = [component1(aux2:-1:aux1)' component1(aux2+1)];
     add2       = [component2(aux2:-1:aux1)' component2(aux2+1)];
     add3       = [component3(aux2:-1:aux1)' component3(aux2+1)];
     
     new_image1 = [new_image1 add1'];
     new_image2 = [new_image2 add2'];
     new_image3 = [new_image3 add3'];
     
 
   
     % To speed up:
    e_matrix           = zeros(2,2);
    e_I_vector_u_iter  = zeros(size(new_image1,2),size(new_image1,1));
    e_I_vector_v_iter  = zeros(size(new_image1,2),size(new_image1,1));
    e_II_vector_u_iter = zeros(size(new_image1,2),size(new_image1,1));
    e_II_vector_v_iter = zeros(size(new_image1,2),size(new_image1,1));
    
    %%%%%%%%%%%%% Calculate principal components (strain rate)%%%%%%%%%%%%%
    
    for i = 1:size(new_image1,2)            % rows
            for j = 1:size(new_image1,1)    % columns
                
                % Create matrix to calculate principal components of each
                % pixel
                e_xy            = new_image3(i,j)/2; 
                e_matrix        = [new_image1(i,j) e_xy; e_xy  new_image2(i,j)];
                [eigvec,eigval] = eigs(e_matrix);
                
                aux1  = eigval(1,1); % eigen value 1
                aux2  = eigval(2,2); % eigen value 2
                
                % eigval does not sort the values:
                
                if aux1 <= aux2
                    e_I(i,j)   = eigval(2,2); % eigen value 1
                    e_II(i,j)  = eigval(1,1); % eigen value 2
                    n1_x       = eigvec(1,2); % nI  component x
                    n1_y       = eigvec(2,2); % nI  component y
                    n2_x       = eigvec(1,1); % nII component x
                    n2_y       = eigvec(2,1); % nII component y
                else
                    e_I(i,j)   = eigval(1,1); % eigen value 1
                    e_II(i,j)  = eigval(2,2); % eigen value 2
                    n1_x       = eigvec(1,1); % nI  component x
                    n1_y       = eigvec(2,1); % nI  component y
                    n2_x       = eigvec(1,2); % nII component x
                    n2_y       = eigvec(2,2); % nII component y
                end
                

                e_I_vector_u_iter  (i,j) = abs(e_I(i,j))  * n1_x; % eI  component x
                e_I_vector_v_iter  (i,j) = abs(e_I(i,j))  * n1_y; % eI  component y
                e_II_vector_u_iter (i,j) = abs(e_II(i,j)) * n2_x; % nII component x
                e_II_vector_v_iter (i,j) = abs(e_II(i,j)) * n2_y; % nII component y
  
            end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Calculate invariants %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Invariants
    invariant_I   = e_I + e_II;
    invariant_II  = e_I.*e_II;
    J_II          = (  ( e_I - (invariant_I/3)  ).^2   +   ( e_II - (invariant_I/3)  ).^2  ).^(0.5);
    
end
     

   