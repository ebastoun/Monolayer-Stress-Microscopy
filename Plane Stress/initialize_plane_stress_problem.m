function [B_store,D,bcdof,bcval] = initialize_plane_stress_problem(fcal,TFM_resolution,E,nu,num_pix)

% Plane stress analysis of plates: initialize variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%--------------------------------------------------------------------------
% Plane stress code:
% Code written by : Siva Srinivas Kolukula                                |
%                   Senior Research Fellow                                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |
%                   India                                                 |
% E-mail : allwayzitzme@gmail.com    
%--------------------------------------------------------------------------|
% Adapted to Monolayer Stress Miscroscopy (MSM)
% Intercellular monolayer stress calculation in 2D
% Modified by Raúl Aparicio Yuste
% PhD student, University of Zaragoza
% E-mail : raparicio@unizar.es
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% INPUTS:

%   - fcal             :    calibration factor [um/pixel]
%   - TFM_resolution   :    size of traction field matrix from TFM (number 
%                           of pixels along 1 dimension)                       
%   - E                :    young's modulus of the cell monolayer    [Pa]
%   - nu               :    poisson's ratio of the cell monolayer    [-]
%   - num_pix          :    original number of pixels along 1 dimension


% OUTPUTS:

%   - B_store          :    elemental matrix of kinematic equation for 
%                           plane stress
%   - D                :    matrix of material properties for plane stress
%   - bcdof            :    list of nodes where there are boundary
%                           conditions
%   - bcval            :    list of values of boundaries associated with
%                           bcdof
%   - coordinates.dat  :    coordinates for the mesh
%   - nodes.dat        :    nodes for the mesh
%   - stiffness.mat    :    global stiffness matrix (global assembly done)





%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%  Generate mesh: coordinates and nodes
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    

%--------------------------------------------------------------------------
%  Generate coordinates
%--------------------------------------------------------------------------                        
   pix_res = num_pix/TFM_resolution;           % pixel_resolution
   len     = num_pix/pix_res;
   x_pix   = (0:pix_res:num_pix-pix_res)*fcal; % Apply fcal
   x_pix   = x_pix';

   for i= 1:len
       x_pixel(i,:) = x_pix(:);
   end
   y_pixel     = x_pixel';

   X0          = x_pixel(:);
   Y0          = y_pixel(:);

   last_coord  = ones(len,1)*num_pix*fcal;
   add_x       = [x_pix ; last_coord; num_pix*fcal];    % last row  ,  last column 
   add_y       = [last_coord ; x_pix(2:end); num_pix*fcal ; 0 ]; 

   X0          = [X0;add_x];
   Y0          = [Y0;add_y]; 
   coordinates = [X0 Y0];

   save coordinates.dat coordinates -ascii
   
   
%--------------------------------------------------------------------------
%  Generate nodes
%--------------------------------------------------------------------------
   column_1 = (1:len*len);
   column_2 = column_1+1;
   column_3 = column_1+len+1;
   column_4 = column_1+len;

   count = len*len+1;
   % Modificate: fix first row
   for i = len:len:len*len
      column_2(i) = count;
      column_3(i) = count+1;
      count       = count +1;
   end

   % Fix last column 
   last_column = len*len - len;
   count       = len;
   for i = last_column:1:len*len
      column_3(i) = len*len + count;
      count       = count + 1; 
   end

   for i = last_column+1:1:len*len
      column_4(i) = column_3(i)-1;
   end

   % Fix last value (last column, last row)
   column_4(end-len+1) = len*len+len*2+1;

   nodes = [column_1; column_2; column_3; column_4]';
   save nodes.dat nodes -ascii
   
   
   
   
   
   
   
%--------------------------------------------------------------------------  
%--------------------------------------------------------------------------
%  Generate global stiffness matrix
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Variable descriptions                                                                                                 
%   k            = element matrix for stiffness
%   f            = element vector
%   stiffness    = system matrix                                             
%   force        = system vector                                                 
%   displacement = system nodal displacement vector
%   coordinates  = coordinate values of each node
%   nodes        = nodal connectivity of each element
%   index        = a vector containing system dofs associated with each element     
%   gausspoint   = matrix containing sampling points for bending term
%   gaussweight  = matrix containing weighting coefficients for bending term
%   bcdof        = a vector containing dofs associated with boundary conditions     
%   bcval        = a vector containing boundary condition values associated with    
%                  the dofs in 'bcdof'                                              
%   B            = matrix for kinematic equation for plane stress
%   D            = matrix for material property for plane stress

%----------------------------------------------------------------------------  


   disp('Please wait, assembling global stiffness matrix')


%--------------------------------------------------------------------------
% Input data for nodal connectivity for each element
%--------------------------------------------------------------------------
   nel   = length(nodes) ;                % number of elements
   nnel  = 4;                             % number of nodes per element
   ndof  = 2;                             % number of dofs per node (UX,UY)
   nnode = length(coordinates) ;          % total number of nodes in system
   sdof  = nnode*ndof;                    % total system dofs  
   edof  = nnel*ndof;                     % degrees of freedom per element

   nglx  = 2; ngly = 2;                   % 2x2 Gauss-Legendre quadrature 
   %nglxy = nglx*ngly;                     % number of sampling points per element

%--------------------------------------------------------------------------
%  initialization of matrices and vectors
%--------------------------------------------------------------------------

   stiffness    = zeros(sdof,sdof);      % system stiffness matrix
   index        = zeros(edof,1);         % index vector
   B            = zeros(3,edof);         % kinematic matrix for bending
   D            = zeros(3,3);            % constitutive matrix for bending


%--------------------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%--------------------------------------------------------------------------
   [Gausspoint,Gaussweight]=GaussQuadrature(nglx);     % sampling points & weights
   D = E/(1-nu^2)*[1 nu 0 ; nu 1 0; 0 0 (1-nu)/2] ;    % Constituent Matrix for Plane stress

% Compute element stiffness matrix: the size of the elements is the same
% k same for all the elements, only calculated once

   iel=1;           % id element 1

   for i=1:nnel
      nd(i)=nodes(iel,i);            % extract connected node for (iel)-th element
      xx(i)=coordinates(nd(i),1);    % extract x value of the node
      yy(i)=coordinates(nd(i),2);    % extract y value of the node
   end

   k = zeros(edof,edof);        % initialization of stiffness matrix

%--------------------------------------------------------------------------
%  numerical integration for stiffness matrix
%--------------------------------------------------------------------------
   count   = 1;
   B_store = cell(1,nglx*ngly);
   
   for intx=1:nglx
    
      xi = Gausspoint(intx,1);                       % sampling point in x-axis
      wtx = Gaussweight(intx,1);                     % weight in x-axis
    
      for inty=1:ngly
         eta = Gausspoint(inty,1);                   % sampling point in y-axis
         wty = Gaussweight(inty,1) ;                 % weight in y-axis

         [shape,dhdr,dhds] = shapefunctions(xi,eta); % compute shape functions and derivatives at sampling point
         
         jacobian = Jacobian(nnel,dhdr,dhds,xx,yy);  % compute Jacobian
         detjacob=det(jacobian);                     % determinant of Jacobian
         invjacob=inv(jacobian);                     % inverse of Jacobian matrix
         
         [dhdx,dhdy]=shapefunctionderivatives(nnel,dhdr,dhds,invjacob); % derivatives w.r.t. physical coordinate
         
         B=fekineps(nnel,dhdx,dhdy);                 % kinematic matrix for stiffness
         B_store{count} = B;
         count = count + 1;
        
%--------------------------------------------------------------------------
%  compute element stiffness matrix
%--------------------------------------------------------------------------

         k = k+B'*D*B*wtx*wty*detjacob;
 
      end
   end                      % end of numerical integration loop for bending term



%--------------------------------------------------------------------------
%  Assembly
%--------------------------------------------------------------------------

   for iel=1:nel                               % loop for the total number of elements

      nd        = nodes(iel,:);                % extract connected node for (iel)-th element
      index     = elementdof(nd,nnel,ndof);    % extract system dofs associated with element
      stiffness = assemble(stiffness,k,index); % assemble element stiffness matrices 

   end
   
   
   
%--------------------------------------------------------------------------
% Input data for boundary conditions
%--------------------------------------------------------------------------

%%%%%%% uy = 0 lower and upper edge    ux = 0 left and right edge %%%%%%%%%
   
%    border  = (1:2:127);                      % ux = 0 left edge 
%    border  = [border 8193];                  % ux = 0 left edge 
%    border1 = (2:128:8066);                   % uy = 0 lower edge 
%    border1 = [border1 8450];                 % uy = 0 lower edge
%    border2 = (8321:2:8447);                  % ux = 0 right edge 
%    border2 = [border2 8449];                 % ux = 0 right edge
%    border3 = (8194:2:8320);                  % uy = 0 upper edge 
%    border3 = [border3 8448];                 % uy = 0 upper edge
   
%    n_ele   = TFM_resolution;
%    
%    aux     = (n_ele*2)-1;
%    aux0    = (n_ele*n_ele+1)*2 -1;
%    border  = (1:2:aux);                      % ux = 0 left edge 
%    border  = [border aux0];                  % ux = 0 left edge 
%    
%    aux1    =  n_ele*2;
%    aux2    = (n_ele*(n_ele-1)+1)*2;
%    aux3    = (n_ele*n_ele+n_ele+n_ele+1)*2;
%    border1 = (2:aux1:aux2);                  % uy = 0 lower edge 
%    border1 = [border1 aux3];                 % uy = 0 lower edge
%    
%    aux4    = (n_ele*n_ele+n_ele+1)*2-1;
%    aux5    = (n_ele*n_ele+n_ele+n_ele)*2-1;
%    aux6    = aux5+2;
%    border2 = (aux4:2:aux5);                  % ux = 0 right edge 
%    border2 = [border2 aux6];                 % ux = 0 right edge
%    
%    aux7    = (n_ele*n_ele+1)*2;
%    aux8    = (n_ele*n_ele+n_ele)*2;
%    aux9    = (n_ele*n_ele+n_ele+n_ele)*2;
%    border3 = (aux7:2:aux8);                  % uy = 0 upper edge 
%    border3 = [border3 aux9];                 % uy = 0 upper edge
% 
%    bcdof   = [border border1 border2 border3];
    
%%%%%%%%%%%%%%  Encastre: 4 edges  %%%%%%%%%%%%%%%%%%%%%

%    border  = (1:2:127);                      % ux = 0 left edge 
%    border  = [border 8193];                  % ux = 0 left edge 
%    border1 = (2:128:8066);                   % uy = 0 lower edge 
%    border1 = [border1 8450];                 % uy = 0 lower edge
%    border3 = (8194:2:8320);                  % uy = 0 upper edge 
%    border3 = [border3 8448];                 % uy = 0 upper edge
%  
%    bcdof   = [border border1 border3];

   n_ele   = TFM_resolution;
   
   aux     = (n_ele*2)-1;
   aux0    = (n_ele*n_ele+1)*2 -1;
   border  = (1:2:aux);                      % ux = 0 left edge 
   border  = [border aux0];                  % ux = 0 left edge 
   bordera =  border+1;
   
   aux1    =  n_ele*2;
   aux2    = (n_ele*(n_ele-1)+1)*2;
   aux3    = (n_ele*n_ele+n_ele+n_ele+1)*2;
   border1 = (2:aux1:aux2);                  % uy = 0 lower edge 
   border1 = [border1 aux3];                 % uy = 0 lower edge
   borderb =  border1-1;
   
   aux4    = (n_ele*n_ele+n_ele+1)*2-1;
   aux5    = (n_ele*n_ele+n_ele+n_ele)*2-1;
   aux6    = aux5+2;
   border2 = (aux4:2:aux5);                  % ux = 0 right edge 
   border2 = [border2 aux6];                 % ux = 0 right edge
   borderc =  border2+1;
   
   aux7    = (n_ele*n_ele+1)*2;
   aux8    = (n_ele*n_ele+n_ele)*2;
   aux9    = (n_ele*n_ele+n_ele+n_ele)*2;
   border3 = (aux7:2:aux8);                  % uy = 0 upper edge 
   border3 = [border3 aux9];                 % uy = 0 upper edge
   borderd =  border3-1;

   bcdof   = [border bordera border1 borderb border2 borderc border3 borderd];








   bcval   = zeros(1,length(bcdof)) ;
   
   
   % Apply boundary conditions to the global stiffness matrix
   
   [stiffness] = constraints_stiffness(stiffness,bcdof);
   
   %save('stiffness.mat', 'stiffness')
   
   save('stiffness.mat', 'stiffness', '-v7.3') % NEW 

  


end