function [displacement,UX,UY] = cell_monolayer_plane_stress(fcal,a,b,elementsalongX,elementsalongY,tx,ty,pixel_resolution,coordinates,nodes,stiffness,bcdof,bcval)

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
% Adapted to Traction Force Microscopy
% Intercellular stress calculation in 2D
% Modified by Raúl Aparicio Yuste
% PhD student, University of Zaragoza
% E-mail : raparicio@unizar.es
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% INPUTS:

%   - fcal             :    Calibration factor                [um/pixel]
%   - a                :    Length of the plate (along X-axes)[um]
%   - b                :    Length of the plate (along Y-axes)[um]
%   - elementsalongX   :    Number of elements along X-axes
%   - elementsalongY   :    Number of elements along Y-axes
%   - tx               :    Traction forces exerted on the susbtrate from
%                           the monolayer [Pa] = [pN/um], X direction
%   - ty               :    Traction forces exerted on the susbtrate from
%                           the monolayer [Pa] = [pN/um], Y direction
%   - pixel_resolution :    pixel size (number of pixels along 1 pixel
%                           dimension)
%   - coordinates.dat  :    coordinates for the mesh
%   - nodes.dat        :    nodes for the mesh
%   - stiffness.mat    :    global stiffness matrix (assembly done)
%   - bcdof            :    list of nodes where there are boundary
%                           conditions
%   - bcval            :    list of values of boundaries associated with
%                           bcdof

% OUTPUTS:

%   - displacement:         nodal displacement of the cell monolayer
%   - UX   :                displacement of the cell monolayer along X
%                           nodes
%   - UY   :                displacement of the cell monolayer along Y
%                           nodes




%--------------------------------------------------------------------------
%  input data 
%--------------------------------------------------------------------------
%disp('Please wait Programme is under Run')
%--------------------------------------------------------------------------
% Input data for nodal coordinate values
%--------------------------------------------------------------------------
%load coordinates.dat ;
%--------------------------------------------------------------------------
% Input data for nodal connectivity for each element
%--------------------------------------------------------------------------
%load nodes.dat ;
%
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

force        = zeros(sdof,1);         % system force vector
%stiffness    = zeros(sdof,sdof);      % system stiffness matrix
displacement = zeros(sdof,1);         % system displacement vector
index        = zeros(edof,1);         % index vector
B            = zeros(3,edof);         % kinematic matrix for bending
D            = zeros(3,3);            % constitutive matrix for bending




%--------------------------------------------------------------------------
% force vector
%--------------------------------------------------------------------------

% Notes:
% Traction forces exerted on the monolayer from the substrate = -tx or
% -ty!!
% Force = Stress * Area ([pN/um]*[um]=[pN])
% Here area = size of the pixel^2

%tx   = ones(64)*150;%   [Pa]=[pN/um2]
%ty   = ones(64)*-322;%   [Pa]=[pN/um2]
tx   =  -tx*fcal^(2)*pixel_resolution^(2);   %[pN]      
ty   =  -ty*fcal^(2)*pixel_resolution^(2);   %[pN]
%tx   = tx*a*b/(elementsalongX*elementsalongY);   %[pN]
%ty   = ty*a*b/(elementsalongX*elementsalongY);   %[pN]
t_x  = zeros(1,elementsalongX*elementsalongX);
t_y  = zeros(1,elementsalongY*elementsalongY);

% Match tx and ty with the mesh element
count = 1;
for j = 1:elementsalongX
    
    for i = elementsalongY:-1:1
        
        t_x(count) = tx(i,j);
        t_y(count) = ty(i,j);
        count = count + 1;
    end
    
    
end

% Distribute each force of the element in 4 nodes
for iel=1:nel                           % loop for the total number of elements


    nd      = nodes(iel,:);             % extract connected node for (iel)-th element
    index   = elementdof(nd,nnel,ndof); % extract system dofs associated with element
    edof    = length(index);
    
    % tx
    for i=1:2:edof
         ii=index(i);
         force(ii)= force(ii)+t_x(iel)/4;
    end
    % ty
    for i=2:2:edof
         ii=index(i);
         force(ii)= force(ii)+t_y(iel)/4;
    end
 
end





%--------------------------------------------------------------------------
% stiffness matrix
%--------------------------------------------------------------------------

%load stiffness.mat

%--------------------------------------------------------------------------
%  apply boundary conditions
%--------------------------------------------------------------------------

[force] = constraints_force(force,bcdof,bcval);


%--------------------------------------------------------------------------
% Solve the matrix equation 
%--------------------------------------------------------------------------
displacement = stiffness\force ;

UX = displacement(1:2:sdof) ;
UY = displacement(2:2:sdof) ;

end