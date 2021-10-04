function PlotMesh(coordinates,nodes)
%--------------------------------------------------------------------------
% Code written by : Siva Srinivas Kolukula                                |
%                   Senior Research Fellow                                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |
%                   India                                                 |
% E-mail : allwayzitzme@gmail.com                                         |
%          http://sites.google.com/site/kolukulasivasrinivas/             |    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Finite Element Method Mesh
% Synopsis :
%           PlotMesh(coordinates,nodes)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [node X Y] 
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]    
%--------------------------------------------------------------------------

nel = length(nodes) ;                  % number of elements
nnode = length(coordinates) ;          % total number of nodes in system
nnel = size(nodes,2);                % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;

for iel=1:nel   
     for i=1:nnel
     nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
     X(i,iel)=coordinates(nd(i),1);    % extract x value of the node
     Y(i,iel)=coordinates(nd(i),2);    % extract y value of the node
     end
end
    
% Plotting the FEM mesh, diaplay Node numbers and Element numbers
     
     f1 = figure ;
     set(f1,'name','Mesh','numbertitle','off') ;
     plot(X,Y,'k')
     fill(X,Y,'w')
     hold on
     
     title('Finite Element Mesh') ;
     axis off ;
     k = nodes(:,1:end);
     nd = k' ;
    for i = 1:nel
     %text(X(:,i),Y(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
     %text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
    end 
    
    
    
   %%%%% NEW %%%%%%%%
   %%% LOWER EDGE
%     border  = [2:128:8066];                    % uy = 0 lower edge 
%     %border  = [border 8450];                  % uy = 0 lower edge 
%     nodesss = border/2;
%     
%     plot(X(1,nodesss),Y(1,nodesss),'r*')
    
    %%% LEFT EDGE
%     border  = [1:2:127];                      % ux = 0 left edge 
%     %border  = [border 8193];                  % ux = 0 left edge 
%     nodesss = (border +1)/2;
%     
%     plot(X(1,nodesss),Y(1,nodesss),'r*')
    
    %%% RIGHT EDGE
%     border_p  = [8321:2:8447];
%     %border_p1 = [8449 border_p];
%     nodesss = (border_p +1)/2 - 128;
%     plot(X(3,nodesss),Y(3,nodesss),'r*')
%     
    hold off
    
    
    
    
    
    