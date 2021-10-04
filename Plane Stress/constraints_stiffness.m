function [stiffness]=constraints_stiffness(stiffness,bcdof)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%
%  Synopsis:
%     [stiffness]=constraints(kk,bcdof,bcval)
%
%  Variable Description:
%     stiffness - system matrix before applying constraints 
%     bcdof - a vector containging constrained d.o.f
%-----------------------------------------------------------

 n=length(bcdof);
 sdof=size(stiffness);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       stiffness(c,j)=0;
       stiffness(j,c)=0;  % ONLY WHEN DISPLACEMENT = 0
    end
    
    stiffness(c,c)=1;
 end

end