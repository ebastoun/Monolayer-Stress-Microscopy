function [force]=constraints_force(force,bcdof,bcval)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%
%  Synopsis:
%     [force]=constraints(ff,bcdof,bcval)
%
%  Variable Description:
%     force - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing contained value 
%-----------------------------------------------------------
 
 n=length(bcdof);

 for i=1:n
    c=bcdof(i);
    force(c)=bcval(i);
 end
