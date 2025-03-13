function RTB = RTBprop(tspan,RTB0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RTBprop - Restricted three-body (RTB) propagation
%   This function takes an initial state in the standard CR3BP frame and
%   propagates it according to the CR3BP equations of motion and returns
%   the propagated state at user-specified times.
% 
%   Syntax:
%       RTB = RTBprop(tspan,RTB0)
% 
%   Input Arguments:
%       tspan - Vector of evaluation times (Nx1 or 1xN)
%       RTB0 - Initial restricted three-body state (6x1 or 1x6)
%               [x,y,z,xdot,ydot,zdot]
% 
%   Output Arguments:
%       RTB - Propagated state at specified times (Nx6)
%               [x,y,z,xdot,ydot,zdot;...]
% 
%   Note: All values should be in non-dimensional units
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
    tspan {mustBeVector}
    RTB0 {mustBeVector}
end
opts = odeset('RelTol',1e-13,'AbsTol',1e-22);
mu = getNFparams();
[~,RTB] = ode113(@(t,x)CR3BP(t,x,mu),tspan,reshape(RTB0,[1 6]),opts);
end