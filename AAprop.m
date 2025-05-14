function AA = AAprop(tspan,AA0,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% AAprop - Action-angle (AA) propagation
%   This function takes an initial action-angle state and propagates it
%   according to the simplified canonical equations in the action-angle
%   space. The output is the propagated state evaluated at user-specified
%   times.
% 
%   Syntax:
%       AA = AAprop(tspan,AA0,Lpt=1,nftype='Birkhoff')
% 
%   Input Arguments:
%       tspan - Vector of evaluation times (Nx1 or 1xN)
%       AA0 - Initial action-angle state (6x1 or 1x6)
%               [I1,I2,I3,phi1,phi2,phi3] (Birkhoff)
%                           or
%               [I1hat,I2hat,I3hat,theta1,theta2,theta3] (Resonant)
%       Lpt - Desired libration point (1, 2, or 3)
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
% 
%   Output Arguments:
%       AA - Propagated state at specified times (Nx6)
%               [I1,I2,I3,phi1,phi2,phi3;...] (Birkhoff)
%                           or
%               [I1hat,I2hat,I3hat,theta1,theta2,theta3;...] (Resonant)
% 
%   Note: If Lpt or nftype are unspecified, default values of 1 and
%         'Birkhoff' will be taken, respectively. Ensure that the form of
%         AA0 and nftype match.
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
        tspan {mustBeVector}
        AA0 {mustBeVector}
        options.Lpt int8 = 1
        options.nftype string = 'Birkhoff'
end
Lpt = options.Lpt;
nftype = options.nftype;
checkNFloaded();
[mu,N] = getNFparams();
data = getNFdata(mu,N,0);

opts = odeset('RelTol',1e-13,'AbsTol',1e-22);

if strcmp(nftype,'Resonant')
    isResonant = 1;
    nf = 2;
else
    isResonant = 0;
    nf = 1;
end

if isResonant
    im = imag(AA0(4));
    [~,AA] = ode113(@(t,x)resAApartials(x,data{nf}{Lpt}.AApartialscell),tspan,real(AA0),opts);
    AA(:,4) = AA(:,4)+1i*im;
else
    partials = partialfunc(AA0,data{nf}{Lpt}.AApartialscell);
    AA = partials'.*reshape(tspan,[length(tspan) 1])+reshape(AA0,[1 6]);
end


end