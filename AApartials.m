function partials = AApartials(AA,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% AApartials - Action-angle (AA) partials
%   This function takes in an action-angle state and returns the time
%   derivative of the state.
% 
%   Syntax:
%       AA = AApartials(AA,Lpt=1,nftype='Birkhoff')
% 
%   Input Arguments:
%       AA - Action-angle state (6x1 or 1x6)
%               [I1,I2,I3,phi1,phi2,phi3] (Birkhoff)
%                           or
%               [I1hat,I2hat,I3hat,theta1,theta2,theta3] (Resonant)
%       Lpt - Desired libration point (1, 2, or 3)
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
% 
%   Output Arguments:
%       partials - Time derivative of action-angle state (6x1 or 1x6)
%               d/dt[I1,I2,I3,phi1,phi2,phi3;...] (Birkhoff)
%                           or
%               d/dt[I1hat,I2hat,I3hat,theta1,theta2,theta3;...] (Resonant)
% 
%   Note: If Lpt or nftype are unspecified, default values of 1 and
%         'Birkhoff' will be taken, respectively. Ensure that the form of
%         AA and nftype match.
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
        AA {mustBeVector}
        options.Lpt int8 = 1
        options.nftype string = 'Birkhoff'
end
Lpt = options.Lpt;
nftype = options.nftype;
checkNFloaded();
[mu,N] = getNFparams();
data = getNFdata(mu,N,0);

if strcmp(nftype,'Resonant')
    isResonant = 1;
    nf = 2;
else
    isResonant = 0;
    nf = 1;
end

temp = size(AA);
flag = 1;
if temp(1)==6
    flag = 0;
end

if isResonant
    partials = resAApartials(AA,data{nf}{Lpt}.AApartialscell);
else
    partials = partialfunc(AA,data{nf}{Lpt}.AApartialscell);
end

if flag
    partials = partials';
end

end