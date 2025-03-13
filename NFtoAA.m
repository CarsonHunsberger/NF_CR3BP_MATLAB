function AA = NFtoAA(NF,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NFtoAA - Normal form (NF) to action-angle (AA) transform
%   This function takes normal form states and transforms them into
%   action-angle states.
% 
%   Syntax:
%       AA = NFtoAA(NF,nftype='Birkhoff')
% 
%   Input Arguments:
%       NF - Normal form states (6xN or Nx6)
%               [xtilde,ytilde,ztilde,pxtilde,pytilde,pztilde;...]
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
% 
%   Output Arguments:
%       AA - Action-angle state (6xN or Nx6)
%               [I1,I2,I3,phi1,phi2,phi3] (Birkhoff)
%                           or
%               [I1hat,I2hat,I3hat,theta1,theta2,theta3] (Resonant)
% 
%   Note: If nftype is unspecified, a default value of 'Birkhoff'
%         will be taken. Ensure that the correct nftype is used.
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
        NF {mustBeNonempty}
        options.nftype string = 'Birkhoff'
end

nftype = options.nftype;

if strcmp(nftype,'Resonant')
    isResonant = 1;
else
    isResonant = 0;
end
[len,width] = size(NF);
flag = 0;
if width == 6 && ~(len==6)
    flag = 1;
    NF = NF';
    width = len;
end
AA = zeros(width,6);
for n=1:width
    AA(n,:) = tildetoAA(NF(:,n),isResonant)';
end
if ~flag
    AA = AA';
end
end