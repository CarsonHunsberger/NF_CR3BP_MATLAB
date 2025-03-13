function NF = AAtoNF(AA,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% AAtoNF - Action-angle (AA) to normal form (NF) transform
%   This function takes action-angle states and transforms them into
%   normal form states.
% 
%   Syntax:
%       NF = AAtoNF(AA,nftype='Birkhoff')
% 
%   Input Arguments:
%       AA - Action-angle states (6xN or Nx6)
%               [I1,I2,I3,phi1,phi2,phi3] (Birkhoff)
%                           or
%               [I1hat,I2hat,I3hat,theta1,theta2,theta3] (Resonant)
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
% 
%   Output Arguments:
%       NF - Normal form states (6xN or Nx6)
%               [xtilde,ytilde,ztilde,pxtilde,pytilde,pztilde;...]
% 
%   Note: If nftype is unspecified, a default value of 'Birkhoff'
%         will be taken. Ensure that the form of AA0 and nftype match.
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
        AA {mustBeNonempty}
        options.nftype string = 'Birkhoff'
end

nftype = options.nftype;

if strcmp(nftype,'Resonant')
    isResonant = 1;
else
    isResonant = 0;
end
[len,width] = size(AA);
flag = 0;
if width == 6 && ~(len==6)
    flag = 1;
    AA = AA';
    width = len;
end
NF = zeros(width,6);
for n=1:width
    NF(n,:) = AAtotilde(AA(:,n),isResonant)';
end
if ~flag
    NF = NF';
end
end