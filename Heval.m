function h = Heval(AA,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Heval - Hamiltonian evaluation
%   This function returns the value of the action-angle Hamiltonian at the
%   specified action-angle state
% 
%   Syntax:
%       h = Heval(AA,Lpt=1,nftype='Birkhoff')
% 
%   Input Arguments:
%       AA - Action-angle state (1xN or Nx1)
%               [I1,I2,I3,phi1,phi2,phi3] (Birkhoff)
%                           or
%               [I1hat,I2hat,I3hat,theta1,theta2,theta3] (Resonant)
%       Lpt - Desired libration point (1, 2, or 3)
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
% 
%   Output Arguments:
%       h - Value of action-angle Hamiltonian (scalar)
% 
%   Note: If Lpt or nftype are unspecified, default values of 1 and 
%         'Birkhoff' will be taken, respectively. Ensure that 
%         the correct nftype is used.
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
H = data{nf}{Lpt}.HAA;
if isResonant
    x = reshape(AA(1:3),[1 3]);
    trigarr = ones(length(H{1}),1);
    for k=1:length(H{1})
        if H{3}(k,1)==1
            trigarr(k,1) = cos(H{4}(k,1)*AA(5));
        end
    end
    h = sum((H{1}.*prod(x.^H{2},2)).*trigarr);
else
    x = reshape(AA(1:3),[1 3]);
    h = sum(H{1}.*prod(x.^H{2},2));
end

end