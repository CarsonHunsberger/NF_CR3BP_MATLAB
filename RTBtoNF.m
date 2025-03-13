function NF = RTBtoNF(RTB,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RTBtoNF - Restricted three-body (RTB) to normal form (NF) transform
%   This function takes non-dimensional CR3BP states and transforms them
%   into normal form states.
% 
%   Syntax:
%       NF = RTBtoNF(RTB,Lpt=1,nftype='Birkhoff',method='anl')
% 
%   Input Arguments:
%       RTB - Restricted three-body states (6xN or Nx6)
%               [x,y,z,xdot,ydot,zdot;...]
%       Lpt - Desired libration point (1, 2, or 3)
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
%       method - Transformation method ('anl' or 'num')
%                   'anl': analytical
%                   'num': numerical
% 
%   Output Arguments:
%       NF - Normal form states (6xN or Nx6)
%               [xtilde,ytilde,ztilde,pxtilde,pytilde,pztilde;...]
% 
%   Notes: The numerical transformation is slower but more accurate.
% 
%          If Lpt, nftype, or method are unspecified, default values of 1,
%         'Birkhoff', and 'anl' will be taken, respectively.
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
        RTB {mustBeNonempty}
        options.Lpt int8 = 1
        options.nftype string = 'Birkhoff'
        options.method string = 'anl'
end
Lpt = options.Lpt;
nftype = options.nftype;
method = options.method;
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

[len,width] = size(RTB);
flag = 0;
if width == 6 && ~(len==6)
    flag = 1;
    RTB = RTB';
    width = len;
end
NF = zeros(width,6);

Vinv = data{nf}{Lpt}.Vinv;
T1inv = data{nf}{Lpt}.T1inv;
Cinv = data{nf}{Lpt}.Cinv;
aLpt = data{nf}{Lpt}.aLpt;
anlqp0toqpN = data{nf}{Lpt}.anlqp0toqpN;

GenFuncEOMseval = @(t,x,EOM)[celleval(EOM{1},x);celleval(EOM{2},x);celleval(EOM{3},x);
    celleval(EOM{4},x);celleval(EOM{5},x);celleval(EOM{6},x)];
temparr = 3:N;
GenFuncEOMs = data{nf}{Lpt}.GenFuncEOMs;
optsnum = odeset('RelTol',1e-12,'AbsTol',1e-16);
numqp0toqpN = @(x)numqptransform(x,GenFuncEOMs,GenFuncEOMseval,-1,temparr,optsnum);

if strcmp(method,'num')
    for n=1:width
        NF(n,:) = numqp0toqpN(RTBtoqp(RTB(:,n),Cinv,T1inv,Vinv,aLpt,mu))';
    end
else
    for n=1:width
        NF(n,:) = sixdeval(RTBtoqp(RTB(:,n),Cinv,T1inv,Vinv,aLpt,mu),anlqp0toqpN)';
    end
end
if ~flag
    NF = NF';
end

end