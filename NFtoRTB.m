function RTB = NFtoRTB(NF,options)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NFtoRTB - Normal form (NF) to restricted three-body (RTB) transform
%   This function takes normal form states and transforms them into
%   non-dimensional CR3BP states.
% 
%   Syntax:
%       RTB = NFtoRTB(NF,Lpt=1,nftype='Birkhoff',method='anl')
% 
%   Input Arguments:
%       NF - Normal form states (6xN or Nx6)
%               [xtilde,ytilde,ztilde,pxtilde,pytilde,pztilde;...]
%       Lpt - Desired libration point (1, 2, or 3)
%       nftype - Normal form type ('Birkhoff' or 'Resonant')
%       method - Transformation method ('anl' or 'num')
%                   'anl': analytical
%                   'num': numerical
% 
%   Output Arguments:
%       RTB - Restricted three-body states (6xN or Nx6)
%               [x,y,z,xdot,ydot,zdot;...]
% 
%   Notes: The analytical and numerical transformations are similar in
%          accuracy, but the numerical approach is significantly slower.
% 
%          If Lpt, nftype, or method are unspecified, default values of 1,
%         'Birkhoff', and 'anl' will be taken, respectively. Ensure that
%          the correct nftype is used.
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
        NF {mustBeNonempty}
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

[len,width] = size(NF);
flag = 0;
if width == 6 && ~(len==6)
    flag = 1;
    NF = NF';
    width = len;
end
RTB = zeros(width,6);

V = data{nf}{Lpt}.V;
T1 = data{nf}{Lpt}.T1;
C = data{nf}{Lpt}.C;
aLpt = data{nf}{Lpt}.aLpt;
anlqpNtoqp0 = data{nf}{Lpt}.anlqpNtoqp0;

GenFuncEOMseval = @(t,x,EOM)[celleval(EOM{1},x);celleval(EOM{2},x);celleval(EOM{3},x);
    celleval(EOM{4},x);celleval(EOM{5},x);celleval(EOM{6},x)];
temparr = 3:N;
GenFuncEOMs = data{nf}{Lpt}.GenFuncEOMs;
optsnum = odeset('RelTol',1e-12,'AbsTol',1e-16);
numqpNtoqp0 = @(x)numqptransform(x,GenFuncEOMs,GenFuncEOMseval,1,flip(temparr),optsnum);

if strcmp(method,'num')
    for n=1:width
        RTB(n,:) = qptoRTB(numqpNtoqp0(NF(:,n)),V,T1,C,aLpt,mu)';
    end
else
    for n=1:width
        RTB(n,:) = qptoRTB(sixdeval(NF(:,n),anlqpNtoqp0),V,T1,C,aLpt,mu)';
    end
end
if ~flag
    RTB = RTB';
end

end