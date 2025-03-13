function data = getNFdata(muval,Nval,flag)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% getNFdata - Get normal form data
%   This is an internal function that returns the normal form data for the
%   specified values of CR3BP mass parameter and normal form truncation
%   degree.
% 
%   Syntax:
%       data = getNFdata(mu,N)
% 
%   Input Arguments:
%       mu - Mass parameter of the CR3BP (m2/(m1+m2))
%       N - Normal form truncation degree (terms of order > N removed)
% 
%   Output Arguments:
%       data - Normal form data for all available libration points and
%              nftypes.
%           Format: ('Birkhoff': 1, 'Resonant': 2)
%           data{Lpt}{nftype} - All data for a libration point + nftype
%           data{Lpt}{nftype}.C - Diagonalizing transformation
%           data{Lpt}{nftype}.gamma - Distance from Lpt to nearest primary
%           data{Lpt}{nftype}.GenFuncEOMs - EOMs used for numerical
%                       transformations between qp0 and qpN (NF) states
%           data{Lpt}{nftype}.aLpt - Constant used in transformation
%                       between CR3BP and qp0 states
%           data{Lpt}{nftype}.gammascale - Constant used in transformation
%                       between CR3BP and qp0 states
%           data{Lpt}{nftype}.anlqpNtoqp0 - Polynomials used to transform
%                       from qpN (NF) to qp0 states
%           data{Lpt}{nftype}.anlqp0toqpN - Polynomials used to transform
%                       from qp0 to qpN (NF) states
%           data{Lpt}{nftype}.AApartialscell - Polynomials that provide the
%                       partial derivatives of the action-angle Hamiltonian
%           data{Lpt}{nftype}.T1 - Matrix used in transformation
%                       between CR3BP and qp0 states
%           data{Lpt}{nftype}.Cinv- Matrix used in transformation
%                       between CR3BP and qp0 states
%           data{Lpt}{nftype}.Vinv - Matrix used in transformation
%                       between CR3BP and qp0 states
%           data{Lpt}{nftype}.T1inv - Matrix used in transformation
%                       between CR3BP and qp0 states
%           data{Lpt}{nftype}.HAA - Action-angle Hamiltonian
% 
%   Note: If data for the specified mu and N values has not been generated
%         or placed in the correct subfolder within ./data, it will not be
%         loaded in.
%               
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
persistent datacell
persistent mu
if isempty(mu) || flag == 1
    mu = muval;
end
persistent N
if isempty(N) || flag == 1
    N = Nval;
end
if isempty(datacell) || flag == 1
    digits(40);
    if ispc
        slash = '\';
    else
        slash = '/';
    end
    temp =num2str(mu,17);
    temp = temp(3:end);
    filepath = fileparts(mfilename('fullpath'));
    addpath(strcat(filepath,slash,'functions'));
    nftypes = {'Birkhoff','Resonant'};
    datacell = cell(1,2);
    for n = 1:2
        datacell{n} = cell(1,3);
        nftype = nftypes{n};
        for Lpt=1:3
            folderstr = strcat(filepath,slash,'data',slash,temp,slash,'N',num2str(N),slash,nftype,slash,'L',num2str(Lpt));
            foldercheck = exist(folderstr,'dir');
            if foldercheck==7
                datastr = strcat(folderstr,slash);
                datacell{n}{Lpt}.C = load(strcat(datastr,"C.mat")).C;
                datacell{n}{Lpt}.gamma = load(strcat(datastr,"gamma.mat")).gamma;
                datacell{n}{Lpt}.GenFuncEOMs = load(strcat(datastr,"GenFuncEOMs.mat")).GenFuncEOMs;
                if Lpt==1
                    aLpt=1-datacell{n}{Lpt}.gamma;
                    gammascale = 1;
                end
                if Lpt==2
                    aLpt = 1+datacell{n}{Lpt}.gamma;
                    gammascale = 1;
                end
                if Lpt==3
                    aLpt=-datacell{n}{Lpt}.gamma;
                    gammascale = -1;
                end
                datacell{n}{Lpt}.aLpt = aLpt;
                datacell{n}{Lpt}.gammascale = gammascale;
                
                datacell{n}{Lpt}.anlqpNtoqp0 = load(strcat(datastr,"anlqpNtoqp0.mat")).anlqpNtoqp0; 
                datacell{n}{Lpt}.anlqp0toqpN = load(strcat(datastr,"anlqp0toqpN.mat")).anlqp0toqpN;

                datacell{n}{Lpt}.AApartialscell = load(strcat(datastr,"AApartialscell.mat")).AApartialscell;
                J6 = [zeros(3) eye(3); -eye(3) zeros(3)];
                datacell{n}{Lpt}.T1 = datacell{n}{Lpt}.gammascale*datacell{n}{Lpt}.gamma*eye(6);
                datacell{n}{Lpt}.T1(3,3) = datacell{n}{Lpt}.gamma;
                datacell{n}{Lpt}.T1(6,6) = datacell{n}{Lpt}.gamma;
                
                datacell{n}{Lpt}.V = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 1 0 1 0 0;-1 0 0 0 1 0; 0 0 0 0 0 1];
                datacell{n}{Lpt}.Cinv = -J6*datacell{n}{Lpt}.C'*J6;
                datacell{n}{Lpt}.Vinv = inv(datacell{n}{Lpt}.V);
                datacell{n}{Lpt}.T1inv = inv(datacell{n}{Lpt}.T1);

                datacell{n}{Lpt}.HAA = load(strcat(datastr,"HAA.mat")).HAA;
            end
        end
    end
end
data = datacell;
end