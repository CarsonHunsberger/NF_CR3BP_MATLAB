function [mu,N] = getNFparams(flag,opts)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% getNFparams - Get normal form parameters
%   This is an internal function that returns the CR3BP mass parameter (mu)
%   and the normal form truncation degree (N) used to generate the current
%   normal form data. While it is an internal function, the user can use it
%   freely with the following syntax.
% 
%   Syntax:
%       [mu,N] = getNFparams();
% 
%   Output Arguments:
%       mu - Mass parameter of the CR3BP (m2/(m1+m2))
%       N - Normal form truncation degree (terms of order > N removed)
% 
%   Note: If this function is called prior to loading in normal form data,
%   NaN values will be returned.
%               
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

arguments
    flag int8 = 0
    opts.mu = nan
    opts.N = nan
end

persistent muval
persistent Nval

if flag
    muval = opts.mu;
    Nval = opts.N;
end

mu = muval;
N = Nval;