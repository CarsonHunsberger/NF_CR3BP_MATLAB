function loadNF(mu,N)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% loadNF - Load in normal form data
%   This is an internal function that calls several other internal
%   functions to properly load in and store normal form data for the
%   desired CR3BP mass parameter and normal form truncation degree.
% 
%   Syntax:
%       loadNF(mu,N);
% 
%   Input Arguments:
%       mu - Mass parameter of the CR3BP (m2/(m1+m2))
%       N - Normal form truncation degree (terms of order > N removed)
% 
%   Notes: If mu or N are unspecified, default values of 0.012154 and 11
%          will be taken.
%          
%          If loadNF is not called prior to using the other functions, data
%          corresponding to the default values will be loaded in.
% 
%          The user is encouraged to manually set the default values within
%          this file if they frequently use a particular configuration.
%               
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
    mu double = 0.012154
    N uint8 = 11
end
getNFdata(mu,N,1);
getNFparams(1,mu=mu,N=N);
checkNFloaded(0); %sets flag to true
end