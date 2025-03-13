function flag = checkNFloaded(temp)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% checkNFloaded - Check if normal form data has already been loaded
%   This is an internal function that returns a 1 upon checking that normal
%   form data has been loaded in. If called with a 0 as input, the
%   persistent variable is set to 1 without calling loadNF.
%               
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
arguments
    temp int8 = 1
end
persistent val

if isempty(val)
    if temp
        loadNF();
    end
    val = 1;
end
flag = val;
end