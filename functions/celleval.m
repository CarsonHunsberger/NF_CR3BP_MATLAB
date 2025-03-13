function result = celleval(cell,x)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Inputs: 
%       Cell:   polynomial to be evaluated ({n x 1}{n x numVariables})
%       x:      input ([m x numVariables])
% Output:
%       result: numerical value of polynomial at desired points ([m x 1])  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ~isempty(cell{1})
    result = sum(cell{1}.*prod(reshape(x,1,6).^cell{2},2));
else
    result = 0;
end