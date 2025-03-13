function result = sixdeval(x,cell)
x = x';
result = zeros(size(x));
for k=1:6
    result(:,k) = celleval(cell{k},x);
end
result = result';