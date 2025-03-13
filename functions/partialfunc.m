function  EOMs = partialfunc(x,cell)
EOMs = zeros(6,1);
EOMs(4,1) = celleval(cell{1},x);
EOMs(5,1) = celleval(cell{2},x);
EOMs(6,1) = celleval(cell{3},x);

end