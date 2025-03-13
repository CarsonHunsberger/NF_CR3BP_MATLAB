function result = resAApartials(AA,cell)
    result = zeros(6,1);
    x = reshape(AA(1:3),1,3);
    for nn=1:6
        if nn ~= 1 && nn ~= 3
            trigarr = ones(length(cell{nn}{1}),1);
            for k=1:length(cell{nn}{1})
                if cell{nn}{3}(k,1)==1
                    trigarr(k,1) = cos(cell{nn}{4}(k,1)*AA(5));
                end
                if cell{nn}{3}(k,1)==2
                    trigarr(k,1) = sin(cell{nn}{4}(k,1)*AA(5));
                end
            end
            result(nn,1) = sum((cell{nn}{1}.*prod(x.^cell{nn}{2},2)).*trigarr);
        end
    end
end