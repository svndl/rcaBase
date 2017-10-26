function mOut = Zero2NaN(m,dim)
    mOut = m;
    for d=1:length(dim)
        if dim(d)>3
            error('dimensions above 3 not allowed');
        else
        end
        repDims = [1,1,1];
        repDims(dim(d)) = size(m,dim(d));
        nanMat = repmat(nansum(m,dim(d)),repDims);
        idx=find(nanMat==0);
        if ~isempty(idx)
            [k{1},k{2},k{3}] = ind2sub(size(m),idx);
            for n=1:length(k{1}); 
                mOut(k{1}(n),k{2}(n),k{3}(n)) = NaN;
            end
        else
        end
    end
end

