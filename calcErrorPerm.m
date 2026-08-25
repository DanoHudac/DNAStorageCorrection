function uniques = calcErrorPerm(errorCount)
%CALCERRORPERM Summary of this function goes here
%   Detailed explanation goes here
nucs = {'A', 'C', 'G', 'T'};

nucsTogether = repmat(nucs, 1, errorCount);
comb = nchoosek(nucsTogether, errorCount);
[~,idx] = unique(cell2mat(comb), 'rows');
uniques = comb(idx,:);

end

