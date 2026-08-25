function [results,error] = CRCDetector(dnaBlock,mapping,nucleotids)
% CRCDetector function
%   This function will provide result from checksum

[~, index] = ismember(dnaBlock, nucleotids);
result = mapping(index);
result = cell2mat(join(result));
result = result(~isspace(result));

results = zeros(1,length(result));
results(result=='1')=1;

crcdetector = comm.CRCDetector([1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1]); %{\displaystyle x^{8}+x^{7}+x^{5}+x^{2}+x+1}
[~,error] = crcdetector(results');

end

