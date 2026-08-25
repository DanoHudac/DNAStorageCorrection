function err = calcErrors(errorCount,n)
%CALCERRORS Summary of this function goes here
%   Detailed explanation goes here
probD = 0.0295;
probI = 0.0216;
probS = 0.0175;


err = MyErrors.empty;

perms = [];
for i = 0:errorCount
    for j = 0:(errorCount-i)
        k = errorCount - i - j;
        perms = [perms; i j k];
    end
end

for i = 1:length(perms)

    tempVar = ([probD^perms(i,1) probI^perms(i,2) probS^perms(i,3)]);
    tempProb = prod(tempVar,'all');
    tempProb = tempProb * nchoosek(n,errorCount);
    err(size(err,2) + 1) = MyErrors(tempProb,errorCount,perms(i,1),perms(i,2),perms(i,3));

end    
end

