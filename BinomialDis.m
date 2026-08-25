function [errors] = BinomialDis(probOfType,errorType,sizeOfCodeword,thresh)
%BINOMICALDIV contains a method for calculating the binomial distribution for individual errors

errors = MyError.empty;

 for k = 1 : 1 : sizeOfCodeword
       prob =  nchoosek(sizeOfCodeword,k) * probOfType^k * (1 - probOfType)^(sizeOfCodeword-k);
       
       if prob > thresh
          errors(size(errors,2) + 1) = MyError(prob,k,errorType);
          
       elseif prob <= thresh
          break;    

       end  
 end  

end

