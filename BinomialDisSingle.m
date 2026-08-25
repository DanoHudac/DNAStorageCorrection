function [error] = BinomialDisSingle(probOfType,errorType,sizeOfCodeword)

       prob =  nchoosek(sizeOfCodeword,1) * probOfType^1 * (1 - probOfType)^(sizeOfCodeword-1);
       error = MyError(prob,1,errorType);
          
end