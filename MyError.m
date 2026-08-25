classdef MyError
    
    properties
        prob {mustBeFloat}
        errorCount {mustBeNumeric}
        errorType
    end
    
    methods
        function obj = MyError (prob,errorCount,errorType)

                obj.prob = prob;
                obj.errorCount = errorCount;
                obj.errorType = errorType;
        end
    end    
end



