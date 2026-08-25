classdef MyErrors
    %ERROR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        prob {mustBeFloat}
        errorCount {mustBeNumeric}
        errDel  {mustBeNumeric}
        errIns  {mustBeNumeric}
        errSub  {mustBeNumeric}
    end
    
    methods
        function obj = MyErrors (prob,errorCount,errDel,errIns,errSub)

                obj.prob = prob;
                obj.errorCount = errorCount;
                obj.errDel = errDel;
                obj.errIns  = errIns;
                obj.errSub  = errSub;
        end
    end    
end



