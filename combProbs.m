classdef combProbs
    %COMBPROBS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        prob {mustBeFloat}
        nucleotids
    end
    
    methods
        function obj = combProbs (prob,nucleotids)
                obj.prob = prob;
                obj.nucleotids = nucleotids;
        end
    end    
end
