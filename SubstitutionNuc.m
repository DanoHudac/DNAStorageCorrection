function [nuc,index] = SubstitutionNuc(rolledDnaBlock,mainNuc,probMin,probMax,nuc1,nuc2,nuc3)
% SUBSTITUTIONNUC 
% Substitute one nucleotid to another

indOfBlock = find(rolledDnaBlock == mainNuc);
indForSub = randsample(indOfBlock,1);


    prob = rand;

    if prob <= probMin
           nuc = nuc1;
           index = indForSub;
  

    elseif prob > probMin && prob <= probMax + probMin
           nuc = nuc2;
           index = indForSub;
           

    elseif prob > probMin + probMax
           nuc = nuc3;
           index = indForSub;
           

    end    

end

