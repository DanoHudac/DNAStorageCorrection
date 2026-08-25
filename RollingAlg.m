function [rolledDnaBlock,errors] = RollingAlg(dnaBlock)
%ROLLINGALG 
% Function for generate errors on string of nucleotids
errorNum = 0;
delNum = 0;
insNum = 0;
subNum = 0;

subATOC = 0.160447;
subATOG = 0.708955;
subATOT = 0.123134;

subCTOA = 0.26291;
subCTOG = 0.11737;
subCTOT = 0.610328;

subGTOA = 0.625;
subGTOC = 0.1198;
subGTOT = 0.24485;

subTTOA = 0.1070;
subTTOC = 0.79510;
subTTOG = 0.091743;

i = 1;

while i < size(dnaBlock,2) + 1
    prob = rand;
        if  prob  <= 0.0687
           errorNum = errorNum + 1;
           prob = rand;

           if prob <= 0.42940 
              delNum = delNum +1;
              rolledDnaBlock(i) = 'D';

           elseif prob > 0.42940 && prob <= 0.74381
              insNum = insNum +1;
              insProb = GenProbRandi(1,4,1,1);

              switch insProb

                  case 1
                      rolledDnaBlock(i) = 'A';
                      rolledDnaBlock(i+1) = dnaBlock(i);
                      dnaBlock(:,i+2:end+1) = dnaBlock(:,i+1:end);
                      i = i+1;
                  case 2
                      rolledDnaBlock(i) = 'C';
                      rolledDnaBlock(i+1) = dnaBlock(i);
                      dnaBlock(:,i+2:end+1) = dnaBlock(:,i+1:end);
                      i = i+1;
                  case 3
                      rolledDnaBlock(i) = 'G';
                      rolledDnaBlock(i+1) = dnaBlock(i);
                      dnaBlock(:,i+2:end+1) = dnaBlock(:,i+1:end);
                      i = i+1;
                  case 4
                      rolledDnaBlock(i) = 'T';
                      rolledDnaBlock(i+1) = dnaBlock(i);
                      dnaBlock(:,i+2:end+1) = dnaBlock(:,i+1:end);
                      i = i+1;
              end        
        
           elseif prob > 0.74381
              subNum = subNum +1;
              rolledDnaBlock(i) = dnaBlock(i);
           end        
        
        else
            rolledDnaBlock(i) = dnaBlock(i);

        end    
   i = i+1;
end 


finalBlock = rolledDnaBlock;

    subA = 0;
    subC = 0;
    subG = 0;
    subT = 0;

subCycle = subNum;   

while subCycle > 0
    

    prob = rand;

        if prob <= 0.268
            if ~isempty(find(rolledDnaBlock == 'A'))
            [nuc,index] = SubstitutionNuc(rolledDnaBlock,'A',subATOC,subATOG,'C','G','T');
            finalBlock(index) = nuc;
            subA = subA + 1;
            subCycle = subCycle - 1;
            elseif isempty(find(rolledDnaBlock == 'A'))
            end    
        
        elseif prob > 0.268 && prob <= 0.481
            if ~isempty(find(rolledDnaBlock == 'C'))
            [nuc,index] = SubstitutionNuc(rolledDnaBlock,'C',subCTOA,subCTOG,'A','G','T');
            finalBlock(index) = nuc;
            subC = subC + 1;
            subCycle = subCycle - 1;
            elseif isempty(find(rolledDnaBlock == 'C'))
            end   
    
        elseif prob > 0.481 && prob <= 0.6730
            if ~isempty(find(rolledDnaBlock == 'G'))
            [nuc,index] = SubstitutionNuc(rolledDnaBlock,'G',subGTOA,subGTOC,'A','C','T');
            finalBlock(index) = nuc;
            subG = subG + 1;
            subCycle = subCycle - 1;
            elseif isempty(find(rolledDnaBlock == 'G'))
            end   
    
        elseif prob > 0.6730
            if ~isempty(find(rolledDnaBlock == 'T'))
            [nuc,index] = SubstitutionNuc(rolledDnaBlock,'T',subTTOA,subTTOC,'A','C','G');
            finalBlock(index) = nuc;
            subT = subT + 1;
            subCycle = subCycle - 1;
            elseif isempty(find(rolledDnaBlock == 'T'))
            end   
   
        end    
      
end

finalBlock = strrep(finalBlock,'D','');
finalBlock = strrep(finalBlock," ",'');


disp(errorNum + " errors which is " + (errorNum / (length(dnaBlock)/100)) + "%  of all nucleotids" );
disp(delNum + " of them are deletions that's " + (delNum / (length(dnaBlock)/100)) + "%  of all nucleotids" );
disp(insNum + " of them are insertions that's " + (insNum / (length(dnaBlock)/100)) + "%  of all nucleotids" );
disp(subNum + " of them are substitutions that's " + (subNum / (length(dnaBlock)/100)) + "%  of all nucleotids " + subA + " of As " + subC + " of Cs " + subG + " of Gs " + subT + " of Ts  = " + (subA+subC+subG+subT) );
if subNum == 0 disp("There are no substitutions in the codeword"); elseif subNum > 0 disp((subA/subNum) * 100 + "% sub is for A " + (subC/subNum) * 100 + "% sub is for C "+ (subG/subNum) * 100  + "% sub is for G " + (subT/subNum) * 100 + "% sub is for T"); end

errors = [delNum insNum subNum];
rolledDnaBlock = finalBlock;



end

