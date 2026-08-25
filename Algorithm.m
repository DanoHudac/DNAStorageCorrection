clc;
clearvars -except data errorPackets;
mapping = ["00", "01", "10", "11"];
nucleotids=['A','C','G','T'];
tStart = 0;
errorPackets = 0;
n = 24;
N = 32; 
count = 1250;
data2 = {char(zeros(count, n)), zeros(count ,3) ,char(zeros(count, n)), zeros(100,1),zeros(1,1)};


for j = 1:count

dnaBlock = int2str(zeros(1,N/2));
dataBlock = [round(rand(1,N))];
[remainder,dataBlock] = CRCSender(dataBlock);

ones = size(strfind( int2str(dataBlock), '1'),2 );
zeros = size(strfind( int2str(dataBlock), '0'),2 );

for i = 1 : 2 : length(dataBlock)
        
	    nucleotid = strrep(int2str(dataBlock(:,i:i+1)), ' ','');

       switch nucleotid

           case '00'
                dnaBlock(i) = 'A';
               
           case '01'
                dnaBlock(i) = 'C';
                
           case '10'
                dnaBlock(i) = 'G';
                
           case '11'
                dnaBlock(i) = 'T';
                             
       end
	 
end 

dnaBlock = strrep(dnaBlock,'0','');
dnaBlock = dnaBlock(~isspace(dnaBlock));
%dnaBlock(21) = '';
%dnaBlock(19) = '';
%dnaBlock(15) = '';
%dnaBlock(13) = '';
 

[rolledDnaBlock,errors] = RollingAlg(dnaBlock);
rolledDnaBlock = convertStringsToChars(rolledDnaBlock(1));

[result, error] = CRCDetector(rolledDnaBlock,mapping,nucleotids);

[infoWord,time] = SelfCorrectAlg(error,rolledDnaBlock);

[~, error2] = CRCDetector(infoWord,mapping,nucleotids);
tStart = tStart + time;

data2{j,1} = dnaBlock;
data2{j,2} = errors;
data2{j,3} = infoWord;
data2{j,4} = strcmp(data2{j,1},data2{j,3});
data2{j,5} = time;



if any(errors > 0) errorPackets = errorPackets + 1; else end

clearvars -except N tStart time data data2 errorPackets mapping nucleotids count n
end

berMin = sum([data2{[data2{:,4}] == 0,2}],2) / (count*2*n) ;
berMax = berMin * 2 ;
berActual = berMin * 1.82;

%tEnd = toc(tStart);
