function [remainder,dataBlock] = CRCSender(dnaBlock)
%CRC Sender part
%   CRC for sender part will perform division with polynomial divider. This
%   function return remainder
polynome = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];
crc_8 = comm.CRCGenerator(polynome); %z^8 + z^7 + z^4 + z^3 + z + 1

dataBlock = crc_8(dnaBlock');
remainder = dataBlock(end-16+1:end);
dataBlock = dataBlock';
end

