function combinations = calcSubComb(comb)
%CALCSUBCOMB Summary of this function goes here
%   Detailed explanation goes here

a = ['G','C','T'];
c = ['T','A','G'];
g = ['A','C','T'];
t = ['C','G','A'];

combinations = '';


for uni = 1:length(comb)
    
    switch comb(uni)

        case 'A'
              combinations(uni,:)  = a;

        case 'C'
              combinations(uni,:)  = c;
     
        case 'G'
              combinations(uni,:)  = g;  

        case 'T'    
              combinations(uni,:)  = t;  
  
    end    
    
end
  combinations = cellstr(combinations);
  combinations = char(combvec((combinations{:}))');

end

