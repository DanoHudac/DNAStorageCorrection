function [infoWord,time] = SelfCorrectAlg(error,codeWord)
%SelfCorrectAlg Summary of this function goes here
%   Detailed explanation goes here

thresh = 0.10;
probE = 0.0687;
probD = 0.0295;
probI = 0.0216;
probS = 0.0175;
n = 24;
mapping = ["00", "01", "10", "11"];
nucleotids=['A','C','G','T'];
errors = MyError.empty;
decisionTree = MyError.empty;
initWord=codeWord;
infoWord=codeWord;

tempErrors = BinomialDis(probE,'Error', n, thresh);
decisionTree(size(decisionTree,2) + 1 : size(decisionTree,2) + size(tempErrors,2)) = tempErrors;

tempErrors = BinomialDisSingle(probD,'D',n);
errors(size(errors,2) + 1 : size(errors,2) + size(tempErrors,2)) = tempErrors;
tempErrors = BinomialDisSingle(probI,'I',n);
errors(size(errors,2) + 1 : size(errors,2) + size(tempErrors,2)) = tempErrors;
tempErrors = BinomialDisSingle(probS,'S',n);
errors(size(errors,2) + 1 : size(errors,2) + size(tempErrors,2)) = tempErrors;



tic();
if error == 0 && n == 15
   infoWord = codeWord;
   time = toc();
   toc();
   return

elseif error == 1
             
       while length(decisionTree) >= 1
           [~,probIdxDT] = max([decisionTree.prob]);
           errorCount = decisionTree(probIdxDT).errorCount;
           
              if errorCount == 1
                 [~,probIdx] = max([errors.prob]);
                 tempErrors = errors;

                for count = 1:size(tempErrors,2)
                   if errors(probIdx).errorType == 'D'

                      check = mod(length(initWord),n-1); 
                        if check ~= 0 
                           disp("It was not only 1 Del"); elseif check == 0
                        for i = 1:length(nucleotids)
                            initWord = [nucleotids(i), codeWord];
                            [~,error] = CRCDetector(initWord,mapping,nucleotids);
                            if error ~= 0 
                               for j = 1:length(codeWord)
                                [initWord(j), initWord(j+1)] = deal(initWord(j+1), initWord(j));
                                [~,error] = CRCDetector(initWord,mapping,nucleotids);
                                if error == 0 disp("We found deletion " + initWord(j+1) + " at position " + (j+1) +" initial codeWord is " + initWord ); 
                                    infoWord = initWord;
                                    time = toc();
                                    toc();
                                    return 
                                elseif error ~= 0
                                end
                                end 
                            elseif error == 0 disp("We found deletion " + nucleotids(i) + " at position " + 1 +" initial codeWord is " + initWord ); 
                                   infoWord = initWord; 
                                   time = toc();
                                   toc();
                                   return
                                                  
                            end
                        end 
                        end
                        

                  elseif errors(probIdx).errorType == 'I'
                         check = mod(length(initWord),n+1); if check ~= 0 disp("It was not only 1 Ins"); elseif check == 0
                            for j = 1:length(codeWord)
                                tempLetter = initWord(j);
                                initWord(j) = '';
                                [~,error] = CRCDetector(initWord,mapping,nucleotids);
                                if error ~= 0 
                                   initWord=codeWord; 
                                     
                                elseif error == 0 
                                    disp("We found insertion " + tempLetter + " at position " + j +" initial codeWord is " + initWord ); 
                                    infoWord = initWord;
                                    time = toc();
                                    toc();
                                    return 
                               end
                            end         
                          end
                        initWord=codeWord; 

                  elseif errors(probIdx).errorType == 'S'
                         check = mod(length(initWord),n); if check ~= 0 disp("It was not only 1 Sub"); elseif check == 0
                           for i = 1:length(initWord)
                               initWord=codeWord; 
                               switch initWord(i)
                                    
                                   case 'A'
                                        nucs = ['G','C','T'];
                                        for j = 1:length(nucs)
                                            initWord(i) = nucs(j);
                                            [~,error] = CRCDetector(initWord,mapping,nucleotids);
                                            if error ~= 0 
                                            elseif error == 0 disp("We found substitution at position " + i + " good nuc is: " + nucs(j)); 
                                                infoWord = initWord;
                                                time = toc();
                                                toc();
                                                return 
                                            end
                                        end

                                   case 'C'
                                        nucs = ['T','A','G'];
                                        for j = 1:length(nucs)
                                            initWord(i) = nucs(j);
                                            [~,error] = CRCDetector(initWord,mapping,nucleotids);
                                            if error ~= 0 
                                            elseif error == 0 disp("We found substitution at position " + i + " good nuc is: " + nucs(j)); 
                                                infoWord = initWord; 
                                                time = toc();
                                                toc();
                                                return 
                                            end
                                        end
 
                                   case 'G'
                                        nucs = ['A','C','T'];
                                        for j = 1:length(nucs)
                                            initWord(i) = nucs(j);
                                            [~,error] = CRCDetector(initWord,mapping,nucleotids);
                                            if error ~= 0 
                                            elseif error == 0 disp("We found substitution at position " + i + " good nuc is: " + nucs(j)); 
                                                infoWord = initWord;
                                                time = toc();
                                                toc();
                                                return  
                                            end
                                        end
                                                                   
                                   case 'T'     
                                        nucs = ['C','G','A'];
                                        for j = 1:length(nucs)
                                            initWord(i) = nucs(j);
                                            [~,error] = CRCDetector(initWord,mapping,nucleotids);
                                            if error ~= 0 
                                            elseif error == 0 disp("We found substitution at position " + i + " good nuc is: " + nucs(j)); 
                                                infoWord = initWord;
                                                time = toc();
                                                toc();
                                                return  
                                            end
                                        end                                                                        
                               end                              
                          end 
                          end
                   end
                   initWord = codeWord;   
                   errors(probIdx)=[];
                   [~,probIdx] = max([errors.prob]);
                end
   
              elseif errorCount > 1

                     err = calcErrors(errorCount,n);
                     pool = 1:n;
                     indices = nchoosek(pool,errorCount);
                     permsInd = perms(1:errorCount);
                     

                     while length(err) >= 1
                     [~,probIdxER] = max([err.prob]);
                     errField = ([err(probIdxER).errDel err(probIdxER).errIns err(probIdxER).errSub]);
                     uniques = calcErrorPerm(errorCount);
                    
                     if err(probIdxER).errDel ~= 0 && err(probIdxER).errIns == 0 &&  err(probIdxER).errSub == 0 || err(probIdxER).errDel == 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub == 0 || err(probIdxER).errDel == 0 && err(probIdxER).errIns == 0 &&  err(probIdxER).errSub ~= 0
                        [~,idx] = find(errField ~= 0); 

                        switch idx
                               
                            case 1
                                check = mod(length(initWord), n+(-1*err(probIdxER).errDel + 1*err(probIdxER).errIns));
                                if check ~= 0 elseif check == 0
                                   tempWordDel = char(repmat(0, 1, n));
                                   for  uni = 1:length(uniques) 
                                       for  ids = 1:length(indices)
                                           excludes = ismember(1:length(tempWordDel),indices(ids,:));
                                           tempWordDel(indices(ids,:)) = char(uniques(uni,:))'; 
                                           tempWordDel(~excludes) = initWord;
                                           [~,error] = CRCDetector(tempWordDel,mapping,nucleotids);
                                           if error ~= 0 
                                           elseif error == 0 disp("We have found " + errorCount + " deletions at positions " + indices(ids,:)); 
                                               infoWord = tempWordDel;
                                               time = toc();
                                               toc();
                                               return 
                                           end
                                       end
                                   end
                                                               
                                end    

                            case 2
                                check = mod(length(initWord), n+(-1*err(probIdxER).errDel + 1*err(probIdxER).errIns));
                                if check ~= 0 elseif check == 0
                                   tempWordIns = initWord; 
                                   for ids = 1:length(indices)
                                            tempWordIns(indices(ids,:)) = [];   
                                            [~,error] = CRCDetector(tempWordIns,mapping,nucleotids);
                                            if error ~= 0 
                                            elseif error == 0 disp("We have found " + errorCount + " insertions at positions " + indices(ids,:)); 
                                                infoWord = tempWordIns;
                                                time = toc();
                                                toc();
                                                return  
                                            end
                                        tempWordIns = initWord;
                                   end
                                end 

                            case 3
                                check = mod(length(initWord), n);
                                if check ~= 0 elseif check == 0
                                    for ids = 1:length(indices)
                                        tempWordSub = initWord;
                                        combinations = calcSubComb(tempWordSub(indices(ids,:)));
                                     for uni = 1:length(combinations)                                        
                                            tempWordSub(indices(ids,:)) = combinations(uni,:);      
                                            [~,error] = CRCDetector(tempWordSub,mapping,nucleotids);
                                            if error ~= 0 
                                            elseif error == 0 disp("We have found " + errorCount + " substitutions at positions " + indices(ids,:)); 
                                                infoWord = tempWordSub;
                                                time = toc();
                                                toc();
                                                return 
                                            end                                            
                                     end 
                                   end                                                                
                                end 
                        end    
                        
                     elseif err(probIdxER).errDel ~= 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub == 0 || err(probIdxER).errDel == 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub ~= 0 || err(probIdxER).errDel ~= 0 && err(probIdxER).errIns == 0 &&  err(probIdxER).errSub ~= 0 || err(probIdxER).errDel ~= 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub ~= 0
                            check = mod(length(initWord),n+(-1*err(probIdxER).errDel + 1*err(probIdxER).errIns));

                            if check ~= 0


                            elseif check == 0
                                   uniqueDel = calcErrorPerm(err(probIdxER).errDel);
                                   allIndices= [];
                                   for i = 1:size(indices, 1)
                                       for j = 1:size(permsInd, 1)
                                           allIndices = [allIndices; indices(i, permsInd(j,:))];
                                       end
                                   end
								   
                                 
                                  if err(probIdxER).errDel ~= 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub == 0 %110
                                     for ids = 1:length(allIndices)
                                         tempErrorWord = char(repmat(0, 1, n + errField(2)));                                         
                                         excludes = ismember(1:length(tempErrorWord),allIndices(ids,1:err(probIdxER).errDel));   
                                         
                                         for uni = 1:length(uniqueDel)
                                             tempErrorWord(allIndices(ids,1:err(probIdxER).errDel)) = char(uniqueDel(uni,:))'; 
                                             tempErrorWord(~excludes) = initWord;
                                             tempErrorWord(allIndices(ids,end-err(probIdxER).errIns + 1:end)) = [];
											 [~,error] = CRCDetector(tempErrorWord,mapping,nucleotids);
                                             if error ~= 0 
                                             elseif error == 0 disp("We have found " + errorCount + " combo D + I errors at positions " + allIndices(ids,:)); 
                                                infoWord = tempErrorWord;
                                                time = toc();
                                                toc();
                                                return 
                                             end

                                         end    

                                     end    

                                  elseif err(probIdxER).errDel ~= 0 && err(probIdxER).errIns == 0 &&  err(probIdxER).errSub ~= 0 %101
                                      for ids = 1:length(allIndices)
                                          tempErrorWord = char(repmat(0, 1, n));
										  excludes = ismember(1:length(tempErrorWord),allIndices(ids,1:err(probIdxER).errDel));
										  
                                          for uni = 1:length(uniqueDel)
											  tempErrorWord(allIndices(ids,1:err(probIdxER).errDel)) = char(uniqueDel(uni,:))'; 
                                              tempErrorWord(~excludes) = initWord;
                                              combinations = calcSubComb(tempErrorWord(allIndices(ids,end-err(probIdxER).errSub + 1:end)));
										  
                                              for uniSub = 1:length(combinations)
												  tempErrorWord(allIndices(ids,end-err(probIdxER).errSub + 1:end)) = combinations(uniSub,:); 
												  [~,error] = CRCDetector(tempErrorWord,mapping,nucleotids);
												  if error ~= 0
                                                  elseif error == 0 disp("We have found " + errorCount + " combo D + S errors at positions " + allIndices(ids,:)); 
                                                     infoWord = tempErrorWord;
                                                     time = toc();
                                                     toc();
                                                     return 
                                                  end

                                              end 

                                          end 
                                     
                                      end 

                                  elseif err(probIdxER).errDel == 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub ~= 0 %011
                                      for ids = 1:length(allIndices)
                                          tempErrorWord = initWord;
										  combinations = calcSubComb(tempErrorWord(allIndices(ids,end-err(probIdxER).errSub + 1:end)));
                                          for uniSub = 1:length(combinations)
										      tempErrorWord(allIndices(ids,1:err(probIdxER).errIns)) = [];
											  tempErrorWord(allIndices(ids,end-err(probIdxER).errSub + 1:end)) = combinations(uniSub,:);
											  [~,error] = CRCDetector(tempErrorWord,mapping,nucleotids);
											  if error ~= 0 
                                                 tempErrorWord = initWord;
                                              elseif error == 0 disp("We have found " + errorCount + " combo I + S errors at positions " + allIndices(ids,:)); 
                                                 infoWord = tempErrorWord;
                                                 time = toc();
                                                 toc();
                                                 return                                                      
                                              end

                                          end 

                                      end 

                                  elseif err(probIdxER).errDel ~= 0 && err(probIdxER).errIns ~= 0 &&  err(probIdxER).errSub ~= 0 %111
                                      for ids = 1:length(allIndices)
										  tempErrorWord = char(repmat(0, 1, n));							  
                                          excludes = ismember(1:length(tempErrorWord),allIndices(ids,1:err(probIdxER).errDel));

                                          for uni = 1:length(uniqueDel)
											  tempErrorWord(allIndices(ids,1:err(probIdxER).errDel)) = char(uniqueDel(uni,:))'; 
                                              tempErrorWord(~excludes) = initWord(~excludes);
                                              combinations = calcSubComb(tempErrorWord(allIndices(ids,end-err(probIdxER).errSub + 1:end)));
											  tempErrorWord(allIndices(ids,err(probIdxER).errDel + 1:err(probIdxER).errDel + err(probIdxER).errIns)) = [];
                                              for uniSub = 1:length(combinations)
												  tempErrorWord(allIndices(ids,end - err(probIdxER).errSub + 1:end)) = combinations(uniSub,:);
												  [~,error] = CRCDetector(tempErrorWord,mapping,nucleotids);
											      if error ~= 0 
                                                  elseif error == 0 disp("We have found " + errorCount + " combo D + S + I errors at positions " + allIndices(ids,:)); 
                                                     infoWord = tempErrorWord;
                                                     time = toc();
                                                     toc();
                                                     return 
                                                  end
                                              end 
                                         end     
                                      end 

                                  end 

                            end    
                     end
                     err(probIdxER)=[];
                     end
              end    

         if error == 0 
            time = toc();
            return
         elseif error ~= 0
         decisionTree(probIdxDT)=[];        
         end
      end    
end
time = toc();
disp(toc);



