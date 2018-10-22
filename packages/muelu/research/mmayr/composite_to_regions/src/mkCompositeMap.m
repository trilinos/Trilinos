% mkCompositeMap.m
%
% Write composite map to file 'myCompositeMap_*' on each processor.
%
% Input:
%   myNodes   list of all nodes owned by this processor
%   myRank    rank of this processor
%   outDir    path to output directory
%
% Output: 
%   [none]
%
function [] = mkCompositeMap(myNodes,myRank,outDir)

  % open file
  filename = sprintf('%s/myCompositeMap_', outDir);
  fp = fopen(sprintf('%s%d',filename,myRank),'w');
  if fp == -1
    error('mkCompositeMap: cannot open myData_%d\n',myRank);
  end
  
  % write to file
%   fprintf(fp,'LoadCompositeMap\n');
  myGIDs = getCompositeIDs(myNodes,myRank);
  for i=1:length(myGIDs)
    fprintf(fp,'%d\n',myGIDs(i));
  end
  
  % close file
  fclose(fp);
  
end
