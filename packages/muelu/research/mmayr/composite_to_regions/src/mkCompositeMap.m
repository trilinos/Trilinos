function [] = mkCompositeMap(myNodes,myRank)

  % open file
  filename = 'myCompositeMap_';
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
