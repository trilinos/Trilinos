function [] = mkRegionsPerGID(myNodes,myRank,maxRegPerGID)

   if exist(sprintf('myData_%d',myRank),'file'),
      fprintf('mkRegionsPerGID: myData_%d already exists\n',myRank);
      keyboard;
   end;
   fp = fopen(sprintf('myData_%d',myRank),'w');
   if fp == -1,
      fprintf('mkRegionsPerGID: cannot open myData_%d\n',myRank);
      keyboard;
   end;
   fprintf(fp,'LoadAndCommRegAssignments\n');
   myGIDs = getCompositeIDs(myNodes,myRank);
   count = 1;
   for i=1:length(myGIDs)
      while myNodes(count).ID ~= myGIDs(i), count = count+1; end;
      temp = -ones(maxRegPerGID,1);
      temp(1:length(myNodes(count).gRegions)) = myNodes(count).gRegions;
      fprintf(fp,'%d ',temp);
      fprintf(fp,'\n');
   end;
   fclose(fp);
