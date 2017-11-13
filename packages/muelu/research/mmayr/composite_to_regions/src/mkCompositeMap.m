function [] = mkCompositeMap(myNodes,myRank,maxRegPerGID)

   if exist(sprintf('myData_%d',myRank),'file'),
      fprintf('mkCompositeMap: myData_%d already exists\n',myRank);
      keyboard;
   end;
   fp = fopen(sprintf('myData_%d',myRank),'w');
   if fp == -1,
      fprintf('mkCompositeMap: cannot open myData_%d\n',myRank);
      keyboard;
   end;
   fprintf(fp,'LoadCompositeMap\n');
   myGIDs = getCompositeIDs(myNodes,myRank);
   count = 1;
   for i=1:length(myGIDs)
      fprintf(fp,'%d\n',myGIDs(i));
   end;
   fclose(fp);
