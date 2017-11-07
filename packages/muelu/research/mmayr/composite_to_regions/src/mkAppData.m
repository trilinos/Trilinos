function [] = mkAppData(allMyNodes,allMyRegions,nProcs)

for myRank=0:nProcs-1,
   if exist(sprintf('myData_%d',myRank),'file'),
      fprintf('mkAppData: myData_%d already exists\n',myRank);
      keyboard;
   end;
   fp = fopen(sprintf('myData_%d',myRank),'w');
   if fp == -1,
      fprintf('mkAppData: cannot open myData_%d\n',myRank);
      keyboard;
   end;
   fprintf(fp,'LoadAppDataForLIDregion()\n');

   curRegionList = allMyRegions{myRank+1}.myRegions;

   for k=1:length(curRegionList),
      maxGID = -1; minGID = 10000000;
      curRegion = curRegionList(k);

      % go through all GIDs on all procs looking for curRegion,

      for j=1:length(allMyNodes{myRank+1}),
         if ~isempty(find(curRegion == allMyNodes{myRank+1}(j).gRegions)),
            theID = allMyNodes{myRank+1}(j).ID;
            if theID < minGID, minGID = theID; end;
            if theID > maxGID, maxGID = theID; end;
         end
      end
      fprintf(fp,'%d %d\n',minGID,maxGID);
   end;
   fclose(fp);
end;
