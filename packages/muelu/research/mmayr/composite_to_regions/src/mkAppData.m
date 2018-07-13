function [] = mkAppData(allMyNodes,allMyRegions,nProcs,whichCase)

for myRank=0:nProcs-1
   
   % open file
   fp = fopen(sprintf('myAppData_%d',myRank),'w');
   if fp == -1
     error('mkAppData: cannot open myAppData_%d\n',myRank);
   end
%    fprintf(fp,'LoadAppDataForLIDregion()\n');

   curRegionList = allMyRegions{myRank+1}.myRegions;

   nRegions = length(curRegionList);
   gDim   = ones(3,nRegions);
   for k=1:nRegions
     maxGID = -1; minGID = 10000000;
     curRegion = curRegionList(k);

     % go through all GIDs on all procs looking for curRegion,

     for j=1:length(allMyNodes{myRank+1})
       if ~isempty(find(curRegion == allMyNodes{myRank+1}(j).gRegions))
         theID = allMyNodes{myRank+1}(j).ID;
         if theID < minGID, minGID = theID; end
         if theID > maxGID, maxGID = theID; end
       end
     end
     gDim(1,k) = maxGID - minGID + 1;
     fprintf(fp,'%d %d\n',minGID,maxGID);
   end
      
   lDim   = ones(3,nRegions);
   lowInd = zeros(3,nRegions);
   if whichCase(1) == 'M'
     lDim = gDim;
   elseif whichCase(1) == 'R'
     temp = allMyNodes{myRank+1};
     maxi = 0; mini = temp(1).ID;
     lowerCorner = mini;
     for i=1:length(allMyNodes)
       temp = allMyNodes{i};
       for j=1:length(temp)
         if find(temp(j).gRegions == allMyRegions{myRank+1}.myRegions)
           if temp(j).ID > maxi, maxi = temp(j).ID; end;
           if temp(j).ID < mini, mini = temp(j).ID; end;
           if (i==myRank+1) && (temp(j).ID < lowerCorner), lowerCorner = temp(j).ID; end;
         end
       end
     end
     gDim(1,1) = maxi-mini+1;
     lDim(1,1) = length(allMyNodes{myRank+1});
     lowInd(1,1) = lowerCorner-mini;
   else
     fprintf('whichCase not right\n'); keyboard;
   end
   for ii=1:nRegions
     fprintf(fp,'%d %d %d    %d %d %d   %d %d %d\n',gDim(1,ii),gDim(2,ii),gDim(3,ii),...
                                                 lDim(1,ii),lDim(2,ii),lDim(3,ii),...
                                                 lowInd(1,ii),lowInd(2,ii),lowInd(3,ii));
   end
   
   % close file
   fclose(fp);
end
