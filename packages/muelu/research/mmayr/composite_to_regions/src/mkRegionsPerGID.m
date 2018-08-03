% mkRegionsPerGID.m
%
% Input:
%   myNodes         list of nodes owned by this processor
%   myRank          rank of this processor
%   maxRegPerGID    max number of regions that a single GID belongs to
%   outDir          path to output directory
%
% Outpu:
%   [none]
%
function [] = mkRegionsPerGID(myNodes, myRank, maxRegPerGID, outDir)

   fp = fopen(sprintf('%s/myRegionAssignment_%d', outDir, myRank), 'w');
   if fp == -1
      error('mkRegionsPerGID: cannot open myRegAssignment_%d\n',myRank);
   end
%    fprintf(fp,'LoadAndCommRegAssignments\n');
   myGIDs = getCompositeIDs(myNodes,myRank);
   count = 1;
   for i=1:length(myGIDs)
      while myNodes(count).ID ~= myGIDs(i), count = count+1; end;
      temp = -ones(maxRegPerGID,1);
      temp(1:length(myNodes(count).gRegions)) = myNodes(count).gRegions;
      fprintf(fp,'%d ',temp);
      fprintf(fp,'\n');
   end
   fclose(fp);
