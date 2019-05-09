function [] = mk2DAppData(allMyNodes,allMyRegions,nProcs,whichCase,...
                       globalDims,localDims,relcorners,abscorners, outDir)
%
% Most of this information is computed in the mk2DRegionfile() and so
% here we basically need to just print it.
%

for myRank=0:nProcs-1
   % open file
   fp = fopen(sprintf('%s/myAppData_%d', outDir, myRank), 'w');
   if fp == -1
     error('mkAppData: cannot open myAppData_%d\n',myRank);
   end
   
   curRegionList = allMyRegions{myRank+1}.myRegions;

   nRegions = length(curRegionList);
   gDim   = ones(3,nRegions);
   for k=1:nRegions
     reg = curRegionList(k)+1;
     if localDims(reg,myRank+1,1) ~= -1
       fprintf(fp,'%d %d 1  %d %d 1  %d %d 0 %d %d 0\n',...
           globalDims(reg,1),globalDims(reg,2),...
           localDims(reg,myRank+1,1),localDims(reg,myRank+1,2),...
           relcorners(reg,myRank+1,1),relcorners(reg,myRank+1,2),...
           abscorners(reg,1),abscorners(reg,2));
     end
   end
   
   % close file
   fclose(fp);
end
