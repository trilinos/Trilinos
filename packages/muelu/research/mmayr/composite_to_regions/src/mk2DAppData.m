function [] = mk2DAppData(allMyNodes,allMyRegions,nProcs,whichCase,...
                       globalDims,localDims,relcorners,abscorners)
%
% Most of this information is computed in the mk2DRegionfile() and so
% here we basically need to just print it.
%

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

   nRegions = length(curRegionList);
   gDim   = ones(3,nRegions);
   for k=1:nRegions
      reg = curRegionList(k)+1;
      if localDims(reg,myRank+1,1) ~= -1, 
         fprintf(fp,'%d %d 1  %d %d 1  %d %d 0 %d %d 0\n',...
             globalDims(reg,1),globalDims(reg,2),...
             localDims(reg,myRank+1,1),localDims(reg,myRank+1,2),...
             relcorners(reg,myRank+1,1),relcorners(reg,myRank+1,2),...
             abscorners(reg,1),abscorners(reg,2));
       end
   end
   fclose(fp);
end;
