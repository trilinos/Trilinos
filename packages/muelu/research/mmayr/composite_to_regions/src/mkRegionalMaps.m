function [] = mkRegionalMaps(myNodes,myRegions,myRank)

% create file with the following format for each region (but do this serially
% waiting for the file to be destroyed by someone else before creating the
% next file).
%
% MakeSubCommAndRegionalMap
%   Region GID
%   #procs in sub-communicator
%   1st Proc in sub-communicator
%   2nd Proc in sub-communicator
%   ...
%   1st GID that k owns within this region
%   2nd GID that k owns within this region
%   ...
%
   pause on;
   if exist(sprintf('myData_%d',myRank),'file'),
      fprintf('mkRegionalMaps: myData_%d already exists\n',myRank);
      keyboard;
   end;
   for i=1:length(myRegions.myRegions)
      fp = fopen(sprintf('myData_%d',myRank),'w');
      if fp == -1,
         fprintf('mkRegionalMaps: cannot open myData_%d\n',myRank);
         keyboard;
      end;
      fprintf(fp,'MakeSubCommAndRegionalMap\n');
      curRegion = myRegions.myRegions(i);
      fprintf(fp,'%d\n',curRegion);
      nProcs = myRegions.startPtr(i+1) - myRegions.startPtr(i);
      fprintf(fp,'%d\n',nProcs);
      temp = sort(myRegions.myProcs(myRegions.startPtr(i):myRegions.startPtr(i+1)-1));
      for j=0:nProcs-1,
         fprintf(fp,'%d\n',temp(j+1));
      end;
%      for j=myRegions.startPtr(i):myRegions.startPtr(i+1)-1
%         fprintf(fp,'%d\n',myRegions.myProcs(j));
%      end;
      for j=1:length(myNodes)
         z = find(myNodes(j).gRegions == curRegion);
         if ~isempty(z) && (myNodes(j).procIds(z)==myRank),
            fprintf(fp,'%d\n',myNodes(j).ID);
         end;
      end;
      fclose(fp);
fprintf('changing pause here ..... \n');
pause(1.);
%      while exist(sprintf('myData_%d',myRank),'file') ~= 0, 
%         pause(.5);
%      end;
   end;
             
   pause off;


