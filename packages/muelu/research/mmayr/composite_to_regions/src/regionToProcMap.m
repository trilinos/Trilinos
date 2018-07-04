function [RegionToProc] = regionToProcMap(myNodes,myRank,filename)

%
% This function fills
%
%   RegionToProc.myRegions with region GIDs that I at least partially own
%
%   RegionToProc.startPtr 
%   RegionToProc.myProcs   so that myProcs(startPtr(k):startPtr(k+1)-1)
%                          gives all the procs that at least partially own a
%                          piece of myRegions(k).
%
% Normally, this calculation might be a bit involved in parallel with some
% communication or fancy Trilinos calls. In this case, we just compute it
% from the file.


%
% figure out the number of regions that I at least partially own
%

nMyRegions = 0;
nMyNodes = length(myNodes);
regionList= zeros(nMyNodes,1); % upper bound on # of local regions I might own

for i=1:nMyNodes,
   id      = myNodes(i).ID;
   regions = myNodes(i).gRegions;
   for j=1:length(regions)
      curRegion = regions(j);
      if myNodes(i).procIds(j) == myRank,
        notFound = 1;
        for k=1:nMyRegions,
          if  regionList(k) == curRegion, notFound=0; break; end; 
        end  % for k=1:nMyRegions,
        if notFound, nMyRegions = nMyRegions+1; regionList(nMyRegions)=curRegion;end;
      end; % if myNodes(i).procIds(j) == myRank,
   end;   % for j=1:length(regions)
end  % for i=1:nMyNodes
regionList = regionList(1:nMyRegions);

startPtr = zeros(nMyRegions+1,1);
procsLeng = 100*nMyRegions;
procs = -ones(procsLeng,1);
for i=1:nMyRegions, startPtr(i) = (i-1)*100+1; end;
regionList = sort(regionList);

% normally this part would require some nasty communication or high level
% Trilinos function to compute. Here, we just read it from the file and 
% figure things out in a slow clumsy way.

fp = fopen(filename,'r');
if (fp == -1)
  error('Cannot read %s \n',filename);
end

% read header information 
fgetl(fp);   % read heading line
temp = fscanf(fp,'%d');
nNodes   = temp(1);
nRegions = temp(2);
nProcs   = temp(3);
fgetl(fp);   % read rest of 2nd line
fgetl(fp);   % read heading line

newLine = 'dummy';
count   = 0;
oldNodeId = -1;

while newLine ~= -1,
   newLine = fgetl(fp);
   if newLine ~= -1,
      temp = sscanf(newLine,'%d');
      curRegion = temp(2);
      curProc   = temp(3);
      found = 0; 
      for i=1:nMyRegions, 
        if regionList(i) == curRegion, found=1; break; end;
      end;
      if found == 1,
         j=(i-1)*100+1;
         procFound = 0;
         while procs(j) ~= -1, 
            if procs(j) == curProc, procFound = 1; break; end;
            j = j+1;
         end
         if procFound ==0, 
            procs(j) = curProc;
            startPtr(i) = startPtr(i)+1; 
         end
      end; % if found == 1,
   end;  %  if newLine ~= -1,
end;  % while newLine ~= -1,
fclose(fp);

%
% now compress procs, which contains lots of -1's and properly set startPtr
%
startPtr(1)  = 1;
nextAvail    = 1;
nextToProcess= 1;

for i=1:nMyRegions,
   while procs(nextToProcess) ~= -1,
      procs(nextAvail) = procs(nextToProcess);
      nextToProcess= nextToProcess+1;
      nextAvail    = nextAvail+1;
   end

   startPtr(i+1)=nextAvail;
   nextToProcess = i*100+1;
end;
RegionToProc.myRegions = regionList;
RegionToProc.myProcs   = procs(1:nextAvail-1);
RegionToProc.startPtr  = startPtr;
