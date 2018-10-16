% readRegionFile.m
%
% Read region information from file.
%
% Input
%   filename  filename of case file to read
%   myRank    rank of this processor
%   inputDir  path to directory with region information
%
% Output
%   myNodes   list of nodes owned by this processor
%
function [myNodes] = readRegionFile(filename,myRank,inputDir)

%
% Construct the list myNodes of ndoes that I own or share for proc 'myRank'.
% Each element in the list has the following fields:
%
%   myNodes(k).ID           global ID of all nodes that I own or share
%
%   myNodes(k).gRegions     list of global region ids that ID belongs to 
%
%   myNodes(k).procIds      procIds(j) gives processor that owns gRegions(j)
%
% 
% IMPORTANT NOTE: myNodes(k).procIds(1), the 1st guy, gives the processor that
% owns a shared node within the composite mesh.  It is also assumed that
% entries corresponding to shared ids are consecutive.


fp = fopen(sprintf('%s/%s', inputDir, filename), 'r');
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

%  read each line of the file and only retain those things that I own
%  this includes recording all other procs and region ids correspond
%  to interface nodes that I own
newLine = 'dummy';
iOwnIt = 0;
myNodes = [];
count   = 0;
oldNodeId = -1;

while newLine ~= -1,
   newLine = fgetl(fp);
   if newLine ~= -1,
      temp = sscanf(newLine,'%d');
      currentNodeId   = temp(1);
      currentRegionId = temp(2);
      currentProcId   = temp(3);
      % store previous node if I owned it.
      if currentNodeId ~= oldNodeId;   % working on a new node so that we need
                                       % to store previous node if I owned it.
         if iOwnIt,
             count = count+1;
             myNodes(count).ID      = oldNodeId;
             myNodes(count).gRegions= globalRegionList;
             myNodes(count).procIds = procList;
         end
         % For the new node, record that I don't yet know if I own it
         % and haven't created the associated lists yet.
         iOwnIt = 0; globalRegionList=[]; procList=[];
      end

      % now work on current node
      if (currentProcId == myRank)  % order so that composite owning proc is 1st
                                    % in procList
         iOwnIt = 1;
         if isempty(procList) % first entry for this node indicating owning
                              % proc in composite list
            globalRegionList = [currentRegionId globalRegionList];
            procList         = [currentProcId         procList  ];
         else
            globalRegionList = [globalRegionList currentRegionId];
            procList         = [procList           currentProcId];
         end
      else
         globalRegionList = [globalRegionList currentRegionId];
         procList         = [procList           currentProcId];
      end
      oldNodeId = currentNodeId;
   else
      % last node in list which must be stored if I own it.
         if iOwnIt,
             count = count+1;
             myNodes(count).ID = oldNodeId;
             myNodes(count).gRegions = globalRegionList;
             myNodes(count).procIds  = procList;
         end
   end
end
fclose(fp);
