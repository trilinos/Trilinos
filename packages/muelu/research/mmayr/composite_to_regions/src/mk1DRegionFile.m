% mk1DRegionFile.m
%
% Assume a one-to-one mapping of regions to procs.
%
% Define a case by setting the vector nNodesPerRegion = [n1 n2 ... nN] with n1 ... nN
% being the number of nodes in region 1...N. The number of regions is given by
% this vector.
%
% Note: 
% - To allow for a coarsening rate of 3 with coarse nodes ending up on
%   region interfaces, n1...nN needs to be 3k+1
% - Number of nodes per region needs to be (k*3^L + 1) with k=1...N and L being
%   the numberOfLevels-1.
%

function [] = mk1DRegionFile(filename)

%% User-defined cases

if (strcmp(filename, 'caseTen') == true)
  nNodesPerRegion = [4 7 4];
elseif (strcmp(filename, 'caseEleven') == true)
  nNodesPerRegion = [7 4 7 4 4 7];
elseif (strcmp(filename, 'caseTwelve') == true)
  nNodesPerRegion = [28 37 25];
elseif (strcmp(filename, 'caseThirteen') == true)
%   nNodesPerRegion = [901 901 874 901 874 928 874 928 901];
%   nNodesPerRegion = [901 901 901 901 901 901 901 901 901 901 901 901];
%   nNodesPerRegion = [901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901 901];
  nNodesPerRegion = [901 901 874 901 874 928 874 928 847 901 901 874 982 874 928 874 928 901 766 901 874 901 874 928 874 928 901];
elseif (strcmp(filename, 'caseFourteen') == true)
  numLevels = 5;
  numRegions = 8;
  baseSizeFac = 15; 
  sizeDeviation = [0 1 1 -2 0 -2 2 0]; % 3 0 -1 0 1 1 -2 0 -2 2 0 3 0 -1 2];
  if (numRegions ~= length(sizeDeviation))
    error('Number of regions does not match length of "sizeDeviation".');
  end
  
  nNodesPerRegion = zeros(numRegions,1);
  for i = 1:numRegions
    nNodesPerRegion(i) = (baseSizeFac + sizeDeviation(i)) * 3^(numLevels-1) + 1;
  end
  disp(nNodesPerRegion');
else
  error('Unknown case "%s".', filename);
end

%% Process data

% test for coarsening rate
for i = 1:length(nNodesPerRegion)
  if mod(nNodesPerRegion(i),3) ~= 1
    error('Number of nodes in region %d is not 3k+1. Perfect coarsening by 4 is not possible.', i);
  end
end

% compute list of node GIDs, region, and proc info
data = [];
gidOffset = 0;
for i = 1:length(nNodesPerRegion)
  for j = 1:nNodesPerRegion(i)
%     if (j < nNodesPerRegion(i) || i == length(nNodesPerRegion)) % assign duplicated nodes to next proc
      currData = [gidOffset+j-i i-1 i-1];
%     else
%       currData = [gidOffset+j-i i-1 i];
%     end
    data = [data; currData];
  end
  gidOffset = sum(nNodesPerRegion(1:i));
end

% Compute summary data for file header
nNodes = max(data(:,1)) + 1;
nRegions = length(nNodesPerRegion);
nProcs = length(nNodesPerRegion);
% whichCase = 'RegionsSpanprocs';
whichCase = 'MultipleRegionsPerProc';


%% Print data to file
fp = fopen(filename,'w');
if fp ~= -1
  fprintf(fp,'     #nodes  #regions   #procs   #which case\n');
  fprintf(fp,'%8d %8d %8d       %s\n',nNodes,nRegions,nProcs,whichCase);
  fprintf(fp,'     nodeID   regionID   procID\n');
  
  for i = 1:size(data,1)
    fprintf(fp,'%8d  %8d %8d\n', data(i,1), data(i,2), data(i,3));
  end

end

end