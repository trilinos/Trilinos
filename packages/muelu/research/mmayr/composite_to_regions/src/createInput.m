% createInput.m
%
% Create input files with composite matrix and region information and write to
% disk. This will be used as starting point for the c++-based driver.
%

!rm -f compX.mm map_compX.mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% USER-PROVIDED INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDims = 1; file = 'caseTen';
% nDims = 2; file = 'caseTwenty';

outDir = sprintf('%s/regDataDir', pwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup the output directory
eval(sprintf('!mkdir -p %s', outDir));
eval(sprintf('!rm -f %s/*', outDir));
eval(sprintf('!cp case* %s', outDir));

if (nDims == 2)
  [globalDims,localDims,relcorners,abscorners] = mk2DRegionFile(file, outDir);
elseif (nDims == 1) && ...
    ((strcmp(file,'caseOne') == false) && (strcmp(file,'caseTwo')) == false && ...
    (strcmp(file,'caseThree') == false) && (strcmp(file,'caseFour') == false))
  
  mk1DRegionFile(file, outDir);
end

%
%  read in some of header type information 
fp = fopen(sprintf('%s/%s', outDir, file), 'r'); fgetl(fp); 
t = fscanf(fp,'%d'); nNodes = t(1); nProcs = t(3); 
whichCase = fgets(fp,16); fclose(fp);

%
% make a global tridiagonal matrix corresponding to the discretization
%
% A = spdiags(rand(nNodes,3), -1:1, nNodes, nNodes);
% A = spdiags(-1*ones(nNodes,3), -1:1, nNodes, nNodes);
if (nDims == 1)
  A = oneDimensionalLaplace(nNodes);
elseif (nDims == 2)
  A = twoDimensionalLaplace(nNodes);
else
  error('nDims is wrong\n');
end

% Each 'processor' reads regional file. 
for myRank=0:nProcs-1
eval(sprintf('myNodes%d = readRegionFile(file,myRank,outDir);',myRank));
eval(sprintf('RegionToProc%d=regionToProcMap(myNodes%d,myRank,file,outDir);',...
               myRank,myRank));
end

maxRegPerProc = 0;
maxRegPerGID = 0;
for myRank=0:nProcs-1
  eval(sprintf('allMyNodes{%d} = myNodes%d;',myRank+1,myRank));
  eval(sprintf('allMyRegions{%d} = RegionToProc%d;',myRank+1,myRank));
  if length(allMyRegions{myRank+1}.myRegions) > maxRegPerProc 
    maxRegPerProc = length(allMyRegions{myRank+1}.myRegions);
  end
  for k=1:length(allMyNodes{myRank+1})
    if length(allMyNodes{myRank+1}(k).gRegions) > maxRegPerGID 
      maxRegPerGID = length(allMyNodes{myRank+1}(k).gRegions);
    end
  end
end

% send header to C++
for myRank=0:nProcs-1
  fp=fopen(sprintf('%s/myRegionInfo_%d',outDir,myRank),'w');
  if (nDims == 1)
    fprintf(fp,'%d %d %d %d %s\n',maxRegPerGID,maxRegPerProc, nNodes, 1, whichCase);
  elseif (nDims == 2)
    fprintf(fp,'%d %d %d %d %s\n',maxRegPerGID,maxRegPerProc,sqrt(nNodes),sqrt(nNodes),whichCase);
  else
    error("Problems of spatial dimension %d are not implemented, yet.", nDims);
  end
  fclose(fp); 
end

% make and send composite map to C++
for myRank=0:nProcs-1
  eval(sprintf('mkCompositeMap(myNodes%d,%d,outDir);',myRank,myRank));
end

% write matrix to file using MatrixMarket format
mkMatrixFile(A, outDir);

for myRank=0:nProcs-1
  eval(sprintf('mkRegionsPerGID(myNodes%d,%d,maxRegPerGID,outDir);',myRank,myRank));
end

for myRank=0:nProcs-1
  fp=fopen(sprintf('%s/myRegions_%d',outDir,myRank),'w');
  for j=1:length(allMyRegions{myRank+1}.myRegions)
    fprintf(fp,'%d\n', allMyRegions{myRank+1}.myRegions(j));
  end
  fclose(fp); 
end

if (nDims == 1)
  mkAppData(allMyNodes,allMyRegions,nProcs,whichCase,outDir);
elseif (nDims == 2)
  mk2DAppData(allMyNodes,allMyRegions,nProcs,whichCase,globalDims,localDims,relcorners,abscorners,outDir);
else
  error("Problems of spatial dimension %d are not implemented, yet.", nDims);
end

str = sprintf('Run now with #procs = %d', nProcs);
disp(str);

