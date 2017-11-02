%
% Because I'm lazy, I'm going to write files in matlab. These files will be
% used by a Trilinos driver running concurrent to my matlab stuff. The
% Trilinos driver on proc k will have a big loop that waits while for a file 
% called myData_k. The first line of this file will have a command that 
% determines which part of a case/switch statement is executed within the loop. 
% The remainder of the file may have more data relevant to the command.
% When the proc is finished with myData_k, it will remove this file. This
% will be a signal to the matlab program that it can proceed to generate
% a new myData_k file with a different command. 
% 
% Here are some possibilities for myData_k:
%
%
% ------------------------------------------------------------------------
% The very first set of myData_k files contain just one line with two
% numbers and then a string (each separated by a single space!).
%
% The first number is the maxRegionsPerGID
% The second number is the maxRegionsPerProc
% The string must be either
%
%      RegionsSpanProcs
%             or
%      MultipleRegionsP
%
% indicating which of the two cases we are in.
% ------------------------------------------------------------------------
% LoadCompositeMap
%   1st GID assigned to k
%   2nd GID assigned to k
%      ...
% ------------------------------------------------------------------------
% LoadAndCommRegAssignments
% 1stRegionSharedByMy1stCompositeNode 2ndRegionSharedByMy1stCompositeNode ...
% 1stRegionSharedByMy2ndCompositeNode 2ndRegionSharedByMy2ndCompositeNode ...
%      .
%      .
%      .
% Here, the number of columns should be equal to maxRegionsPerGID and -1's
% would appear in column j if the composite node associated with this
% line in the file is not shared by j regions.
% ------------------------------------------------------------------------
% LoadRegions
% My1stRegion
% My2ndRegion
%      .
%      .
%      .
% ------------------------------------------------------------------------
% LoadAppDataForLIDregion()
% My1stRegionsMinGID  My1stRegionsMaxGID 
% My2ndRegionsMinGID  My2ndRegionsMaxGID 
%      .
%      .
%      .
% ------------------------------------------------------------------------
% ReadMatrix 
% ------------------------------------------------------------------------
% MakeGrpRegRowMaps
% ------------------------------------------------------------------------
% MakeGrpRegColMaps
% ------------------------------------------------------------------------
% MakeExtendedGrpRegMaps
% ------------------------------------------------------------------------
% MakeRegionMatrices
% ------------------------------------------------------------------------
% PrintCompositeMap
% ------------------------------------------------------------------------
% PrintRegionAssignments
% ------------------------------------------------------------------------
% PrintGrpRegColMaps
% ------------------------------------------------------------------------
% PrintRevisedColMaps
% ------------------------------------------------------------------------
% PrintGrpRegDomMaps
% ------------------------------------------------------------------------
% PrintGrpRegRowMaps
% ------------------------------------------------------------------------
% PrintRegionMatrices
% ------------------------------------------------------------------------
% Terminate
% ------------------------------------------------------------------------
%
%
%
% Read in region information so that we can setup a 1D Poisson problem
% on a composite mesh as well as on each of the individual region meshes
%
!rm -f myData_* 

file='caseOne'; % mkRegionFile(file);

%
%  read in some of header type information 
fp = fopen(file,'r'); fgetl(fp); 
t = fscanf(fp,'%d'); nNodes = t(1); nProcs = t(3); 
whichCase = fgets(fp,16); fclose(fp);

%
% make a global tridiagonal matrix corresponding to the discretization
%
A = spdiags(rand(nNodes,3), -1:1, nNodes, nNodes);
% A = spdiags(ones(nNodes,3), -1:1, nNodes, nNodes);


% Each 'processor' reads regional file. 

for myRank=0:nProcs-1
   eval(sprintf('myNodes%d = readRegionalFile(file,myRank);',myRank));
   eval(sprintf('RegionToProc%d=regionToProcMap(myNodes%d,myRank,file);',...
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

%
% send header to C++
%
for myRank=0:nProcs-1
   fp=fopen(sprintf('myData_%d',myRank),'w');
   fprintf(fp,'%d %d %s\n',maxRegPerGID,maxRegPerProc,whichCase);
   fclose(fp); 
end
waitForRmDataFiles(nProcs);

%
% make and send composite map to C++
%
for myRank=0:nProcs-1
   eval(sprintf('mkCompositeMap(myNodes%d,%d,maxRegPerGID);',myRank,myRank));
end
waitForRmDataFiles(nProcs);
%send('PrintCompositeMap',nProcs);

%
% now send matrix to C++
%
fp = fopen('Amat.mm','w');
fprintf(fp,'%%%%MatrixMarket matrix coordinate real general\n');
fprintf(fp,'%d %d %d\n',size(A,1),size(A,2),nnz(A));
[aaa,bbb,ccc] = find(A);
[dummy,ii] = sort(aaa);
aaa = aaa(ii); bbb = bbb(ii); ccc = ccc(ii);
row = aaa(1);
startRow = 1;  endRow = startRow;
for i=2:length(aaa)
   if aaa(i) ~= row % sort and print out the previous row
      aa = aaa(startRow:i-1);
      bb = bbb(startRow:i-1);
      cc = ccc(startRow:i-1);
      [dummy,ii] = sort(bb);
      aa = aa(ii);
      bb = bb(ii);
      cc = cc(ii);
      for k=1:length(aa)
        fprintf(fp,'%d %d %20.13e\n',aa(k),bb(k),cc(k));
      end
      startRow=i;   % indicates that we are now storing a new row
   end
end
% print the last row
aa = aaa(startRow:end);
bb = bbb(startRow:end);
cc = ccc(startRow:end);
[dummy,ii] = sort(bb);
aa = aa(ii);
bb = bb(ii);
cc = cc(ii);
for k=1:length(aa)
  fprintf(fp,'%d %d %20.13e\n',aa(k),bb(k),cc(k));
end
fclose(fp);

% now signal that the real program that the matrix is ready
send('ReadMatrix',nProcs);

%
% make and send to C++ regionsPerGID info
%
waitForRmDataFiles(nProcs);
for myRank=0:nProcs-1
   eval(sprintf('mkRegionsPerGID(myNodes%d,%d,maxRegPerGID);',myRank,myRank));
end
%send('PrintRegionAssignments',nProcs);

waitForRmDataFiles(nProcs);
for myRank=0:nProcs-1
   fp=fopen(sprintf('myData_%d',myRank),'w');
   fprintf(fp,'LoadRegions\n');
   for j=1:length(allMyRegions{myRank+1}.myRegions)
      fprintf(fp,'%d\n', allMyRegions{myRank+1}.myRegions(j));
   end
   fclose(fp); 
end

waitForRmDataFiles(nProcs);
mkAppData(allMyNodes,allMyRegions,nProcs);
waitForRmDataFiles(nProcs);

send('MakeGrpRegRowMaps',nProcs);
send('MakeGrpRegColMaps',nProcs);
send('MakeQuasiRegionMatrices',nProcs);
% send('PrintQuasiRegionMatrices',nProcs);
send('MakeExtendedGrpRegMaps',nProcs);
% send('PrintGrpRegDomMaps',nProcs);
% send('PrintRevisedRowMaps',nProcs);
% send('PrintRevisedColMaps',nProcs);
% send('PrintGrpRegColMaps',nProcs);
send('MakeRegionMatrices',nProcs);
% send('PrintRegionMatrices',nProcs);
% send('PrintRegionMatrixRowMap',nProcs);
% send('PrintRegionMatrixColMap',nProcs);
% send('PrintRegionMatrixRangeMap',nProcs);
% send('PrintRegionMatrixDomainMap',nProcs);
send('ComputeMatVecs',nProcs);
% send('PrintCompositeVectorX',nProcs);
% send('PrintCompositeVectorY',nProcs);
% send('PrintQuasiRegVectorX',nProcs);
% send('PrintRegVectorX',nProcs);
% send('PrintRegVectorY',nProcs);
% send('PrintRegVectorYComp',nProcs);
send('Terminate',nProcs);

