function mkRegionFile(filename)
%
%  writes mapping to a file based on the information in the 
%  'edit this part' section of the code below. In the region file,
%  it is assumed that the first entry for any shared node corresponds
%  to the owner in the corresponding composite mesh.
%
%  Note: regionStart(nRegions+1) should be equal to the total number of nodes
%        in the mesh. More generally, the kth region starts at regionStart(k)
%        and includes all nodes up to regionStart(k+1)
%
%  Further Note:  procAssignment generally gives the assignment of nodes to
%                 processors in any regional sense. Of course, different
%                 processors might own the same node if it belongs to two
%                 regions and each region is associated with a different
%                 processor. As there is only one procAssignment entry for 
%                 each node, we need some additional information to handle
%                 this case. Specifically, procAssignment(k) gives the 
%                 processor assignment of any kth point that corresponds to the 
%                 rightmost boundary of a region. However, the processor 
%                 assignment of any kth point that is at a leftmost region 
%                 boundary is determined by procAssignment(k+1).
%                 THIS MEANS THAT procAssignment(k) should be set to RegionA when k corresponds to the 
%                 rightmost point of regionA and the leftmost point of regionB.
%                 
% Even Further Note: 'readRegionalFile.m' assigns the 1st shared node
% encountered to the composite map. As this generator always orders 
% first nodes within a region that is to the left of another region,
% this implies that shared nodes are always associated with the 
% the region to the left in the composite map. It might be possible to 
% change this (I haven't checked) by editing the file produced by this
% function so that two interface nodes appear in reverse order.
% ---------------------------------------------------
% ------------   edit this part  --------------------
% ---------------------------------------------------

example = 1;

if example == 1,
% case 1 : regions may cross processors but no processor owns a piece of more
%          than one region
% 
%
%P 0000000000001111111111222222222222222222222222333333333 44444 55566666
%
%  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
%
%
%  |--------------------<-----------------------<-----------------------|
%          region 0             region 1               region 2
%
%  region 0 crosses 2 procs     proc 0 has 1 region
%  region 1 crosses 1 procs     proc 1 has 1 region
%  region 2 crosses 4 procs     proc 2 has 1 region
%                               proc 3 has 1 region
%                               proc 4 has 1 region
%                               proc 5 has 1 region
%                               proc 6 has 1 region
%
   whichCase = 'RegionsSpanProcs';
   nRegions = 3;
   regionStart(1) = 0;   
   regionStart(2) = 7;   
   regionStart(3) = 15;  
   regionStart(4) = 23;  
   procAssignment = [ 0  0  0  0  1  1  1  1 ...
                      2  2  2  2  2  2  2  2 ...
                      3  3  3  4  4  5  5  6];
end;
if example == 2,
% case 2 : no regions cross processors but processors may own several regions 
%
%P 0000000000001111111111222222222222222222222222333333333344444444444444
%
%  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
%
%
%  |-----<-----<--------<-----<--------<--------<----<---<-----<--<-----|
%
%     r0    r1     r2      r3     r4       r5     r6   r7   r8  r9   r10
%
%
%  no regions cross procs        proc 0 has 2 regions
%                                proc 1 has 1 region
%                                proc 2 has 3 region
%                                proc 3 has 2 regions
%                                proc 4 has 3 regions
%
%round 1 (r0, r2, r3, r6, r8)
%              row map                | col map is always
%P0  0,  1,  2                        | identical to row map
%P1  4,  5,  6,  7                    | on all procs
%P2  7,  8,  9                        |
%P3 15, 16, 17                        |
%P4 19, 20                            |
%=============================================================
%round 2 (r1,  -, r4, r7, r9)
%P0  2,  3,  4                        |
%P1                                   |
%P2  9, 10, 11, 12                    |
%P3 17, 18                            |
%P4 20, 21                            |
%=============================================================
%round 3 ( -,  -, r5,  -, r10)
%P0                                   |
%P1                                   |
%P2 12, 13, 14, 15                    |
%P3                                   |
%P4 21, 22, 23                        |
   whichCase = 'MultipleRegionsPerProc';
   nRegions = 11;
   regionStart(1) = 0;   
   regionStart(2) = 2;   
   regionStart(3) = 4;  
   regionStart(4) = 7;  
   regionStart(5) = 9;  
   regionStart(6) =12;  
   regionStart(7) =15;  
   regionStart(8) =17;  
   regionStart(9) =18;  
   regionStart(10)=20;  
   regionStart(11)=21;  
   regionStart(12)=23;  
   procAssignment = [ 0  0  0  0  0  1  1  1 ...
                      2  2  2  2  2  2  2  2 ...
                      3  3  3  4  4  4  4  4];
end

% ---------------------------------------------------
% ------------   end of edit this part  -------------
% ---------------------------------------------------

nNodes = regionStart(nRegions+1)+1;
if length(procAssignment) ~= nNodes, fprintf('proc Assignment is wrong length\n');keyboard;end;

if (whichCase ~= 'RegionsSpanProcs') & (whichCase ~= 'MultipleRegionsPerProc'),
   fprintf('which Case not set properly\n'); keyboard;
end;
procsRegionAssignment = [];
regionsProcAssignment = [];
if whichCase == 'RegionsSpanProcs', procsRegionAssignment = -ones(max(procAssignment)+1,1); end;
if whichCase == 'MultipleRegionsPerProc', regionsProcAssignment = -ones(nRegions,1); end;
fp = fopen(filename,'w');
if fp ~= -1, 
  fprintf(fp,'     #nodes  #regions   #procs   #which case\n');
  fprintf(fp,'%8d %8d %8d       %s\n',nNodes,nRegions,max(procAssignment)+1,whichCase);
  fprintf(fp,'     nodeID   regionID   procID\n');
  priorRegion = -1;
  for i=0:nRegions-1,
     % leftmost node can be assigned to procAssignment(j+1   or j+2)
     j = regionStart(i+1); 
     pAssign = procAssignment(j+1);
     if i ~=0, pAssign= procAssignment(j+2); end;
     [procsRegionAssignment,regionsProcAssignment] = checkAssumptions(whichCase,procsRegionAssignment,regionsProcAssignment,pAssign,i);
     fprintf(fp,'%8d  %8d  %8d\n',j,i,pAssign);
     for j=regionStart(i+1)+1:regionStart(i+2)
        [procsRegionAssignment,regionsProcAssignment] = checkAssumptions(whichCase,procsRegionAssignment,regionsProcAssignment,procAssignment(j+1),i);
        fprintf(fp,'%8d  %8d  %8d\n',j,i,procAssignment(j+1));
     end;
  end;
  fclose(fp);
else
   fprintf('could not open file %s\n',filename);
   keyboard;
end;

% an old file description that is no longer allowed as a valid use case
% nRegions = 3;
% regionStart(1) = 0;    
% regionStart(2) = 7;    
% regionStart(3) = 15;   
% regionStart(4) = 23;   
% procAssignment = [ 0  0  0  0  1  1  1  1 ...
%                    2  2  1  1  1  1  1  1 ...
%                    3  3  3  4  4  5  5  6];
