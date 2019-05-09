The main files in this directory are


getCompositeIds :  returns GIDs of nodes that myRank owns in composite map.

mkRegionFile    :  writes a mapping file based on information that is hardwired
                   in the function. Thus, user's are expected to edit this
                   to define different regional meshes.

readRegionalFile:  fills structure myNodes with a list of all node owned or 
                   shared by myRank based on information in file. 

caseOne,caseTwo : in addition to obvious header information, each line of this
                   file consists of 3 numbers (nodeID,regionID,procID). There
                   should be one entry for each regional node. That is, nodes
                   that are shared between 2 regions should have 2 entries.
                   The entry that appears 1st for shared nodes corresponds to
                   the 'owner' in the composite mesh. It is assumed that
                   entries corresponding to shared ids are consecutive.


regionToProcMap :  fills a structure indicating all procs that partially own
                   each region that myRank partially owns.


waitForRmDataFiles : sleeps until there are no myData_* files in the directory

mkCompositeMap : writes myData_* files associated with composite map


send :          writes myData_* files that consist of a single command line
                passed in as an argument

mkRegionsPerGID : writes myData_* files associated with regionsPerGID

mkAppData :      writes myData_* files associated with min/max GID per region

caseOne         : describes a 1D problem with regions corresponding to the following

P 0000000000001111111111222222222222222222222222333333333 44444 55555666

  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23


  |--------------------<-----------------------<-----------------------|
          region 0             region 1               region 2

  region 0 crosses 2 procs     proc 0 has 1 region
  region 1 crosses 1 procs     proc 1 has 1 region
  region 2 crosses 4 procs     proc 2 has 1 region
                               proc 3 has 1 region
                               proc 4 has 1 region
                               proc 5 has 1 region
                               proc 6 has 1 region

so all procs do everything at the same time with MPI_COMM_WORLD

              row map                |    col map
P0  0,  1,  2,  3                    | 0, 1, 2, 3, 4
P1  4,  5,  6,  7                    | 3, 4, 5, 6, 7  [not 8]
P2  7,  8,  9, 10, 11, 12, 13, 14, 15| 7, 8, 9,10,11, 12, 13, 14, 15 [not 6,16]
P3 15, 16, 17, 18                    |15,16,17,18,19 [not 14]
P4 19, 20                            |18,19,20,21
P5 21, 22                            |20,21,22,23
P6 23                                |22,23

caseTwo         : describes a 1D problem with regions corresponding to the following

P 0000000000001111111111222222222222222222222222333333333444444444444444

  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23


  |-----<----->--------<----->-------->--------<---->--->----->--<-----|

     r0    r1     r2      r3     r4       r5     r6   r7   r8  r9   r10


  no regions cross procs        proc 0 has 2 regions
                                proc 1 has 1 region
                                proc 2 has 3 region
                                proc 3 has 2 regions
                                proc 4 has 3 regions

round 1 (r0, r2, r3, r6, r8)
              row map                | col map is always
P0  0,  1,  2                        | identical to row map
P1  4,  5,  6,  7                    | on all procs
P2  7,  8,  9                        |
P3 15, 16, 17                        |
P4 18, 19, 20                        |
=============================================================
round 2 (r1,  -, r4, r7, r9)
P0  2,  3,  4                        |
P1                                   |
P2  9, 10, 11, 12                    |
P3 17, 18                            |
P4 20, 21                            |
=============================================================
round 3 ( -,  -, r5,  -, r10)
P0                                   |
P1                                   |
P2 12, 13, 14, 15                    |
P3                                   |
P4 21, 22, 23                        |

