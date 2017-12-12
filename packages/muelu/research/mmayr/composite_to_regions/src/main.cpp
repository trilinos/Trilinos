#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <iostream>

#define RegionsSpanProcs  1
#define MultipleRegionsPerProc  2
#include "ml_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_MatrixMatrix.h"

#include "Teuchos_Assert.hpp"

int LIDregion(void *ptr, int LIDcomp, int whichGrp);

// this little widget handles application specific data
// used to implement LIDregion()
struct widget {
   int *minGIDComp;
   int *maxGIDComp;
   int *myRegions;
   Epetra_Map *colMap;
   int maxRegPerGID;
   Epetra_MultiVector *regionsPerGIDWithGhosts;
};

void stripTrailingJunk(char *command);
void printGrpMaps(std::vector < Epetra_Map* > &mapPerGrp, int maxRegPerProc, char *str);

void compositeToRegional(Epetra_Vector* compVec,
    std::vector<Epetra_Vector*>& quasiRegVecs,
    std::vector<Epetra_Vector*>& regVecs, const int maxRegPerProc,
    std::vector<Epetra_Map*> rowMapPerGrp,
    std::vector<Epetra_Map*> revisedRowMapPerGrp,
    std::vector<Epetra_Import*> rowImportPerGrp)
{
  // quasiRegional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create empty vectors and fill it by extracting data from composite vector
    quasiRegVecs[j] = new Epetra_Vector(*(rowMapPerGrp[j]), true);
    int err = quasiRegVecs[j]->Import(*compVec, *(rowImportPerGrp[j]), Insert);
    TEUCHOS_ASSERT(err == 0);
  }

  // regional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create regVecs vector (copy from quasiRegVecs and swap the map)
    regVecs[j] = new Epetra_Vector(*(quasiRegVecs[j]));
    int err = regVecs[j]->ReplaceMap(*(revisedRowMapPerGrp[j]));
    TEUCHOS_ASSERT(err == 0);
  }

  return;
}

void createRegionalVector(std::vector<Epetra_Vector*>& regVecs,
    const int maxRegPerProc, std::vector<Epetra_Map*> revisedRowMapPerGrp)
{
  for (int j = 0; j < maxRegPerProc; j++)
    regVecs[j] = new Epetra_Vector(*(revisedRowMapPerGrp[j]), true);
  return;
}

// Interact with a Matlab program through a bunch of myData_procID files.
// Basically, the matlab program issues a bunch of commands associated with
// either composite or regional things (via the files). This program reads
// those files and performs the commands. See unpublished latex document
// ``Using Trilinos Capabilities to Facilitate Composite-Regional Transitions''

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  MPI_Group world_group;
  MPI_Comm_group(Comm.Comm(),&world_group);
#else
  Epetra_SerialComm Comm;
#endif
  int  myRank;
  FILE *fp;
  char fileName[20];
  char command[40];

  myRank = Comm.MyPID();

  // read maxRegPerGID (maximum # of regions shared by any one node) and
  //      maxRegPerProc (maximum # of partial regions owned by any proc)
  //      whichCase (MultipleRegionsPerProc: No regions spans across multiple
  //                                procs but a proc might own multiple regions
  //                 RegionsSpanProcs: processors own only a piece of 1 region
  //                                but regions may span across many procs

  int maxRegPerGID, maxRegPerProc, whichCase, iii = 0;
  sprintf(fileName,"myData_%d",myRank);
  while ( (fp = fopen(fileName,"r") ) == NULL) sleep(1);
  fgets(command,80,fp);
  sscanf(command,"%d",&maxRegPerGID);
  while ( command[iii] != ' ') iii++;
  sscanf(&(command[iii+1]),"%d",&maxRegPerProc);
  while ( command[iii] == ' ') iii++;
  iii++;
  while ( command[iii] != ' ') iii++;
  if      (command[iii+1] == 'M') whichCase = MultipleRegionsPerProc;
  else if (command[iii+1] == 'R') whichCase = RegionsSpanProcs;
  else {fprintf(stderr,"%d: head messed up %s\n",myRank,command); exit(1);}
  sprintf(command,"/bin/rm -f myData_%d",myRank); system(command); sleep(1);

  // ******************************************************************
  // Application Specific Data for LIDregion()
  // ******************************************************************
  struct widget appData;                            // ****************
  std::vector<int>  minGIDComp(maxRegPerProc);      // ****************
  std::vector<int>  maxGIDComp(maxRegPerProc);      // ****************
  // ******************************************************************
  // ******************************************************************


  std::vector<int> myRegions; // regions that myRank owns
  Epetra_CrsMatrix* AComp = NULL; // composite form of matrix
  Epetra_CrsMatrix* ACompSplit = NULL; // composite form of matrix
  Epetra_Map* mapComp= NULL; // composite map used to build AComp

  // regionsPerGIDWithGhosts[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix col Map.
  Epetra_MultiVector *regionsPerGIDWithGhosts = NULL;

  /* rowMapPerGrp[i] and colMapPerGrp[i] are based on the composite maps, i.e.
   * they don't include duplication of interface DOFs.
   *
   * revisedRowMapPerGrp[i] and revisedColMapPerGrp[i] are build manually to
   * account for duplicated, but distinct interface nodes, i.e. they inlucde
   * new GIDs for the duplicated nodes.
   *
   * We make some assumptions:
   * - For the composite as well as the region maps , we assume that
   *   row, column, and range map to be the same. So, we only need a
   *   row and a column map. Range and domain map can be taken as the row map.
   * - For the quasiRegional maps, we only need row and column maps. FillComplete()
   *   will generate some range and domain maps, though we will never use them
   *   since the quasiRegional matrices are only an intermediate step and are
   *   never used for actual calculations.
   */

  std::vector < Epetra_Map* > rowMapPerGrp(maxRegPerProc); // row map associated with myRank's ith region in composite layout
  std::vector < Epetra_Map* > colMapPerGrp(maxRegPerProc); // column map associated with myRank's ith region in composite layout
  std::vector < Epetra_Map* > revisedRowMapPerGrp(maxRegPerProc); // revised row map associated with myRank's ith region for regional layout
  std::vector < Epetra_Map* > revisedColMapPerGrp(maxRegPerProc); // revised column map associated with myRank's ith region for regional layout

  std::vector< Epetra_Import* > rowImportPerGrp(maxRegPerProc); // row importers per group
  std::vector< Epetra_Import* > colImportPerGrp(maxRegPerProc); // column importers per group
  std::vector< Epetra_Export* > rowExportPerGrp(maxRegPerProc); // row exporters per group

  std::vector< Epetra_CrsMatrix * > quasiRegionGrpMats(maxRegPerProc); // region-wise matrices with quasiRegion maps (= composite GIDs)
  std::vector< Epetra_CrsMatrix * > regionGrpMats(maxRegPerProc); // region-wise matrices in true region layout with unique GIDs for replicated interface DOFs

  std::vector< Epetra_Map* > coarseRowMapPerGrp(maxRegPerProc); // region-wise row map in true region layout with unique GIDs for replicated interface DOFs
  std::vector< Epetra_Map* > coarseColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector< Epetra_CrsMatrix * > regionGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs

  Epetra_Vector* compX = NULL; // initial guess for truly composite calculations
  Epetra_Vector* compY = NULL; // result vector for truly composite calculations
  Epetra_Vector* regYComp = NULL; // result vector in composite layout, but computed via regional operations

  std::vector<Epetra_Vector*> quasiRegX(maxRegPerProc); // initial guess associated with myRank's ith region in quasiRegional layout
  std::vector<Epetra_Vector*> quasiRegY(maxRegPerProc); // result vector associated with myRank's ith region in quasiRegional layout
  std::vector<Epetra_Vector*> regX(maxRegPerProc); // initial guess associated with myRank's ith region in regional layout
  std::vector<Epetra_Vector*> regY(maxRegPerProc); // result vector associated with myRank's ith region in regional layout

  std::vector<int> intIDs; // LIDs of interface DOFs

  // loop forever (or until the matlab program issues a terminate command
  //               which causes exit() to be invoked).

  while ( 1 == 1 ) {
    Comm.Barrier();

    // wait for file from Matlab program
    while ((fp = fopen(fileName, "r")) == NULL)
      sleep(1);

    fgets(command,40,fp); stripTrailingJunk(command);

    printf("%d: Matlab command is __%s__\n",myRank,command);
    if (strcmp(command,"LoadCompositeMap") == 0) {

      std::vector<int> fileData;        // composite GIDs
      int i;
      while ( fscanf(fp,"%d",&i) != EOF) fileData.push_back(i);

      mapComp= new Epetra_Map(-1,(int) fileData.size(),fileData.data(),0,Comm);
    }
    else if (strcmp(command,"PrintCompositeMap") == 0) {
       sleep(myRank);
       if (AComp == NULL) {
         printf("%d: printing owned part of composite map\n",myRank);
         const int *rowGIDs = mapComp->MyGlobalElements();
         for (int i = 0; i < mapComp->NumMyElements(); i++)
           printf("%d: compMap(%d) = %d\n",myRank,i,rowGIDs[i]);
       }
       else {
         printf("%d: printing extended composite map\n",myRank);
         const int *colGIDs= AComp->ColMap().MyGlobalElements();
         for (int i = 0; i < AComp->ColMap().NumMyElements(); i++)
           printf("%d: compMap(%d) = %d\n",myRank,i,colGIDs[i]);
       }
    }
    else if (strcmp(command,"LoadRegions") == 0) {
      while ( fgets(command,80,fp) != NULL ) {
        int i;
        sscanf(command,"%d",&i);
        myRegions.push_back(i);
      }
    }
    else if (strcmp(command,"ReadMatrix") == 0) {
       EpetraExt::MatrixMarketFileToCrsMatrix("Amat.mm", *mapComp, AComp);
    }
    else if (strcmp(command,"LoadAndCommRegAssignments") == 0) {

      Epetra_MultiVector *regionsPerGID= new Epetra_MultiVector(AComp->RowMap(),maxRegPerGID);

      int k;
      double *jthRegions; // pointer to jth column in regionsPerGID
      for (int i = 0; i < mapComp->NumMyElements(); i++) {
        for (int j = 0; j < maxRegPerGID; j++) {
          jthRegions = (*regionsPerGID)[j];

          if ( fscanf(fp,"%d",&k) == EOF) {
            fprintf(stderr,"not enough region assignments\n"); exit(1);
          }
          jthRegions[i] = (double) k;

          // identify interface DOFs. An interface DOF is assinged to at least two regions.
          if (j > 0 and k != -1)
            intIDs.push_back(i);
        }
      }

      // make extended Region Assignments
      Epetra_Import Importer(AComp->ColMap(), *mapComp);
      regionsPerGIDWithGhosts = new Epetra_MultiVector(AComp->ColMap(),maxRegPerGID);
      regionsPerGIDWithGhosts->Import(*regionsPerGID, Importer, Insert);
    }
    else if (strcmp(command,"PrintRegionAssignments") == 0) {
      double *jthRegions;
      sleep(myRank);
      printf("%d: printing out extended regionsPerGID\n",myRank);
      const int *colGIDsComp = AComp->ColMap().MyGlobalElements();
      for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {
        printf("compGID(%d,%d) = %d:",myRank,i,colGIDsComp[i]);
        for (int j = 0; j <  maxRegPerGID; j++) {
          jthRegions = (*regionsPerGIDWithGhosts)[j];
          printf("%d ",(int) jthRegions[i]);
        }
        printf("\n");
      }
    }
    else if (strcmp(command,"LoadAppDataForLIDregion()") == 0) {
      int minGID,maxGID;
      for (int i = 0; i < (int) myRegions.size(); i++) {
        fscanf(fp,"%d%d",&minGID,&maxGID);
        minGIDComp[i] = minGID;
        maxGIDComp[i] = maxGID;
      }
      appData.minGIDComp   = minGIDComp.data();
      appData.maxGIDComp   = maxGIDComp.data();
      appData.myRegions    = myRegions.data();
      appData.colMap       = (Epetra_Map *) &(AComp->ColMap());
      appData.maxRegPerGID = maxRegPerGID;
      appData.regionsPerGIDWithGhosts = regionsPerGIDWithGhosts;
    }
    else if (strcmp(command,"MakeGrpRegColMaps") == 0) {
      if (whichCase == MultipleRegionsPerProc) {
        // clone rowMap
        for (int j=0; j < maxRegPerProc; j++) {
          if (j < (int) myRegions.size()) {
            colMapPerGrp[j] = new Epetra_Map(-1,
                rowMapPerGrp[j]->NumMyElements(),
                rowMapPerGrp[j]->MyGlobalElements(), 0, Comm);
          }
          else colMapPerGrp[j] = new Epetra_Map(-1,0,NULL,0,Comm);
        }
      }
      else if (whichCase == RegionsSpanProcs) {//so maxRegPerProc = 1
        std::vector<int> colIDsReg;

        // copy the rowmap
        int *rowGIDsReg = rowMapPerGrp[0]->MyGlobalElements();
        for (int i = 0; i < rowMapPerGrp[0]->NumMyElements(); i++)
          colIDsReg.push_back(rowGIDsReg[i]);

        // append additional ghosts who are in my region and
        // for whom I have a LID
        int LID;
        double *jthRegions;
        int *colGIDsComp =  AComp->ColMap().MyGlobalElements();
        for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {
          LID = LIDregion(&appData, i, 0);
          if (LID == -1) {
            for (int j = 0; j < maxRegPerGID; j++) {
              jthRegions = (*regionsPerGIDWithGhosts)[j];
              if  ( ((int) jthRegions[i]) == myRegions[0]) {
                colIDsReg.push_back(colGIDsComp[i]);
                break;
              }
            }
          }
        }
        if ((int) myRegions.size() > 0) {
         colMapPerGrp[0] = new Epetra_Map(-1, (int) colIDsReg.size(),
             colIDsReg.data(), 0, Comm);
        }
        else colMapPerGrp[0] = new Epetra_Map(-1,0,NULL,0,Comm);
      }
      else { fprintf(stderr,"whichCase not set properly\n"); exit(1); }
    }
    else if (strcmp(command,"PrintGrpRegColMaps") == 0) {
      sleep(myRank);
      char str[80]; sprintf(str,"%d: colMap ",myRank);
      printGrpMaps(colMapPerGrp, maxRegPerProc, str);
    }
    else if (strcmp(command,"MakeExtendedGrpRegMaps") == 0) {
      int nLocal   = AComp->RowMap().NumMyElements();
      int nExtended= AComp->ColMap().NumMyElements();
      int nTotal = 0;
      Comm.SumAll(&nLocal,&nTotal,1);

      // first step of NewGID calculation just counts the number of NewGIDs
      // and sets firstNewGID[k] such that it is equal to the number of
      // NewGIDs that have already been counted in (0:k-1) that myRank
      // is responsible for.

      Epetra_Vector firstNewGID(*mapComp);
      firstNewGID[0] = 0.;
      double *jthRegions = NULL;
      for (int k = 0; k < nLocal-1; k++) {
        firstNewGID[k+1] = firstNewGID[k]-1.;
        for (int j = 0; j < maxRegPerGID; j++) {
          jthRegions = (*regionsPerGIDWithGhosts)[j];
          if (jthRegions[k] != -1) (firstNewGID[k+1]) += 1.;
        }
      }
      // So firstNewGID[nLocal-1] is number of NewGIDs up to nLocal-2
      // To account for newGIDs associated with nLocal-1, we'll just
      // use an upper bound (to avoid the little loop above).
      // By adding maxRegPerGID-1 we account for the maximum
      // number of possible newGIDs due to last composite id
      int upperBndNumNewGIDs = (int) firstNewGID[nLocal-1] + maxRegPerGID-1;
      int upperBndNumNewGIDsAllProcs;
      Comm.MaxAll(&upperBndNumNewGIDs,&upperBndNumNewGIDsAllProcs,1);

      // Now that we have an upper bound on the maximum number of
      // NewGIDs over all procs, we sweep through firstNewGID again
      // to assign ids to the first NewGID associated with each row of
      // regionsPerGIDWithGhosts (by adding an offset)
      for (int k = 0; k < nLocal; k++)
        firstNewGID[k] += upperBndNumNewGIDsAllProcs*myRank+nTotal;

      Epetra_Import Importer(AComp->ColMap(), *mapComp);
      Epetra_Vector firstNewGIDWithGhost(AComp->ColMap());
      firstNewGIDWithGhost.Import(firstNewGID, Importer, Insert);

      std::vector<int> revisedGIDs;

      for (int k = 0; k < (int) myRegions.size(); k++) {
        revisedGIDs.resize(0);
        int curRegion = myRegions[k];
        double *jthRegions;
        int *colGIDsComp =  AComp->ColMap().MyGlobalElements();
        std::vector<int> tempRegIDs(nExtended);

        // must put revisedGIDs in application-provided order by
        // invoking LIDregion() and sorting
        for (int i = 0; i < nExtended; i++)
          tempRegIDs[i] = LIDregion(&appData, i, k);
        std::vector<int> idx(tempRegIDs.size());
        std::iota(idx.begin(),idx.end(),0);
        sort(idx.begin(),idx.end(),[tempRegIDs](int i1,int i2){return tempRegIDs[i1] < tempRegIDs[i2];});

        // Now sweep through regionsPerGIDWithGhosts looking for those
        // associated with curRegion and record the revisedGID
        int j;
        for (int i = 0; i < nExtended; i++) {

          if (tempRegIDs[idx[i]] != -1) {// if a valid LID for this region
            for (j = 0; j < maxRegPerGID; j++) {
              jthRegions = (*regionsPerGIDWithGhosts)[j];
              if (((int) jthRegions[idx[i]]) == curRegion) break;
            }

            // (*regionsPerGIDWithGhosts)[0] entries keep original GID
            // while others use firstNewGID to determine NewGID

            if (j == 0) revisedGIDs.push_back(colGIDsComp[idx[i]]);
            else if (j < maxRegPerGID) {
               revisedGIDs.push_back((int) firstNewGIDWithGhost[idx[i]]+j-1);
            }
          }
        }

        revisedRowMapPerGrp[k] = new Epetra_Map(-1,(int) revisedGIDs.size(),
                                            revisedGIDs.data(),0,Comm);
        // now append more stuff to handle ghosts ... needed for
        // revised version of column map

        for (int i = 0; i < nExtended; i++) {
          if (tempRegIDs[i] == -1) {// only in revised col map
                     // note: since sorting not used when making the
                     // original regional colmap, we can't use
                     // it here either ... so no idx[]'s.
            for (j = 0; j < maxRegPerGID; j++) {
              jthRegions = (*regionsPerGIDWithGhosts)[j];
              if  ( ((int) jthRegions[i]) == curRegion) break;
            }
            // (*regionsPerGIDWithGhosts)[0] entries keep original GID
            // while others use firstNewGID to determine NewGID

            if (j == 0) revisedGIDs.push_back(colGIDsComp[i]);
            else if (j < maxRegPerGID) {
              revisedGIDs.push_back((int) firstNewGIDWithGhost[i]+j-1);
            }
          }
        }
        revisedColMapPerGrp[k] = new Epetra_Map(-1, (int) revisedGIDs.size(),
            revisedGIDs.data(), 0, Comm);
      }
      for (int k = (int) myRegions.size(); k < maxRegPerProc; k++) {
        revisedRowMapPerGrp[k] = new Epetra_Map(-1,0,NULL,0,Comm);
        revisedColMapPerGrp[k] = new Epetra_Map(-1,0,NULL,0,Comm);
      }

    }
    else if (strcmp(command,"MakeGrpRegRowMaps") == 0) {
      std::vector<int> rowGIDsReg;
      int *colGIDsComp = AComp->ColMap().MyGlobalElements();
      for (int k = 0; k < (int) myRegions.size(); k++) {
        rowGIDsReg.resize(0);
        std::vector<int> tempRegIDs(AComp->ColMap().NumMyElements());
        for (int i = 0; i < AComp->ColMap().NumMyElements(); i++)
          tempRegIDs[i] = LIDregion(&appData, i, k);

        std::vector<int> idx(tempRegIDs.size());
        std::iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(), [tempRegIDs](int i1,int i2) { return tempRegIDs[i1] < tempRegIDs[i2];});

        for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {
          if (tempRegIDs[idx[i]] != -1)
            rowGIDsReg.push_back(colGIDsComp[idx[i]]);
        }
        rowMapPerGrp[k] = new Epetra_Map(-1, (int) rowGIDsReg.size(),
            rowGIDsReg.data(), 0, Comm);
      }

      for (int k=(int) myRegions.size(); k < maxRegPerProc; k++) {
        rowMapPerGrp[k] = new Epetra_Map(-1,0,NULL,0,Comm);
      }
    }
    else if (strcmp(command,"PrintRevisedRowMaps") == 0) {
      sleep(myRank);
      char str[80]; sprintf(str,"%d: revisedRowMap ",myRank);
      printGrpMaps(revisedRowMapPerGrp, maxRegPerProc, str);
    }
    else if (strcmp(command,"PrintRevisedColMaps") == 0) {
      sleep(myRank);
      char str[80]; sprintf(str,"%d: revisedColMap ",myRank);
      printGrpMaps(revisedColMapPerGrp, maxRegPerProc, str);
    }
    else if (strcmp(command,"PrintGrpRegDomMaps") == 0) {
      sleep(myRank);
      char str[80]; sprintf(str,"%d: domMap ",myRank);
      printGrpMaps(revisedRowMapPerGrp, maxRegPerProc, str);
    }
    else if (strcmp(command,"PrintGrpRegRowMaps") == 0) {
      sleep(myRank);
      char str[80]; sprintf(str,"%d: rowMap ",myRank);
      printGrpMaps(rowMapPerGrp, maxRegPerProc, str);
    }
    else if (strcmp(command,"MakeQuasiRegionMatrices") == 0) {
      // copy and modify the composite matrix. Extract diagonal and divide interface entries by 2
      ACompSplit = new Epetra_CrsMatrix(*AComp);
      Epetra_Vector* diagAComp = new Epetra_Vector(*mapComp, true);
      ACompSplit->ExtractDiagonalCopy(*diagAComp);

      for (int i = 0; i < intIDs.size(); ++i)
        (*diagAComp)[intIDs[i]] *= 0.5;

      ACompSplit->ReplaceDiagonalValues(*diagAComp);

      // Import data from ACompSplit into the quasiRegion matrices
      for (int j = 0; j < maxRegPerProc; j++) {
        rowImportPerGrp[j] = new Epetra_Import(*(rowMapPerGrp[j]),*mapComp);
        colImportPerGrp[j] = new Epetra_Import(*(colMapPerGrp[j]),*mapComp);
        quasiRegionGrpMats[j]     = new Epetra_CrsMatrix(Copy,*(rowMapPerGrp[j]),
                                                     *(colMapPerGrp[j]),-1);
        quasiRegionGrpMats[j]->Import(*ACompSplit, *(rowImportPerGrp[j]), Insert);
        quasiRegionGrpMats[j]->FillComplete();
      }
    }
    else if (strcmp(command,"PrintQuasiRegionMatrices") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: Matrix in Grp %d\n",myRank,j);
        std::cout << *(quasiRegionGrpMats[j]) << std::endl;
      }
    }
    else if (strcmp(command,"MakeRegionMatrices") == 0) {
      // We work on a copy. Just for safety.
      for (int j = 0; j < maxRegPerProc; j++) {
        // create empty matrix with correct row and column map
        regionGrpMats[j] = new Epetra_CrsMatrix(Copy, *(revisedRowMapPerGrp[j]), *(revisedColMapPerGrp[j]), 3);

        // extract pointers to crs arrays from quasiRegion matrix
        Epetra_IntSerialDenseVector & qRowPtr = quasiRegionGrpMats[j]->ExpertExtractIndexOffset();
        Epetra_IntSerialDenseVector & qColInd = quasiRegionGrpMats[j]->ExpertExtractIndices();
        double *& qVals = quasiRegionGrpMats[j]->ExpertExtractValues();

        // extract pointers to crs arrays from region matrix
        Epetra_IntSerialDenseVector & rowPtr = regionGrpMats[j]->ExpertExtractIndexOffset();
        Epetra_IntSerialDenseVector & colInd = regionGrpMats[j]->ExpertExtractIndices();
        double *& vals = regionGrpMats[j]->ExpertExtractValues();

        // assign array values from quasiRegional to regional matrices
        rowPtr.Resize(qRowPtr.Length());
        colInd.Resize(qColInd.Length());
        delete [] vals;
        vals = qVals;
        for (int i = 0; i < rowPtr.Length(); ++i) rowPtr[i] = qRowPtr[i];
        for (int i = 0; i < colInd.Length(); ++i) colInd[i] = qColInd[i];
      }

      /* add domain and range map to region matrices (pass in revisedRowMap, since
       * we assume that row map = range map = domain map)
       */
      for (int j = 0; j < maxRegPerProc; j++) {
        regionGrpMats[j]->ExpertStaticFillComplete(*(revisedRowMapPerGrp[j]),*(revisedRowMapPerGrp[j]));
      }
    }
    else if (strcmp(command,"PrintRegionMatrices") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: Matrix in Grp %d\n",myRank,j);
        std::cout << myRank << *(regionGrpMats[j]) << std::endl;
      }
    }
    else if (strcmp(command,"PrintRegionMatrixRowMap") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: Row Map Grp %d\n", myRank, j);
        regionGrpMats[j]->RowMap().Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegionMatrixColMap") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: Col Map Grp %d\n", myRank, j);
        regionGrpMats[j]->ColMap().Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegionMatrixRangeMap") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: Range Map Grp %d\n", myRank, j);
        regionGrpMats[j]->OperatorRangeMap().Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegionMatrixDomainMap") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: Domain Map Grp %d\n", myRank, j);
        regionGrpMats[j]->OperatorDomainMap().Print(std::cout);
      }
    }
    else if (strcmp(command,"ComputeMatVecs") == 0) {
      // create input and result vector in composite layout
      compX = new Epetra_Vector(*mapComp, true);
      compY = new Epetra_Vector(*mapComp, true);
      compX->Random();
//      compX->PutScalar(myRank);

      // perform matvec in composite layout
      AComp->Apply(*compX, *compY);

      // transform composite vector to regional layout
      compositeToRegional(compX, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);
      createRegionalVector(regY, maxRegPerProc, revisedRowMapPerGrp);

      // perform matvec
      for (int j = 0; j < maxRegPerProc; j++) {
        int err = regionGrpMats[j]->Apply(*(regX[j]), *(regY[j]));
        TEUCHOS_ASSERT(err == 0);
      }

      // transform regY back to composite layout
      regYComp = new Epetra_Vector(*mapComp, true);
      for (int j = 0; j < maxRegPerProc; j++) {
        // copy vector and replace map
        quasiRegY[j] = new Epetra_Vector(*(regY[j]));
        quasiRegY[j]->ReplaceMap(*(rowMapPerGrp[j]));

        int err = regYComp->Export(*quasiRegY[j], *(rowImportPerGrp[j]), Epetra_AddLocalAlso);
        TEUCHOS_ASSERT(err == 0);
      }

      Epetra_Vector* diffY = new Epetra_Vector(*mapComp, true);
      diffY->Update(1.0, *compY, -1.0, *regYComp, 0.0);
      double diffNormTwo = 0.0;
      double diffNormInf = 0.0;
      diffY->Norm2(&diffNormTwo);
      diffY->NormInf(&diffNormInf);
      std::cout << myRank << ": diffNormTwo: " << diffNormTwo << "\tdiffNormInf: " << diffNormInf << std::endl;

    }
    else if (strcmp(command,"MakeRegionTransferOperators") == 0) {
      /* Populate a fine grid vector with 1 at coarse nodes and 0 at fine nodes.
       * Then, transform to regional layout and find GIDs with entry 1
       */
      int numNodes = mapComp->NumGlobalPoints() / 3 + 1;
      double vals[numNodes];
      int ind[numNodes];
      for (int i = 0; i < numNodes; ++i)
      {
        vals[i] = 1.0;
        ind[i] = 3*i;
      }
      Epetra_Vector* coarseGridToggle = new Epetra_Vector(*mapComp, true);
      coarseGridToggle->ReplaceGlobalValues(numNodes, vals, ind);

      // transform composite vector to regional layout
      std::vector<Epetra_Vector*> quasiRegCoarseGridToggle(maxRegPerProc);
      std::vector<Epetra_Vector*> regCoarseGridToggle(maxRegPerProc);
      compositeToRegional(coarseGridToggle, quasiRegCoarseGridToggle,
          regCoarseGridToggle, maxRegPerProc, rowMapPerGrp, revisedRowMapPerGrp,
          rowImportPerGrp);

      // Print regCoarseGridToggle
//      sleep(myRank);
//      for (int j = 0; j < maxRegPerProc; j++) {
//        printf("%d: regCoarseGridToggle %d\n", myRank, j);
//        regCoarseGridToggle[j]->Print(std::cout);
//      }

      // create coarse grid row maps in region layout
      for (int j = 0; j < maxRegPerProc; j++) {
        std::vector<int> coarseRowGIDsReg;
        for (int i = 0; i < regCoarseGridToggle[j]->Map().NumMyElements(); ++i) {
          if ((*regCoarseGridToggle[j])[i] == 1.0)
            coarseRowGIDsReg.push_back(regCoarseGridToggle[j]->Map().GID(i));
        }

        coarseRowMapPerGrp[j] = new Epetra_Map(-1, coarseRowGIDsReg.size(), coarseRowGIDsReg.data(), 0, Comm);
//        printf("%d: coarseRowMapPerGrp %d\n", myRank, j);
//        coarseRowMapPerGrp[j]->Print(std::cout);
      }

      // create coarse grid column map
      Epetra_Vector* compGIDVec = new Epetra_Vector(*mapComp, true);
      for (int i = 0; i < mapComp->NumMyElements(); ++i)
      {
        int myGID = mapComp->GID(i);
        if (myGID % 3 == 0) // that is the coarse point
          (*compGIDVec)[i] = myGID;
        else if (myGID % 3 == 1) // this node is right of a coarse point
          (*compGIDVec)[i] = myGID - 1;
        else if (myGID % 3 == 2) // this node is left of a coarse point
          (*compGIDVec)[i] = myGID + 1;
        else
          TEUCHOS_ASSERT(false)
      }
//      compGIDVec->Print(std::cout);

      // transform composite vector to regional layout
      std::vector<Epetra_Vector*> quasiRegGIDVec(maxRegPerProc);
      std::vector<Epetra_Vector*> regGIDVec(maxRegPerProc);
      compositeToRegional(compGIDVec, quasiRegGIDVec, regGIDVec, maxRegPerProc,
          rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

      // deal with duplicated nodes and replace with new GID when necessary
      for (int j = 0; j < maxRegPerProc; j++) {
        for (int i = 0; i < regGIDVec[j]->MyLength(); ++i) {
          if (regGIDVec[j]->Map().GID(i) >= mapComp->NumGlobalElements()) {
            // replace 'aggregate ID' with duplicated 'aggregate ID'
            (*regGIDVec[j])[i] = (revisedRowMapPerGrp[j])->GID(i);

            // look for other nodes with the same fine grid aggregate ID that needs to be replaced with the duplicated one
            for (int k = 0; k < regGIDVec[j]->MyLength(); ++k) {
              if ((*regGIDVec[j])[k] == (rowMapPerGrp[j])->GID(i))
                (*regGIDVec[j])[k] = (revisedRowMapPerGrp[j])->GID(i);
            }
          }
        }
      }

//      // Print regGIDVec
//      sleep(myRank);
//      for (int j = 0; j < maxRegPerProc; j++) {
//        printf("%d: regGIDVec %d\n", myRank, j);
//        regGIDVec[j]->Print(std::cout);
//      }

      for (int j = 0; j < maxRegPerProc; j++) {
        std::vector<int> regColGIDs;
        for (int i = 0; i < revisedRowMapPerGrp[j]->NumMyElements(); ++i) {
          const int currGID = revisedRowMapPerGrp[j]->GID(i);

          // loop over existing regColGIDs and see if current GID is already in that list
          bool found = false;
          for (int k = 0; k < regColGIDs.size(); ++k) {
            if ((*regGIDVec[j])[regGIDVec[j]->Map().LID(currGID)] == regColGIDs[k]) {
              found = true;
              break;
            }
          }

          if (not found)
            regColGIDs.push_back((*regGIDVec[j])[regGIDVec[j]->Map().LID(currGID)]);
        }

        coarseColMapPerGrp[j] = new Epetra_Map(-1, regColGIDs.size(), regColGIDs.data(), 0, Comm);
      }

//      sleep(myRank);
//      for (int j = 0; j < maxRegPerProc; j++) {
//        printf("%d: coarseColMapPerGrp %d\n", myRank, j);
//        coarseColMapPerGrp[j]->Print(std::cout);
//      }

      // Build the actual prolongator
      for (int j = 0; j < maxRegPerProc; j++) {
        regionGrpProlong[j] = new Epetra_CrsMatrix(Copy, *revisedRowMapPerGrp[j], *coarseColMapPerGrp[j], 1, false);
        for (int c = 0; c < regionGrpProlong[j]->NumMyCols(); ++c) {
          for (int r = 0; r < regionGrpProlong[j]->NumMyRows(); ++r) {
            if (regionGrpProlong[j]->ColMap().GID(c) == (*regGIDVec[j])[r]) {
              double vals[1];
              int inds[1];
              vals[0] = 1.0/3.0;
              inds[0] = c;
              regionGrpProlong[j]->InsertMyValues(r, 1, &*vals, &*inds);
            }
          }
        }
        regionGrpProlong[j]->FillComplete();
//        regionGrpProlong[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"RunTwoLevelMethod") == 0) {
      // intial guess for solution
      compX = new Epetra_Vector(*mapComp, true);

      // forcing vector
      Epetra_Vector* compB = new Epetra_Vector(*mapComp, true);
      compB->ReplaceGlobalValue(compB->GlobalLength() - 1, 0, 1.0);

      // residual vector
      Epetra_Vector* compRes = new Epetra_Vector(*mapComp, true);
      AComp->Multiply(false, *compX, *compRes);
      compRes->Update(1.0, *compB, -1.0);

      // -----------------------------------------------------------------------
      // Compute coarse grid operators (RAP) in region layout
      // -----------------------------------------------------------------------
      // Compute P'*A*P
      std::vector<Epetra_CrsMatrix*> regCoarseMatPerGrp(maxRegPerProc); // store coarse RAP
      std::vector<Epetra_CrsMatrix*> regCoarseTmpMatPerGrp(maxRegPerProc); // store intermediate result P'*A
      for (int j = 0; j < maxRegPerProc; j++) {
        regCoarseTmpMatPerGrp[j] = new Epetra_CrsMatrix(Copy, regionGrpProlong[j]->RowMap(), regionGrpMats[j]->ColMap(), 3, false);

        int err = EpetraExt::MatrixMatrix::Multiply(*regionGrpProlong[j], true, *regionGrpMats[j], false, *regCoarseTmpMatPerGrp[j], false);
        TEUCHOS_ASSERT(err == 0);
        regCoarseTmpMatPerGrp[j]->FillComplete();

        regCoarseMatPerGrp[j] = new Epetra_CrsMatrix(Copy, regCoarseTmpMatPerGrp[j]->RowMap(), regionGrpProlong[j]->ColMap(), 3, false);

        err = EpetraExt::MatrixMatrix::Multiply(*regCoarseTmpMatPerGrp[j], false, *regionGrpProlong[j], false, *regCoarseMatPerGrp[j], false);
        TEUCHOS_ASSERT(err == 0);
        regCoarseMatPerGrp[j]->FillComplete();
      }

      // transform composite vectors to regional layout
      compositeToRegional(compX, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);
      std::vector<Epetra_Vector*> quasiRegRes(maxRegPerProc);
      std::vector<Epetra_Vector*> regRes(maxRegPerProc);
      compositeToRegional(compRes, quasiRegRes, regRes, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);

//      // pre-smoothing on fine level
//      const int maxIter = 3;
//      for (int iter = 0; iter < maxIter; ++iter) {
//        for (int j = 0; j < maxRegPerProc; j++) {
//          Epetra_Vector* diag = new Epetra_Vector(regionGrpMats[j]->RowMap(), true);
//          Epetra_Vector* invDiag = new Epetra_Vector(regionGrpMats[j]->RowMap(), true);
//          regionGrpMats[j]->ExtractDiagonalCopy(diag);
//          diag->Reciprocal(invDiag);
//          for (int i = 0; i < regX[j]->MyLength(); ++i) {
//
//          }
//        }
//      }


    }
    else if (strcmp(command,"PrintCompositeVectorX") == 0) {
      sleep(myRank);
      compX->Print(std::cout);
    }
    else if (strcmp(command,"PrintCompositeVectorY") == 0) {
      sleep(myRank);
      compY->Print(std::cout);
    }
    else if (strcmp(command,"PrintQuasiRegVectorX") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: quasiRegion vector X Grp %d\n", myRank, j);
        quasiRegX[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegVectorX") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: region vector X Grp %d\n", myRank, j);
        regX[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegVectorY") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: region vector Y Grp %d\n", myRank, j);
        regY[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegVectorYComp") == 0) {
      sleep(myRank);
      regYComp->Print(std::cout);
    }
    else if (strcmp(command,"Terminate") == 0) {
      sprintf(command,"/bin/rm -f %s",fileName); system(command); exit(1);
    }
    else { printf("%d: Command (%s) not recognized !\n",myRank,command);  exit(1); }
    sprintf(command,"/bin/rm -f myData_%d",myRank); system(command); sleep(1);
  }
}

// returns local ID (within region curRegion) for the LIDcomp^th composite grid
// point that myRank owns. If this grid point is not part of curRegion, then
// -1 is returned.
//
// Code currently hardwired for 1D only.

int LIDregion(void *ptr, int LIDcomp, int whichGrp)
{
   struct widget * myWidget = (struct widget *) ptr;

   int        *minGIDComp  = myWidget->minGIDComp;
   int        *maxGIDComp  = myWidget->maxGIDComp;
   int        *myRegions   = myWidget->myRegions;
   Epetra_Map *colMap      = myWidget->colMap;
   int        maxRegPerGID = myWidget->maxRegPerGID;
   Epetra_MultiVector *regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

   int curRegion = myRegions[whichGrp];
   if (LIDcomp >=  colMap->NumMyElements()) return(-1);

   double *jthRegions;
   const int *colGIDsComp = colMap->MyGlobalElements();

   if (colGIDsComp[LIDcomp] < minGIDComp[whichGrp]) return(-1);
   if (colGIDsComp[LIDcomp] > maxGIDComp[whichGrp]) return(-1);

   bool found = false;
   for (int j = 0; j <  maxRegPerGID; j++) {
      jthRegions = (*regionsPerGIDWithGhosts)[j];
      if  ( ((int) jthRegions[LIDcomp]) == curRegion ) {
         found = true;
         break;
      }
   }
   if (found == false) return(-1);

   return( colGIDsComp[LIDcomp] - minGIDComp[whichGrp] );
}



// just some junk code to remove things like space and newlines when we
// read strings from files ... so that strcmp() works properly
void stripTrailingJunk(char *command)
{
   int i  = strlen(command)-1;
   while ( (command[i] == ' ') || (command[i] == '\n')) {
      command[i] = '\0';
      i--;
   }
}
// prints Grp-style maps
void printGrpMaps(std::vector < Epetra_Map* > &mapPerGrp, int maxRegPerProc, char *str)
{
   for (int j = 0; j < maxRegPerProc; j++) {
      printf("%s in grp(%d):",str,j);
      Epetra_Map *jthMap = mapPerGrp[j];
      const int *GIDsReg = jthMap->MyGlobalElements();
      for (int i = 0; i < mapPerGrp[j]->NumMyElements(); i++) {
         printf("%d ",GIDsReg[i]);
      }
      printf("\n");
   }
}
