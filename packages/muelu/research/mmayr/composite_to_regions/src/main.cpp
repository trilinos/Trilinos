#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <numeric>

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
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_Transpose_RowMatrix.h"

#include <MueLu.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_NullspaceFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>

#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>

#include "Teuchos_Assert.hpp"

int LIDregion(void *ptr, int LIDcomp, int whichGrp);
int LID2Dregion(void *ptr, int LIDcomp, int whichGrp);

// this little widget handles application specific data
// used to implement LIDregion()
struct widget {
   int *minGIDComp;
   int *maxGIDComp;
   int *myRegions;
   Epetra_Map *colMap;
   int maxRegPerGID;
   Teuchos::RCP<Epetra_MultiVector> regionsPerGIDWithGhosts;
   int *gDim, *lDim, *lowInd;
   int       *trueCornerx; // global coords of region
   int       *trueCornery; // corner within entire 2D mesh
   int       *relcornerx;  // coords of corner relative
   int       *relcornery;  // to region corner
   int       *lDimx;
   int       *lDimy;
   int        nx;
int myRank;
};

void stripTrailingJunk(char *command);
void printGrpMaps(std::vector<Teuchos::RCP<Epetra_Map> > &mapPerGrp, int maxRegPerProc, char *str);

/*! \brief Transform composite vector to regional layout
 */
void compositeToRegional(Teuchos::RCP<const Epetra_Vector> compVec, ///< Vector in composite layout [in]
    std::vector<Teuchos::RCP<Epetra_Vector> >& quasiRegVecs, ///< Vector in quasiRegional layout [in/out]
    std::vector<Teuchos::RCP<Epetra_Vector> >& regVecs, ///< Vector in regional layout [in/out]
    const int maxRegPerProc, ///< max number of regions per proc [in]
    const std::vector<Teuchos::RCP<Epetra_Map> > rowMapPerGrp, ///< row maps in region layout [in]
    const std::vector<Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp, ///< revised row maps in region layout [in]
    const std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp ///< row importer in region layout [in]
    )
{
  // quasiRegional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create empty vectors and fill it by extracting data from composite vector
    quasiRegVecs[j] = Teuchos::rcp(new Epetra_Vector(*(rowMapPerGrp[j]), true));
    TEUCHOS_ASSERT(!quasiRegVecs[j].is_null());
    int err = quasiRegVecs[j]->Import(*compVec, *(rowImportPerGrp[j]), Insert);
    TEUCHOS_ASSERT(err == 0);
  }

  // regional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create regVecs vector (copy from quasiRegVecs and swap the map)
    regVecs[j] = Teuchos::rcp(new Epetra_Vector(*(quasiRegVecs[j])));
    TEUCHOS_ASSERT(!regVecs[j].is_null());
    int err = regVecs[j]->ReplaceMap(*(revisedRowMapPerGrp[j]));
    TEUCHOS_ASSERT(err == 0);
  }

  return;
}

/*! \brief Transform regional vector to composite layout
 */
void regionalToComposite(const std::vector<Teuchos::RCP<Epetra_Vector> >& regVec, ///< Vector in region layout [in]
    Teuchos::RCP<Epetra_Vector> compVec, ///< Vector in composite layout [in/out]
    const int maxRegPerProc, ///< max number of regions per proc
    const std::vector<Teuchos::RCP<Epetra_Map> > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
    const std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp, ///< row importer in region layout [in]
    const Epetra_CombineMode combineMode, ///< Combine mode for import/export [in]
    const bool addLocalManually ///< perform ADD of local values manually (not via Epetra CombineMode)
    )
{
  if (not addLocalManually) {
    /* Use the Eptra_AddLocalAlso combine mode to add processor-local values.
     * Note that such a combine mode is not available in Tpetra.
     */

    std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegVec(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      // copy vector and replace map
      quasiRegVec[j] = Teuchos::rcp(new Epetra_Vector(*(regVec[j])));
      int err = quasiRegVec[j]->ReplaceMap(*(rowMapPerGrp[j]));
      TEUCHOS_ASSERT(err == 0);

      err = compVec->Export(*quasiRegVec[j], *(rowImportPerGrp[j]), combineMode);
      TEUCHOS_ASSERT(err == 0);
    }
  }
  else {
    /* Let's fake an ADD combine mode that also adds local values by
     * 1. exporting quasiRegional vectors to auxiliary composite vectors (1 per group)
     * 2. add all auxiliary vectors together
     */

    Teuchos::RCP<Epetra_MultiVector> partialCompVec = Teuchos::rcp(new Epetra_MultiVector(compVec->Map(), maxRegPerProc, true));
    TEUCHOS_ASSERT(!partialCompVec.is_null());

    std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegVec(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      // copy vector and replace map
      quasiRegVec[j] = Teuchos::rcp(new Epetra_Vector(*(regVec[j])));
      TEUCHOS_ASSERT(!quasiRegVec[j].is_null());

      int err = quasiRegVec[j]->ReplaceMap(*(rowMapPerGrp[j]));
      TEUCHOS_ASSERT(err == 0);

      err = (*partialCompVec)(j)->Export(*quasiRegVec[j], *(rowImportPerGrp[j]), Add);
      TEUCHOS_ASSERT(err == 0);
    }

    compVec->PutScalar(0.0);
    for (int j = 0; j < maxRegPerProc; j++) {
      int err = compVec->Update(1.0, *(*partialCompVec)(j), 1.0);
      TEUCHOS_ASSERT(err == 0);
    }
  }

  return;
}

void printRegionalVector(const std::string vectorName, ///< string to be used for screen output
    const std::vector<Teuchos::RCP<Epetra_Vector> > regVecs, ///< regional vector to be printed to screen
    const int myRank ///< rank of calling proc
    )
{
//  sleep(myRank);
  for (int j = 0; j < (int) regVecs.size(); j++) {
    printf("%d: %s %d\n", myRank, vectorName.c_str(), j);
    regVecs[j]->Print(std::cout);
  }
}

void printRegionalMap(const std::string mapName, ///< string to be used for screen output
    const std::vector<Teuchos::RCP<Epetra_Map> > regMaps, ///< regional map to be printed to screen
    const int myRank ///< rank of calling proc
    )
{
  sleep(myRank);
  for (int j = 0; j < (int) regMaps.size(); j++) {
    printf("%d: %s %d\n", myRank, mapName.c_str(), j);
    regMaps[j]->Print(std::cout);
  }
}

/*! \brief Sum region interface values
 */
void sumInterfaceValues(std::vector<Teuchos::RCP<Epetra_Vector> >& regVec,
    Teuchos::RCP<const Epetra_Map> compMap,
    const int maxRegPerProc, ///< max number of regions per proc [in]
    const std::vector<Teuchos::RCP<Epetra_Map> > rowMapPerGrp,///< row maps in region layout [in]
    const std::vector<Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp,///< revised row maps in region layout [in]
    const std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp ///< row importer in region layout [in])
    )
{
  Teuchos::RCP<Epetra_Vector> compVec = Teuchos::rcp(new Epetra_Vector(*compMap, true));
  TEUCHOS_ASSERT(!compVec.is_null());

  std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegVec(maxRegPerProc);
  regionalToComposite(regVec, compVec, maxRegPerProc, rowMapPerGrp,
      rowImportPerGrp, Add, true);

  compositeToRegional(compVec, quasiRegVec, regVec, maxRegPerProc,
      rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

  return;
}

void createRegionalVector(std::vector<Teuchos::RCP<Epetra_Vector> >& regVecs,
    const int maxRegPerProc, const std::vector<Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp)
{
  regVecs.resize(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++)
    regVecs[j] = Teuchos::rcp(new Epetra_Vector(*(revisedRowMapPerGrp[j]), true));
  return;
}

std::vector<Teuchos::RCP<Epetra_Vector> > computeResidual(
    std::vector<Teuchos::RCP<Epetra_Vector> >& regRes, ///< residual (to be evaluated)
    const std::vector<Teuchos::RCP<Epetra_Vector> > regX, ///< left-hand side (solution)
    const std::vector<Teuchos::RCP<Epetra_Vector> > regB, ///< right-hand side (forcing term)
    const std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regionGrpMats,
    Teuchos::RCP<const Epetra_Map> mapComp, ///< composite map, computed by removing GIDs > numDofs in revisedRowMapPerGrp
    const std::vector<Teuchos::RCP<Epetra_Map> > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
    const std::vector<Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
    const std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp ///< row importer in region layout [in]
    )
{
  const int maxRegPerProc = regX.size();

  /* Update the residual vector
   * 1. Compute tmp = A * regX in each region
   * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
   * 3. Compute r = B - tmp
   */
  for (int j = 0; j < maxRegPerProc; j++) { // step 1
    int err = regionGrpMats[j]->Multiply(false, *regX[j], *regRes[j]);
    TEUCHOS_ASSERT(err == 0);
    TEUCHOS_ASSERT(regionGrpMats[j]->DomainMap().PointSameAs(regX[j]->Map()));
    TEUCHOS_ASSERT(regionGrpMats[j]->RangeMap().PointSameAs(regRes[j]->Map()));
  }

  sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
      revisedRowMapPerGrp, rowImportPerGrp);

  for (int j = 0; j < maxRegPerProc; j++) { // step 3
    int err = regRes[j]->Update(1.0, *regB[j], -1.0);
    TEUCHOS_ASSERT(err == 0);
    TEUCHOS_ASSERT(regRes[j]->Map().PointSameAs(regB[j]->Map()));
  }

  return regRes;
}

void jacobiIterate(const int maxIter,
    const double omega,
    std::vector<Teuchos::RCP<Epetra_Vector> >& regX, // left-hand side (or solution)
    const std::vector<Teuchos::RCP<Epetra_Vector> > regB, // right-hand side (or residual)
    const std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regionGrpMats, // matrices in true region layout
    const std::vector<Teuchos::RCP<Epetra_Vector> > regionInterfaceScaling, // recreate on coarse grid by import Add on region vector of ones
    const int maxRegPerProc, ///< max number of regions per proc [in]
    Teuchos::RCP<const Epetra_Map> mapComp, ///< composite map
    const std::vector<Teuchos::RCP<Epetra_Map> > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
    const std::vector<Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
    const std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp ///< row importer in region layout [in]
    )
{
  std::vector<Teuchos::RCP<Epetra_Vector> > regRes(maxRegPerProc);
  createRegionalVector(regRes, maxRegPerProc, revisedRowMapPerGrp);

  // extract diagonal from region matrices and recover true diagonal values
  std::vector<Teuchos::RCP<Epetra_Vector> > diag(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    // extract inverse of diagonal from matrix
    diag[j] = Teuchos::rcp(new Epetra_Vector(regionGrpMats[j]->RowMap(), true));
    TEUCHOS_ASSERT(!diag[j].is_null());
    int err = regionGrpMats[j]->ExtractDiagonalCopy(*diag[j]);
    TEUCHOS_ASSERT(err == 0);
    for (int i = 0; i < diag[j]->MyLength(); ++i) // ToDo: replace this by an Epetra_Vector routine
      (*diag[j])[i] *= (*regionInterfaceScaling[j])[i]; // Scale to obtain the true diagonal
  }

  for (int iter = 0; iter < maxIter; ++iter) {

    /* Update the residual vector
     * 1. Compute tmp = A * regX in each region
     * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
     * 3. Compute r = B - tmp
     */
    for (int j = 0; j < maxRegPerProc; j++) { // step 1
      TEUCHOS_ASSERT(regionGrpMats[j]->OperatorDomainMap().PointSameAs(regX[j]->Map()));
      TEUCHOS_ASSERT(regionGrpMats[j]->OperatorRangeMap().PointSameAs(regRes[j]->Map()));

      int err = regionGrpMats[j]->Multiply(false, *regX[j], *regRes[j]);
      TEUCHOS_ASSERT(err == 0);
    }

    sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
        revisedRowMapPerGrp, rowImportPerGrp);

    for (int j = 0; j < maxRegPerProc; j++) { // step 3
      int err = regRes[j]->Update(1.0, *regB[j], -1.0);
      TEUCHOS_ASSERT(err == 0);
    }

    // check for convergence
    {
      Teuchos::RCP<Epetra_Vector> compRes = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
          rowImportPerGrp, Add, true);
      TEUCHOS_ASSERT(compRes->Map().UniqueGIDs());
      double normRes = 0.0;
      int err = compRes->Norm2(&normRes);
      TEUCHOS_ASSERT(err == 0);

      if (normRes < 1.0e-12)
        return;
    }

    for (int j = 0; j < maxRegPerProc; j++) {
      // update solution according to Jacobi's method
      for (int i = 0; i < regX[j]->MyLength(); ++i) {
        (*regX[j])[i] += omega / (*diag[j])[i] * (*regRes[j])[i];
      }
    }
  }

  return;
}

/*! \brief Find common regions of two nodes
 *
 */
std::vector<int> findCommonRegions(const int nodeA, ///< GID of first node
    const int nodeB, ///< GID of second node
    const Epetra_MultiVector& nodesToRegions ///< mapping of nodes to regions
    )
{
  // extract node-to-regions mapping for both nodes A and B
  std::vector<int> regionsA;
  std::vector<int> regionsB;
  {
    const Epetra_BlockMap& map = nodesToRegions.Map();
    for (int i = 0; i < nodesToRegions.NumVectors(); ++i) {
      regionsA.push_back((nodesToRegions[i])[map.LID(nodeA)]);
      regionsB.push_back((nodesToRegions[i])[map.LID(nodeB)]);
    }
  }

//  // Print list of regions for both nodes
//  {
//    int myRank = nodesToRegions.Comm().MyPID();
//    const Epetra_BlockMap& map = nodesToRegions.Map();
//    std::cout << myRank << ": nodeA = " << map.GID(nodeA) << ": ";
//    for (int i = 0; i < nodesToRegions.NumVectors(); ++i)
//      std::cout << ", " << regionsA[i];
//    std::cout << std::endl;
//    std::cout << myRank << ": nodeB = " << map.GID(nodeB) << ": ";
//      for (int i = 0; i < nodesToRegions.NumVectors(); ++i)
//        std::cout << ", " << regionsB[i];
//      std::cout << std::endl;
//  }

  // identify common regions
  std::vector<int> commonRegions(nodesToRegions.NumVectors());
  std::sort(regionsA.begin(), regionsA.end());
  std::sort(regionsB.begin(), regionsB.end());

  std::vector<int>::iterator it = std::set_intersection(regionsA.begin(),
      regionsA.end(), regionsB.begin(), regionsB.end(), commonRegions.begin());
  commonRegions.resize(it-commonRegions.begin());

  // remove '-1' entries
  std::vector<int> finalCommonRegions;
  for (std::size_t i = 0; i < commonRegions.size(); ++i) {
    if (commonRegions[i] != -1)
      finalCommonRegions.push_back(commonRegions[i]);
  }

//  // Print result for debugging purposes
//  {
//    int myRank = nodesToRegions.Comm().MyPID();
//    std::cout << myRank << ": " << nodeA << "/" << nodeB << ": ";
//    for (int i = 0; i < finalCommonRegions.size(); ++i)
//      std::cout << ", " << finalCommonRegions[i];
//    std::cout << std::endl;
//  }

  return finalCommonRegions;
}

void vCycle(const int l, ///< ID of current level
    const int numLevels, ///< Total number of levels
    const int maxFineIter, ///< max. sweeps on fine and intermediate levels
    const int maxCoarseIter, ///< max. sweeps on coarse level
    const double omega, ///< damping parameter for Jacobi smoother
    const int maxRegPerProc, ///< Max number of regions per process
    std::vector<Teuchos::RCP<Epetra_Vector> >& fineRegX, ///< solution
    std::vector<Teuchos::RCP<Epetra_Vector> > fineRegB, ///< right hand side
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_CrsMatrix> > > regMatrices, ///< Matrices in region layout
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_CrsMatrix> > > regProlong, ///< Prolongators in region layout
    Teuchos::Array<Teuchos::RCP<Epetra_Map> > compMaps, ///< composite maps
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Map> > > quasiRegRowMaps, ///< quasiRegional row maps
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Map> > > regRowMaps, ///< regional row maps
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Import> > > regRowImporters, ///< regional row importers
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Vector> > > regInterfaceScalings ///< regional interface scaling factors
    )
{
  using Teuchos::RCP;
  using Teuchos::rcp;

//  const int myRank = compMaps[l]->Comm().MyPID();
//  std::cout << myRank << ": Entering V-cycle on level " << l << std::endl;

  if (l < numLevels - 1) { // fine or intermediate levels

    // pre-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
        regInterfaceScalings[l], maxRegPerProc, compMaps[l],
        quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    std::vector<Teuchos::RCP<Epetra_Vector> > regRes(maxRegPerProc);
    createRegionalVector(regRes, maxRegPerProc, regRowMaps[l]);
    computeResidual(regRes, fineRegX, fineRegB, regMatrices[l], compMaps[l],
        quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    // Transfer to coarse level
    std::vector<Teuchos::RCP<Epetra_Vector> > coarseRegX(maxRegPerProc);
    std::vector<Teuchos::RCP<Epetra_Vector> > coarseRegB(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      coarseRegX[j] = Teuchos::rcp(new Epetra_Vector(*regRowMaps[l+1][j], true));
      coarseRegB[j] = Teuchos::rcp(new Epetra_Vector(*regRowMaps[l+1][j], true));

      for (int i = 0; i < regRes[j]->MyLength(); ++i)
        (*regRes[j])[i] /= (*((regInterfaceScalings[l])[j]))[i];

      int err = regProlong[l+1][j]->Multiply(true, *regRes[j], *coarseRegB[j]);
      TEUCHOS_ASSERT(err == 0);
      TEUCHOS_ASSERT(regProlong[l+1][j]->RangeMap().PointSameAs(regRes[j]->Map()));
      TEUCHOS_ASSERT(regProlong[l+1][j]->DomainMap().PointSameAs(coarseRegB[j]->Map()));
    }
    sumInterfaceValues(coarseRegB, compMaps[l+1], maxRegPerProc,
        quasiRegRowMaps[l+1], regRowMaps[l+1], regRowImporters[l+1]);

    // Call V-cycle recursively
    vCycle(l+1, numLevels, maxFineIter, maxCoarseIter, omega, maxRegPerProc,
        coarseRegX, coarseRegB, regMatrices, regProlong, compMaps,
        quasiRegRowMaps, regRowMaps, regRowImporters, regInterfaceScalings);

    // Transfer coarse level correction to fine level
    std::vector<Teuchos::RCP<Epetra_Vector> > regCorrection(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      regCorrection[j] = Teuchos::rcp(new Epetra_Vector(*regRowMaps[l][j], true));
      int err = regProlong[l+1][j]->Multiply(false, *coarseRegX[j], *regCorrection[j]);
      TEUCHOS_ASSERT(err == 0);
      TEUCHOS_ASSERT(regProlong[l+1][j]->DomainMap().PointSameAs(coarseRegX[j]->Map()));
      TEUCHOS_ASSERT(regProlong[l+1][j]->RangeMap().PointSameAs(regCorrection[j]->Map()));
    }

    // apply coarse grid correction
    for (int j = 0; j < maxRegPerProc; j++) {
      int err = fineRegX[j]->Update(1.0, *regCorrection[j], 1.0);
      TEUCHOS_ASSERT(err == 0);
    }

    // post-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
        regInterfaceScalings[l], maxRegPerProc, compMaps[l], quasiRegRowMaps[l],
        regRowMaps[l], regRowImporters[l]);
  }
  else { // coarse level

    // coarse level solve
    /* For now, we do Jacobi. We need to figure out how to deal with free-floating
     * regions or recombine to a composite matrix on the coarse level
     */
    jacobiIterate(maxCoarseIter, omega, fineRegX, fineRegB, regMatrices[l],
        regInterfaceScalings[l], maxRegPerProc, compMaps[l],
        quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);
  }

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
  bool doing1D = false;
  int  globalNx, globalNy;

  myRank = Comm.MyPID();

  if (myRank == 0)
    std::cout << "Trying to sync with Matlab..." << std::endl
        << "If you're stuck here, check whether numProcs in Matlab and C-executable match" << std::endl;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;


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
  while ( command[iii] != ' ') iii++;
  sscanf(&(command[iii+1]),"%d",&globalNx);
  while ( command[iii] == ' ') iii++;
  while ( command[iii] != ' ') iii++;
  sscanf(&(command[iii+1]),"%d",&globalNy);
  while ( command[iii] == ' ') iii++;
  while ( command[iii] != ' ') iii++;
  if      (command[iii+1] == 'M') whichCase = MultipleRegionsPerProc;
  else if (command[iii+1] == 'R') whichCase = RegionsSpanProcs;
  else {fprintf(stderr,"%d: head messed up %s\n",myRank,command); exit(1);}
  sprintf(command,"/bin/rm -f myData_%d",myRank); system(command); sleep(1);

  // check for 1D or 2D problem
  if (globalNy == 1)
    doing1D = true;
  else
    doing1D = false;

  // ******************************************************************
  // Application Specific Data for LIDregion()
  // ******************************************************************
  struct widget appData;                            // ****************
  std::vector<int>  minGIDComp(maxRegPerProc);      // ****************
  std::vector<int>  maxGIDComp(maxRegPerProc);      // ****************
  // ******************************************************************
  // ******************************************************************


  std::vector<int> myRegions; // regions that myRank owns
  Teuchos::RCP<Epetra_CrsMatrix> AComp = Teuchos::null; // composite form of matrix
  Teuchos::RCP<Epetra_CrsMatrix> ACompSplit = Teuchos::null; // composite form of matrix
  Teuchos::RCP<Epetra_Map> mapComp = Teuchos::null; // composite map used to build AComp

  // regionsPerGID[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix row Map.
  Teuchos::RCP<Epetra_MultiVector> regionsPerGID = Teuchos::null;

  // regionsPerGIDWithGhosts[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix col Map.
  Teuchos::RCP<Epetra_MultiVector> regionsPerGIDWithGhosts = Teuchos::null;

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

  std::vector <Teuchos::RCP<Epetra_Map> > rowMapPerGrp(maxRegPerProc); // row map associated with myRank's ith region in composite layout
  std::vector <Teuchos::RCP<Epetra_Map> > colMapPerGrp(maxRegPerProc); // column map associated with myRank's ith region in composite layout
  std::vector <Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp(maxRegPerProc); // revised row map associated with myRank's ith region for regional layout
  std::vector <Teuchos::RCP<Epetra_Map> > revisedColMapPerGrp(maxRegPerProc); // revised column map associated with myRank's ith region for regional layout

  std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp(maxRegPerProc); // row importers per group
  std::vector<Teuchos::RCP<Epetra_Import> > colImportPerGrp(maxRegPerProc); // column importers per group
  std::vector<Teuchos::RCP<Epetra_Export> > rowExportPerGrp(maxRegPerProc); // row exporters per group

  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > quasiRegionGrpMats(maxRegPerProc); // region-wise matrices with quasiRegion maps (= composite GIDs)
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regionGrpMats(maxRegPerProc); // region-wise matrices in true region layout with unique GIDs for replicated interface DOFs

  Teuchos::RCP<Epetra_Map> coarseCompRowMap; ///< composite row map on the coarse grid
  std::vector<Teuchos::RCP<Epetra_Map> > coarseRowMapPerGrp(maxRegPerProc); // region-wise row map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Epetra_Map> > coarseQuasiRowMapPerGrp(maxRegPerProc); // region-wise row map in quasiRegion layout with original GIDs from fine level
  std::vector<Teuchos::RCP<Epetra_Map> > coarseColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Epetra_Map> > coarseAltColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Epetra_Import> > coarseRowImportPerGrp(maxRegPerProc); // coarse level row importer per group
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regionGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regionAltGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regCoarseMatPerGrp(maxRegPerProc); // coarse level operator 'RAP' in region layout

  std::vector<Teuchos::RCP<Epetra_Vector> > regionInterfaceScaling(maxRegPerProc);
  std::vector<Teuchos::RCP<Epetra_Vector> > coarseRegionInterfaceScaling(maxRegPerProc);

  std::vector<Teuchos::RCP<Epetra_Vector> > regNspViolation(maxRegPerProc); // violation of nullspace property in region layout

  Teuchos::RCP<Epetra_Vector> compX = Teuchos::null; // initial guess for truly composite calculations
  Teuchos::RCP<Epetra_Vector> compY = Teuchos::null; // result vector for truly composite calculations
  Teuchos::RCP<Epetra_Vector> regYComp = Teuchos::null; // result vector in composite layout, but computed via regional operations
  Teuchos::RCP<Epetra_Vector> nspViolation = Teuchos::null; // violation of nullspace property in composite layout

  std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegX(maxRegPerProc); // initial guess associated with myRank's ith region in quasiRegional layout
  std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegY(maxRegPerProc); // result vector associated with myRank's ith region in quasiRegional layout
  std::vector<Teuchos::RCP<Epetra_Vector> > regX(maxRegPerProc); // initial guess associated with myRank's ith region in regional layout
  std::vector<Teuchos::RCP<Epetra_Vector> > regY(maxRegPerProc); // result vector associated with myRank's ith region in regional layout

  std::vector<int> intIDs; // LIDs of interface DOFs
  std::vector<std::vector<int> > regIntIDs(maxRegPerProc); // LIDs of interface DOFs

  /* Stuff for multi-level algorithm
   *
   * To allow for multi-level schemes with more than two levels, we need to store
   * maps, matrices, vectors, and stuff like that on each level. Since we call the
   * multi-level scheme recursively, this should be reflected in the design of
   * variables.
   *
   * We use Teuchos::Array<T> to store each quantity on each level.
   */
  int numLevels = 0;
  Array<RCP<Epetra_Map> > compMaps;
  Array<std::vector<RCP<Epetra_Map> > > regRowMaps;
  Array<std::vector<RCP<Epetra_Map> > > quasiRegRowMaps;
  Array<std::vector<RCP<Epetra_CrsMatrix> > > regMatrices;
  Array<std::vector<RCP<Epetra_CrsMatrix> > > regProlong;
  Array<std::vector<RCP<Epetra_Import> > > regRowImporters;
  Array<std::vector<RCP<Epetra_Vector> > > regInterfaceScalings;

  Array<std::vector<std::vector<int> > > interfaceLIDs; // local IDs of interface nodes on each level in each group
  Array<std::vector<std::vector<std::vector<int> > > > interfaceGIDPairs; // pairs of GIDs of interface nodes on each level in each group

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

      mapComp = Teuchos::rcp(new Epetra_Map(-1,(int) fileData.size(),fileData.data(),0,Comm));
    }
    else if (strcmp(command,"PrintCompositeMap") == 0) {
       sleep(myRank);
       if (AComp == Teuchos::null) {
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
      Epetra_CrsMatrix* ACompPtr;
      EpetraExt::MatrixMarketFileToCrsMatrix("Amat.mm", *mapComp, ACompPtr);
      AComp = Teuchos::rcp(new Epetra_CrsMatrix(*ACompPtr));
    }
    else if (strcmp(command,"PrintCompositeMatrix") == 0) {
      AComp->Print(std::cout);
    }
    else if (strcmp(command,"LoadAndCommRegAssignments") == 0) {

      regionsPerGID = Teuchos::rcp(new Epetra_MultiVector(AComp->RowMap(),maxRegPerGID));

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
      regionsPerGIDWithGhosts = Teuchos::rcp(new Epetra_MultiVector(AComp->ColMap(),maxRegPerGID));
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
      if (doing1D) {
        int minGID,maxGID;
        for (int i = 0; i < (int) myRegions.size(); i++) {
          fscanf(fp,"%d%d",&minGID,&maxGID);
          minGIDComp[i] = minGID;
          maxGIDComp[i] = maxGID;
        }
        appData.gDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lowInd= (int *) malloc(sizeof(int)*3*myRegions.size());
        for (int i = 0; i < (int) myRegions.size(); i++) {
          fscanf(fp,"%d%d%d",&(appData.gDim[3*i]),&(appData.gDim[3*i+1]),&(appData.gDim[3*i+2]));
          fscanf(fp,"%d%d%d",&(appData.lDim[3*i]),&(appData.lDim[3*i+1]),&(appData.lDim[3*i+2]));
          fscanf(fp,"%d%d%d",&(appData.lowInd[3*i]),&(appData.lowInd[3*i+1]),&(appData.lowInd[3*i+2]));
        }
        appData.minGIDComp   = minGIDComp.data();
        appData.maxGIDComp   = maxGIDComp.data();
      }
      else {
        appData.nx = globalNx;
        appData.gDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lDimx= (int *) malloc(sizeof(int)*myRegions.size());
        appData.lDimy= (int *) malloc(sizeof(int)*myRegions.size());
        appData.trueCornerx= (int *) malloc(sizeof(int)*myRegions.size());
        appData.trueCornery= (int *) malloc(sizeof(int)*myRegions.size());
        appData.relcornerx= (int *) malloc(sizeof(int)*myRegions.size());
        appData.relcornery= (int *) malloc(sizeof(int)*myRegions.size());
        int garbage;
        for (int i = 0; i < (int) myRegions.size(); i++) {
          fscanf(fp,"%d%d%d",&(appData.gDim[3*i]),&(appData.gDim[3*i+1]),&(appData.gDim[3*i+2]));
          fscanf(fp,"%d%d%d",&(appData.lDimx[i]),&(appData.lDimy[i]),&garbage);
          fscanf(fp,"%d%d%d",&(appData.relcornerx[i]),&(appData.relcornery[i]),&garbage);
          fscanf(fp,"%d%d%d",&(appData.trueCornerx[i]),&(appData.trueCornery[i]),&garbage);
        }
      }
      appData.maxRegPerGID = maxRegPerGID;
      appData.myRegions    = myRegions.data();
      appData.colMap       = (Epetra_Map *) &(AComp->ColMap());
      appData.regionsPerGIDWithGhosts = regionsPerGIDWithGhosts;
      appData.myRank = myRank;
    }
    else if (strcmp(command,"MakeGrpRegColMaps") == 0) {
      if (whichCase == MultipleRegionsPerProc) {
        // clone rowMap
        for (int j=0; j < maxRegPerProc; j++) {
          if (j < (int) myRegions.size()) {
            colMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(-1,
                rowMapPerGrp[j]->NumMyElements(),
                rowMapPerGrp[j]->MyGlobalElements(), 0, Comm));
          }
          else colMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,Comm));
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
          if (doing1D) LID = LIDregion(&appData, i, 0);
          else LID = LID2Dregion(&appData, i, 0);
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
         colMapPerGrp[0] = Teuchos::rcp(new Epetra_Map(-1, (int) colIDsReg.size(),
             colIDsReg.data(), 0, Comm));
        }
        else colMapPerGrp[0] = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,Comm));
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
        int *colGIDsComp =  AComp->ColMap().MyGlobalElements();
        std::vector<int> tempRegIDs(nExtended);

        // must put revisedGIDs in application-provided order by
        // invoking LIDregion() and sorting
        for (int i = 0; i < nExtended; i++) {
          if (doing1D)
            tempRegIDs[i] = LIDregion(&appData, i, k);
          else
            tempRegIDs[i] = LID2Dregion(&appData, i, k);
        }
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

            // add entry to listDulicatedGIDs
          }
        }

        revisedRowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(-1,(int) revisedGIDs.size(),
                                            revisedGIDs.data(),0,Comm));
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
        revisedColMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(-1, (int) revisedGIDs.size(),
            revisedGIDs.data(), 0, Comm));
      }
      for (int k = (int) myRegions.size(); k < maxRegPerProc; k++) {
        revisedRowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,Comm));
        revisedColMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,Comm));
      }

      // Setup importers
      for (int j = 0; j < maxRegPerProc; j++) {
        rowImportPerGrp[j] = Teuchos::rcp(new Epetra_Import(*(rowMapPerGrp[j]), *mapComp));
        colImportPerGrp[j] = Teuchos::rcp(new Epetra_Import(*(colMapPerGrp[j]), *mapComp));
      }

    }
    else if (strcmp(command,"MakeGrpRegRowMaps") == 0) {
      std::vector<int> rowGIDsReg;
      int *colGIDsComp = AComp->ColMap().MyGlobalElements();
      for (int k = 0; k < (int) myRegions.size(); k++) {
        rowGIDsReg.resize(0);
        std::vector<int> tempRegIDs(AComp->ColMap().NumMyElements());
        for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {
          if (doing1D) tempRegIDs[i] = LIDregion(&appData, i, k);
          else tempRegIDs[i] = LID2Dregion(&appData, i, k);
        }

        std::vector<int> idx(tempRegIDs.size());
        std::iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(), [tempRegIDs](int i1,int i2) { return tempRegIDs[i1] < tempRegIDs[i2];});

        for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {
          if (tempRegIDs[idx[i]] != -1)
            rowGIDsReg.push_back(colGIDsComp[idx[i]]);
        }
        rowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(-1, (int) rowGIDsReg.size(),
            rowGIDsReg.data(), 0, Comm));
      }

      for (int k=(int) myRegions.size(); k < maxRegPerProc; k++) {
        rowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,Comm));
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
    else if (strcmp(command,"TestRegionalToComposite") == 0) {
      /* Create a random vector in regional layout and use both regionalToComposite
       * modes (with automatic and manual ADD of local values). Compare results.
       */

      Teuchos::RCP<Epetra_Vector> startVec = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      startVec->Random();

      std::vector<Teuchos::RCP<Epetra_Vector> > regStartVec(maxRegPerProc);
      std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegStartVec(maxRegPerProc);
      compositeToRegional(startVec, quasiRegStartVec, regStartVec, maxRegPerProc, rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

      // use automatic ADD of local values via Export()
      Teuchos::RCP<Epetra_Vector> autoCompVec = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      regionalToComposite(regStartVec, autoCompVec, maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Epetra_AddLocalAlso, false);

      Teuchos::RCP<Epetra_Vector> manCompVec = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      regionalToComposite(regStartVec, manCompVec, maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Add, true);

      // check for difference in result
      {
        Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
        int err = diff->Update(1.0, *autoCompVec, -1.0, *manCompVec, 0.0);
        TEUCHOS_ASSERT( err == 0);
        double diffNorm2 = 0.0;
        diff->Norm2(&diffNorm2);
        TEUCHOS_TEST_FOR_EXCEPT_MSG(diffNorm2 > 1.0e-15, "regionalToComposite() delivers "
            "different results for adding local values via Export() combine mode or manually.");
      }
    }
    else if (strcmp(command,"MakeQuasiRegionMatrices") == 0) {
      /* We use the edge-based splitting, i.e. we first modify off-diagonal
       * entries in the composite matrix, then decompose it into region matrices
       * and finally take care of diagonal entries by enforcing the nullspace
       * preservation constraint.
       */

      // copy and modify the composite matrix
      ACompSplit = Teuchos::rcp(new Epetra_CrsMatrix(*AComp));

      for (int row = 0; row < ACompSplit->NumMyRows(); ++row) { // loop over local rows of composite matrix
        int rowGID = ACompSplit->RowMap().GID(row);
        int numEntries; // number of nnz
        double* vals; // non-zeros in this row
        int* inds; // local column indices
        int err = ACompSplit->ExtractMyRowView(row, numEntries, vals, inds);
        TEUCHOS_ASSERT(err == 0);

        for (int c = 0; c < numEntries; ++c) { // loop over all entries in this row
          int col = inds[c];
          int colGID = ACompSplit->ColMap().GID(col);
          std::vector<int> commonRegions;
          if (rowGID != colGID) { // Skip the diagonal entry. It will be processed later.
            commonRegions = findCommonRegions(rowGID, colGID, *regionsPerGIDWithGhosts);
          }

          if (commonRegions.size() > 1) {
            vals[c] *= 0.5;
          }
        }
      }

      // Import data from ACompSplit into the quasiRegion matrices
      for (int j = 0; j < maxRegPerProc; j++) {
        quasiRegionGrpMats[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(rowMapPerGrp[j]),
            *(colMapPerGrp[j]), -1));
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
    else if (strcmp(command,"MakeInterfaceScalingFactorsRecursively") == 0) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!numLevels>0, "We require numLevel > 0. Probaly, numLevel has not been set, yet.");

      for (int l = 0; l < numLevels; l++)
      {
        // initialize region vector with all ones.
        for (int j = 0; j < maxRegPerProc; j++) {
          regInterfaceScalings[l][j] = rcp(new Epetra_Vector(*regRowMaps[l][j], true));
          regInterfaceScalings[l][j]->PutScalar(1.0);
        }

        // transform to composite layout while adding interface values via the Export() combine mode
        RCP<Epetra_Vector> compInterfaceScalingSum = rcp(new Epetra_Vector(*compMaps[l], true));
        regionalToComposite(regInterfaceScalings[l], compInterfaceScalingSum, maxRegPerProc, quasiRegRowMaps[l], regRowImporters[l], Add, true);

        /* transform composite layout back to regional layout. Now, GIDs associated
         * with region interface should carry a scaling factor (!= 1).
         */
        std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegInterfaceScaling(maxRegPerProc);
        compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
            regInterfaceScalings[l], maxRegPerProc, quasiRegRowMaps[l],
            regRowMaps[l], regRowImporters[l]);
      }
    }
    else if (strcmp(command,"MakeInterfaceScalingFactors") == 0) {
      // Fine level
      {
        // initialize region vector with all ones.
        for (int j = 0; j < maxRegPerProc; j++) {
          regionInterfaceScaling[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
          regionInterfaceScaling[j]->PutScalar(1.0);
        }

        // transform to composite layout while adding interface values via the Export() combine mode
        Teuchos::RCP<Epetra_Vector> compInterfaceScalingSum = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
        regionalToComposite(regionInterfaceScaling, compInterfaceScalingSum,
              maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Add, true);

        /* transform composite layout back to regional layout. Now, GIDs associated
         * with region interface should carry a scaling factor (!= 1).
         */
        std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegInterfaceScaling(maxRegPerProc);
        compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
            regionInterfaceScaling, maxRegPerProc, rowMapPerGrp,
            revisedRowMapPerGrp, rowImportPerGrp);
      }

      // Coarse level
      {
        // initialize region vector with all ones.
        for (int j = 0; j < maxRegPerProc; j++) {
          coarseRegionInterfaceScaling[j] = Teuchos::rcp(new Epetra_Vector(*coarseRowMapPerGrp[j], true));
          coarseRegionInterfaceScaling[j]->PutScalar(1.0);
        }

        // transform to composite layout while adding interface values via the Export() combine mode
        Teuchos::RCP<Epetra_Vector> compInterfaceScalingSum = Teuchos::rcp(new Epetra_Vector(*coarseCompRowMap, true));
        regionalToComposite(coarseRegionInterfaceScaling, compInterfaceScalingSum,
            maxRegPerProc, coarseQuasiRowMapPerGrp, coarseRowImportPerGrp, Add, true);

        /* transform composite layout back to regional layout. Now, GIDs associated
         * with region interface should carry a scaling factor (!= 1).
         */
        std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegInterfaceScaling(maxRegPerProc);
        compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
            coarseRegionInterfaceScaling, maxRegPerProc, coarseQuasiRowMapPerGrp,
            coarseRowMapPerGrp, coarseRowImportPerGrp);
      }
    }
    else if (strcmp(command,"MakeRegionMatrices") == 0) {

      // Copy data from quasiRegionGrpMats, but into new map layout
      {
        for (int j = 0; j < maxRegPerProc; j++) {
          // create empty matrix with correct row and column map
          regionGrpMats[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(revisedRowMapPerGrp[j]), *(revisedColMapPerGrp[j]), 3));

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

      // enforce nullspace constraint
      {
        // compute violation of nullspace property close to DBCs
        Teuchos::RCP<Epetra_Vector> nspVec = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
        nspVec->PutScalar(1.0);
        nspViolation = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
        int err = AComp->Apply(*nspVec, *nspViolation);
        TEUCHOS_ASSERT(err == 0);

        // move to regional layout
        std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegNspViolation(maxRegPerProc);
        createRegionalVector(quasiRegNspViolation, maxRegPerProc, rowMapPerGrp);
        createRegionalVector(regNspViolation, maxRegPerProc, revisedRowMapPerGrp);
        compositeToRegional(nspViolation, quasiRegNspViolation, regNspViolation,
            maxRegPerProc, rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

        /* The nullspace violation computed in the composite layout needs to be
         * transfered to the regional layout. Since we use this to compute
         * the splitting of the diagonal values, we need to split the nullspace
         * violation. We'd like to use the 'regInterfaceScaling', though this is
         * not setup at this point. So, let's compute it right now.
         *
         * ToDo: Move setup of 'regInterfaceScaling' up front to use it here.
         */
        {
          // initialize region vector with all ones.
          std::vector<Teuchos::RCP<Epetra_Vector> > interfaceScaling(maxRegPerProc);
          for (int j = 0; j < maxRegPerProc; j++) {
            interfaceScaling[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
            interfaceScaling[j]->PutScalar(1.0);
          }

          // transform to composite layout while adding interface values via the Export() combine mode
          Teuchos::RCP<Epetra_Vector> compInterfaceScalingSum = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
          regionalToComposite(interfaceScaling, compInterfaceScalingSum,
              maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Epetra_AddLocalAlso, false);

          /* transform composite layout back to regional layout. Now, GIDs associated
           * with region interface should carry a scaling factor (!= 1).
           */
          std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegInterfaceScaling(maxRegPerProc);
          compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
              interfaceScaling, maxRegPerProc, rowMapPerGrp,
              revisedRowMapPerGrp, rowImportPerGrp);

          // modify its interface entries
          for (int j = 0; j < maxRegPerProc; j++) {
            for (int i = 0; i < regNspViolation[j]->MyLength(); ++i) {
              (*regNspViolation[j])[i] /= (*interfaceScaling[j])[i];
            }
          }
        }
      }

      std::vector<Teuchos::RCP<Epetra_Vector> > regNsp(maxRegPerProc);
      std::vector<Teuchos::RCP<Epetra_Vector> > regCorrection(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) {
        regNsp[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
        regNsp[j]->PutScalar(1.0);

        regCorrection[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
        int err = regionGrpMats[j]->Apply(*regNsp[j], *regCorrection[j]);
        TEUCHOS_ASSERT(err == 0);
      }


      std::vector<Teuchos::RCP<Epetra_Vector> > regDiag(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) {
        regDiag[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
        int err = regionGrpMats[j]->ExtractDiagonalCopy(*regDiag[j]);
        TEUCHOS_ASSERT(err == 0);
        err = regDiag[j]->Update(-1.0, *regCorrection[j], 1.0, *regNspViolation[j], 1.0);
        TEUCHOS_ASSERT(err == 0);

        err = regionGrpMats[j]->ReplaceDiagonalValues(*regDiag[j]);
        TEUCHOS_ASSERT(err == 0);
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
      compX = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      compY = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      compX->Random();
//      compX->PutScalar(myRank);

      // perform matvec in composite layout
      {
        int err = AComp->Apply(*compX, *compY);
        TEUCHOS_ASSERT(err == 0);
      }

      // transform composite vector to regional layout
      compositeToRegional(compX, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);
      createRegionalVector(regY, maxRegPerProc, revisedRowMapPerGrp);

      // perform matvec in regional layout
      for (int j = 0; j < maxRegPerProc; j++) {
        int err = regionGrpMats[j]->Apply(*(regX[j]), *(regY[j]));
        TEUCHOS_ASSERT(err == 0);
      }

      // transform regY back to composite layout
      regYComp = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      regionalToComposite(regY, regYComp, maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Epetra_AddLocalAlso, false);

      Teuchos::RCP<Epetra_Vector> diffY = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      diffY->Update(1.0, *compY, -1.0, *regYComp, 0.0);

//      diffY->Print(std::cout);

      sleep(8);

      double diffNormTwo = 0.0;
      double diffNormInf = 0.0;
      diffY->Norm2(&diffNormTwo);
      diffY->NormInf(&diffNormInf);
      std::cout << myRank << ": diffNormTwo: " << diffNormTwo << "\tdiffNormInf: " << diffNormInf << std::endl;

    }
    else if (strcmp(command,"MakeMueLuTransferOperators") == 0) {

//      printRegionalMap("rowMapPerGrp", rowMapPerGrp, myRank);
//
//      sleep(3);
//
//      printRegionalMap("revisedRowMapPerGrp", revisedRowMapPerGrp, myRank);

      typedef MueLu::HierarchyManager<double,int,int,Xpetra::EpetraNode> HierarchyManager;
      typedef MueLu::Hierarchy<double,int,int,Xpetra::EpetraNode> Hierarchy;
      typedef MueLu::ParameterListInterpreter<double,int,int,Xpetra::EpetraNode> ParameterListInterpreter;
      typedef MueLu::Aggregates<int,int,Xpetra::EpetraNode> Aggregates;
      typedef MueLu::FactoryManagerBase FactoryManagerBase;
      typedef MueLu::Utilities<double,int,int,Xpetra::EpetraNode> Utilities;
      typedef Xpetra::Map<int,int,Xpetra::EpetraNode> Map;
      typedef Xpetra::MultiVector<double,int,int,Xpetra::EpetraNode> MultiVector;
      typedef Xpetra::MultiVectorFactory<double,int,int,Xpetra::EpetraNode> MultiVectorFactory;
      typedef Xpetra::Matrix<double,int,int,Xpetra::EpetraNode> Matrix;

      using Teuchos::RCP; using Teuchos::rcp; using Teuchos::ParameterList;

      Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));

      // convert Epetra region matrices to Xpetra
      std::vector<RCP<Matrix> > regionGrpXMats(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) {
        regionGrpXMats[j] = MueLu::EpetraCrs_To_XpetraMatrix<double,int,int,Xpetra::EpetraNode>(regionGrpMats[j]);
      }

      std::vector<RCP<Hierarchy> > regGrpHierarchy(maxRegPerProc);
      std::vector<RCP<HierarchyManager> > regMueLuFactory(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) {

        // create nullspace vector
        RCP<const Map> map = regionGrpXMats[j]->getRowMap();
        RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
        nullspace->putScalar(1.0);

        // create dummy coordinates vector
        RCP<MultiVector> coordinates = MultiVectorFactory::Build(map, 3);
        coordinates->putScalar(1.0);

//        const std::string xmlFileName = "/home/mmayr/wdir/codes/mayrmt-trilinos/Trilinos/packages/muelu/research/mmayr/composite_to_regions/src/muelu_expert.xml";
//        const std::string xmlFileName = "/home/mmayr/wdir/codes/mayrmt-trilinos/Trilinos/packages/muelu/research/mmayr/composite_to_regions/src/muelu_simple.xml";
        const std::string xmlFileName = "/home/mmayr/wdir/codes/mayrmt-trilinos/Trilinos/packages/muelu/research/mmayr/composite_to_regions/src/muelu_structured.xml";
        RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::DefaultComm<int>::getComm();
        Teuchos::ParameterList paramList;
        Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *teuchosComm);

        regMueLuFactory[j] = rcp(new ParameterListInterpreter(paramList));
        regGrpHierarchy[j] = regMueLuFactory[j]->CreateHierarchy();
        regGrpHierarchy[j]->GetLevel(0)->Set("A", regionGrpXMats[j]);
        regGrpHierarchy[j]->GetLevel(0)->Set("Nullspace", nullspace);
        regGrpHierarchy[j]->GetLevel(0)->Set("Coordinates", coordinates);

        // For structured aggregation
        {
          Array<int> gNodesPerDim(3);
          for (int i = 0; i < gNodesPerDim.size(); ++i)
            gNodesPerDim[i] = -1;

          Array<int> lNodesPerDim(3);
          // One-dimensional problems
          {
//            lNodesPerDim[0] = revisedRowMapPerGrp[j]->NumMyElements();
//            lNodesPerDim[1] = 1;
//            lNodesPerDim[2] = 1;
          }

          // Two-dimensional problems
          {
//            // caseFifteen
//            {
//              lNodesPerDim[0] = 4;
//              lNodesPerDim[1] = 4;
//              lNodesPerDim[2] = 1;
//            }

//            // caseSixteen
//            {
//              lNodesPerDim[0] = 7;
//              lNodesPerDim[1] = 7;
//              lNodesPerDim[2] = 1;
//            }

//            // caseSeventeen
//            {
//              lNodesPerDim[0] = 31;
//              lNodesPerDim[1] = 31;
//              lNodesPerDim[2] = 1;
//            }

//            // caseEightteen / caseNineteen
//            {
//              lNodesPerDim[0] = 16;
//              lNodesPerDim[1] = 16;
//              lNodesPerDim[2] = 1;
//            }

            // caseTwenty
            {
              lNodesPerDim[0] = 10;
              lNodesPerDim[1] = 10;
              lNodesPerDim[2] = 1;
            }
          }

          regGrpHierarchy[j]->GetLevel(0)->Set("gNodesPerDim", gNodesPerDim);
          regGrpHierarchy[j]->GetLevel(0)->Set("lNodesPerDim", lNodesPerDim);
        }

        regMueLuFactory[j]->SetupHierarchy(*regGrpHierarchy[j]);

      }

      // resize Arrays and vectors
      {
        // resize level containers
        numLevels = regGrpHierarchy[0]->GetNumLevels();
        compMaps.resize(numLevels);
        regRowMaps.resize(numLevels);
        quasiRegRowMaps.resize(numLevels);
        regMatrices.resize(numLevels);
        regProlong.resize(numLevels);
        regRowImporters.resize(numLevels);
        regInterfaceScalings.resize(numLevels);

        // resize group containers on each level
        for (int l = 0; l < numLevels; ++l) {
          regRowMaps[l].resize(maxRegPerProc);
          quasiRegRowMaps[l].resize(maxRegPerProc);
          regMatrices[l].resize(maxRegPerProc);
          regProlong[l].resize(maxRegPerProc);
          regRowImporters[l].resize(maxRegPerProc);
          regInterfaceScalings[l].resize(maxRegPerProc);
        }
      }

      // Fill fine level with our data
      {
        compMaps[0] = mapComp;
        quasiRegRowMaps[0] = rowMapPerGrp;
        regRowMaps[0] = revisedRowMapPerGrp;
        regRowImporters[0] = rowImportPerGrp;
        regMatrices[0] = regionGrpMats;

        /* MueLu stores prolongator on coarse level, so there is no prolongator
         * on the fine level. To have level containers of the same size, let's
         * just put in dummy data
         */
        std::vector<RCP<Epetra_CrsMatrix> > fineLevelProlong(maxRegPerProc);
        for (int j = 0; j < maxRegPerProc; ++j)
          fineLevelProlong[j] = Teuchos::null;
        regProlong[0] = fineLevelProlong;
      }

      // Get coarse level matrices and prolongators from MueLu hierarchy (Note: fine level has been dealt with previously)
      for (int l = 1; l < numLevels; ++l) { // Note: we start at level 1 (which is the first coarse level)
        for (int j = 0; j < maxRegPerProc; ++j) {
          RCP<MueLu::Level> level = regGrpHierarchy[j]->GetLevel(l);

          RCP<Matrix> regPXpetra = level->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
          regProlong[l][j] = Utilities::Op2NonConstEpetraCrs(regPXpetra);

          RCP<Matrix> regRAPXpetra = level->Get<RCP<Matrix> >("A", MueLu::NoFactory::get());
          regMatrices[l][j] = Utilities::Op2NonConstEpetraCrs(regRAPXpetra);

          regRowMaps[l][j] = Teuchos::rcp(new Epetra_Map(regMatrices[l][j]->RowMap())); // ToDo (mayr.mt) Do not copy!
        }
      }

      /* Reconstruct coarse-level maps
       *
       * We know the regional map on the coarse levels since they are just the
       * row maps of the coarse level operators. Though, we need to re-construct
       * the quasiRegional and composite maps ourselves.
       */

      {
        /* General strategy to create coarse level maps:
         * =============================================
         *
         * We need three maps on the next coarser level:
         * - regional map: just extract the RowMap() from the prolongator
         * - quasiRegional map: needs to be formed manually
         * - composite map: needs to be formed manually
         *
         * After extracting the RowMap() from the prolongator to be used as the
         * coarse level's regional map, we need to identify pairs/sets of GIDs
         * that represent the same (duplicated) DOF at a region interface in
         * order to form the quasiRegional and composite map for the next
         * coarser level.
         *
         * Identifying this pairs consists of several steps:
         * 1. Identify GIDs that represent the same (duplicated) DOF at a region
         *   interface on the current level
         * 2. Identify GIDs that represent the same (duplicated) DOF at a region
         *   interface on the next coarser level
         *
         * To form the quasiRegional map, we loop over the regional map and
         * - copy each GID that's associated with an interior DOF away from a
         *   region interface.
         * - replace each GID that's associated with an interface DOF by its
         *   quasiRegional counterpart.
         *
         * Note: We choose the quasiRegional (and composite) GID of an interface
         * DOF to be the std::min of all GIDs that represent that same
         * duplicated DOF.
         *
         * To form the composite map, we loop over the quasiRegional map and
         * - every proc accepts every GID that it owns.
         * - every proc ignores every GID that it doesn't own.
         */

        int err = 0;

        ////////////////////////////////////////////////////////////////////////
        // IDENTIFY DUPLICATED GIDS ON FINE LEVEL
        ////////////////////////////////////////////////////////////////////////

        RCP<Epetra_CrsMatrix> summedDuplicateMatrix = Teuchos::null;

        // Find lists of duplicated GIDs locally
        interfaceLIDs.resize(numLevels);
        interfaceGIDPairs.resize(numLevels);
        for (int l = 0; l < numLevels; ++l) {
          interfaceLIDs[l].resize(maxRegPerProc);
          interfaceGIDPairs[l].resize(maxRegPerProc);
        }

        for (int l = 0; l < numLevels - 1; ++l) {
          Comm.Barrier();
          if (Comm.MyPID() == 0) {
            std::cout << std::endl << std::endl
                << "Processing GID pairs on level " << l
                << std::endl << std::endl;
          }
          Comm.Barrier();

//          sleep(1);
//          std::cout << "Prolongator" << std::endl;
//          regProlong[l+1][0]->Print(std::cout);

          // create list of LIDs per group
          for (int j = 0; j < maxRegPerProc; j++) {
            for (int i = 0; i < regRowMaps[l][j]->NumMyElements(); ++i) {
              if (regRowMaps[l][j]->GID(i) != quasiRegRowMaps[l][j]->GID(i)) {
                // This is an interface node
                interfaceLIDs[l][j].push_back(i);
              }
            }
          }

//          for (int j = 0; j < maxRegPerProc; j++) {
//            std::cout << myRank << " | group = " << j << " | LIDs = ";
//            for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
//              std::cout << interfaceLIDs[l][j][i] << ", ";
//            std::cout << std::endl;
//          }

          for (int j = 0; j < maxRegPerProc; j++)
            interfaceGIDPairs[l][j].resize(interfaceLIDs[l][j].size());

          // create list of LIDPairs per group
          for (int j = 0; j < maxRegPerProc; j++)
            for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
              interfaceGIDPairs[l][j][i].push_back(quasiRegRowMaps[l][j]->GID(interfaceLIDs[l][j][i]));
          for (int j = 0; j < maxRegPerProc; j++)
            for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
              if (regRowMaps[l][j]->GID(interfaceLIDs[l][j][i]) != quasiRegRowMaps[l][j]->GID(interfaceLIDs[l][j][i]))
                interfaceGIDPairs[l][j][i].push_back(regRowMaps[l][j]->GID(interfaceLIDs[l][j][i]));

  //          std::cout << myRank << " | Print interfaceGIDPairs:" << std::endl;
  //          for (int j = 0; j < maxRegPerProc; j++) {
  //            std::cout << myRank << " | " << "Level " << l <<" | Group " << j << ":" << std::endl;
  //            for (std::size_t i = 0; i < interfaceGIDPairs[l][j].size(); ++i) {
  //              std::cout << "   " << myRank << " | ";
  //              for (std::size_t k = 0; k < interfaceGIDPairs[l][j][i].size(); ++k)
  //                std::cout << interfaceGIDPairs[l][j][i][k] << ", ";
  //              std::cout << std::endl;
  //            }
  //          }

          std::vector<RCP<Epetra_Vector> > regDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
          std::vector<RCP<Epetra_Vector> > quasiRegDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
          for (int j = 0; j < maxRegPerProc; ++j) {
            regDupGIDVec[j] = rcp(new Epetra_Vector(*regRowMaps[l][j], true));
            quasiRegDupGIDVec[j] = rcp(new Epetra_Vector(*quasiRegRowMaps[l][j], true));

            for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
              (*regDupGIDVec[j])[interfaceLIDs[l][j][i]] = 1.0;
          }

          RCP<Epetra_Vector> compDupGIDVec = rcp(new Epetra_Vector(*compMaps[l]), true);
          regionalToComposite(regDupGIDVec, compDupGIDVec, maxRegPerProc,
              quasiRegRowMaps[l], regRowImporters[l], Add, true);

          compositeToRegional(compDupGIDVec, quasiRegDupGIDVec, regDupGIDVec,
              maxRegPerProc, quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

          // create row/range/domain map for fine level duplicates mapping matrix
          RCP<Epetra_Map> duplicateMap = Teuchos::null;
          {
            std::vector<int> myIntGIDs;
            for (int j = 0; j < maxRegPerProc; j++) {
              for (std::size_t i = 0; i < interfaceGIDPairs[l][j].size(); ++i) {
                for (std::size_t k = 0; k < interfaceGIDPairs[l][j][i].size(); ++k) {
                  if (regRowMaps[l][j]->MyGID(interfaceGIDPairs[l][j][i][k]))
                    myIntGIDs.push_back(interfaceGIDPairs[l][j][i][k]);
                }
              }
            }

            duplicateMap = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myIntGIDs.size(), myIntGIDs.data(), 0, Comm));

  //            sleep(1);
  //            Comm.Barrier();
  //            std::cout << myRank << " | duplicateMap:" << std::endl;
  //            duplicateMap->Print(std::cout);
          }

          // create row/range/domain map for the transpose of the fine level duplicates mapping matrix
          RCP<Epetra_Map> fullDuplicateMap = Teuchos::null;
          {
            std::vector<int> myIntGIDs;
            for (int i = 0; i < regRowMaps[l][0]->NumMyElements(); ++i) {
              if ((*regDupGIDVec[0])[i] != 0)
                myIntGIDs.push_back(regRowMaps[l][0]->GID(i));
            }

  //              std::cout << myRank << " | myIntGIDs = ";
  //              for (auto gid : myIntGIDs)
  //                std::cout << gid << ", ";
  //              std::cout << std::endl;

            fullDuplicateMap = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myIntGIDs.size(), myIntGIDs.data(), 0, Comm));

  //            sleep(1);
  //            Comm.Barrier();
  //            std::cout << myRank << " | fullDuplicateMap:" << std::endl;
  //            fullDuplicateMap->Print(std::cout);
          }

          // Create and fill matrix
          RCP<Epetra_CrsMatrix> duplicateMatrix = rcp(new Epetra_CrsMatrix(Copy, *fullDuplicateMap, 2));
          {
            // Fill diagonal
            for (int i = 0; i < fullDuplicateMap->NumMyElements(); ++i) {
              int globalRow = fullDuplicateMap->GID(i);
              double vals[1];
              int inds[1];
              vals[0] = globalRow;
              inds[0] = globalRow;

              err = duplicateMatrix->InsertGlobalValues(globalRow, 1, &*vals, &*inds);
              TEUCHOS_ASSERT(err >= 0);
            }

            // Fill off-diagonals (if known)
            for (int i = 0; i < duplicateMap->NumMyElements(); ++i) {
              int globalRow = duplicateMap->GID(i);
              double vals[2];
              int inds[2];
              for (std::size_t k = 0; k < interfaceGIDPairs[l][0][i].size(); ++k) {
                vals[k] = globalRow;
                inds[k] = interfaceGIDPairs[l][0][i][k];
              }

              err = duplicateMatrix->InsertGlobalValues(globalRow, 2, &*vals, &*inds);
              TEUCHOS_ASSERT(err >= 0);
            }
          }

          err = duplicateMatrix->FillComplete();
          TEUCHOS_ASSERT(err == 0);

  //          sleep(1);
  //          Comm.Barrier();
  //          std::cout << myRank << " | duplicateMatrix:" << std::endl;
  //          duplicateMatrix->Print(std::cout);
  //
  //          sleep(1);
  //          Comm.Barrier();
  //          std::cout << myRank << " | duplicateMatrix->RangeMap():" << std::endl;
  //          duplicateMatrix->OperatorRangeMap().Print(std::cout);
  //
  //          sleep(1);
  //          Comm.Barrier();
  //          std::cout << myRank << " | duplicateMatrix->DomainMap():" << std::endl;
  //          duplicateMatrix->OperatorDomainMap().Print(std::cout);
  //
  //          sleep(1);
  //          Comm.Barrier();
  //          std::cout << myRank << " | duplicateMatrix->ColMap():" << std::endl;
  //          duplicateMatrix->ColMap().Print(std::cout);

          {
            EpetraExt::RowMatrix_Transpose transposer;
            RCP<Epetra_CrsMatrix> transDuplicateMatrix =
                rcp(new Epetra_CrsMatrix(dynamic_cast<Epetra_CrsMatrix&>(transposer(*duplicateMatrix))));

            TEUCHOS_ASSERT(!transDuplicateMatrix.is_null());

            err = transDuplicateMatrix->FillComplete();
            TEUCHOS_ASSERT(err == 0);

    //          sleep(1);
    //          Comm.Barrier();
    //          std::cout << myRank << " | transDuplicateMatrix:" << std::endl;
    //          transDuplicateMatrix->Print(std::cout);


            summedDuplicateMatrix = rcp(new Epetra_CrsMatrix(Copy, *fullDuplicateMap, 2));
            Epetra_CrsMatrix* summedDuplicateMatrix_ptr = summedDuplicateMatrix.get();
            err = EpetraExt::MatrixMatrix::Add(*transDuplicateMatrix, false, 1.0, *duplicateMatrix, false, 1.0, summedDuplicateMatrix_ptr);
            TEUCHOS_ASSERT(err == 0);
          }

          err = summedDuplicateMatrix->FillComplete(*fullDuplicateMap, *fullDuplicateMap);
  //          err = summedDuplicateMatrix->FillComplete();
          TEUCHOS_ASSERT(err == 0);

  //          sleep(1);
  //          Comm.Barrier();
  //          std::cout << myRank << " | summedDuplicateMatrix:" << std::endl;
  //          summedDuplicateMatrix->Print(std::cout);

          std::vector<std::vector<int> > myDuplicates; // pairs of duplicates with locally owned GID listed first.
          myDuplicates.resize(summedDuplicateMatrix->NumMyRows());
          for (int i = 0; i < summedDuplicateMatrix->NumMyRows(); ++i) {
            int numEntries = 0;
            double* vals;
            int* inds;
  //          std::cout << myRank << " | Extracting my row " << i << ", global row " << summedDuplicateMatrix->RowMap().GID(i) << std::endl;
            err = summedDuplicateMatrix->ExtractMyRowView(i, numEntries, vals, inds);
            TEUCHOS_ASSERT(err == 0);

            std::vector<int> gidsToSort;
            for (int k = 0; k < numEntries; ++k) {
  //            std::cout << myRank << " | inds[" << k << "] = " << inds[k] << std::endl;
              gidsToSort.push_back(summedDuplicateMatrix->ColMap().GID(inds[k]));
            }

            /* Note: At intersections of more than two regions, i.e. at interior
             * vertices of the regional interface, the list of gidsToSort is not
             * identical on all procs. In particular, only the proc owning the
             * original GID of this vertex knows the full list of all duplicated
             * GIDs. Those procs, that own a duplicated GID, only know about the
             * two-member pair of their duplicate and the original GID, but not
             * about the duplicated GIDs of all other procs.
             *
             * However, this does not matter since we are only looking for the
             * mapping of each proc's duplicated GID to the original GID, which
             * every proc knows about.
             *
             * Bottom line, this should be sufficient.
             */

  //          sleep(1);
  //          std::cout << myRank << " | gidsToSort:" << std::endl << "  ";
  //          for (auto gid : gidsToSort)
  //             std::cout << gid << ", ";
  //           std::cout << std::endl;
  //           sleep(1);

            /* sort s.t. my GID is first
             *
             * 1. Find index of the one GID that I own
             * 2. Insert this GID at the beginning of the vector
             * 3. Erase this GID from its initial position in the vector
             *
             * ToDo (mayr.mt) Is this really necessary?
             */
            {
  //            int indOfMyGID = -1;
  //            for (std::size_t k = 0; k < gidsToSort.size(); ++k) {
  //              if (regRowMaps[l][0]->MyGID(gidsToSort[k])) {
  //                indOfMyGID = k;
  //                break;
  //              }
  //            }
  //            TEUCHOS_ASSERT(indOfMyGID >= 0);
  //
  //            int tmpIndex = gidsToSort[indOfMyGID];
  //            gidsToSort.erase(gidsToSort.begin() + indOfMyGID);
  //            gidsToSort.insert(gidsToSort.begin(), tmpIndex);

  //            for (std::size_t k = 0; i < gidsToSort.size(); ++k)
  //              myDuplicates[i].push_back(gidsToSort[i]);

              myDuplicates[i] = gidsToSort;
            }
          }

//          sleep(myRank);
//          std::cout << std::endl << std::endl << std::endl << std::endl;
//          for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
//            std::cout << myRank << " | myDuplicates:" << std::endl << "  ";
//            for (std::size_t k = 0; k < myDuplicates[i].size(); ++k)
//              std::cout << myDuplicates[i][k] << ", ";
//            std::cout << std::endl;
//          }
//          std::cout << std::endl << std::endl << std::endl << std::endl;

          ////////////////////////////////////////////////////////////////////////
          // IDENTIFY DUPLICATED GIDS ON COARSE LEVEL
          ////////////////////////////////////////////////////////////////////////

          std::vector<std::vector<int> > myCoarseDuplicates; // pairs of duplicates with locally owned GID listed first on coarse level.
          myCoarseDuplicates.resize(myDuplicates.size());

          std::vector<int> myCoarseInterfaceGIDs;
          for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
            const int rowGID = myDuplicates[i][0];
            const int rowLID = regProlong[l+1][0]->RowMap().LID(rowGID);
            int numEntries;
            double* vals;
            int* inds;
            err = regProlong[l+1][0]->ExtractMyRowView(rowLID, numEntries, vals, inds);
  //          std::cout << myRank << " | ExtractMyRowView err = " << err << std::endl;
            TEUCHOS_ASSERT(err == 0);
            TEUCHOS_ASSERT(numEntries == 1); // tentative P: only one entry per row

            myCoarseInterfaceGIDs.push_back(regProlong[l+1][0]->ColMap().GID(inds[0]));
          }

  //        std::cout << "myCoarseInterfaceGIDs on proc " << myRank << ": ";
  //        for (auto gid : myCoarseInterfaceGIDs) {
  //          std::cout << gid << ", ";
  //        }
  //        std::cout << std::endl;

          // Build row/range/domain map of duplicate mapping matrix
          RCP<Epetra_Map> interfaceMap = Teuchos::null;
          {
            std::vector<int> interfaceMapGIDs(myDuplicates.size());
            for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
              interfaceMapGIDs[i] = myDuplicates[i][0];
            }

  //          std::cout << "interfaceMapGIDs on proc " << myRank << ":      ";
  //          for (auto gid : interfaceMapGIDs) {
  //            std::cout << gid << ", ";
  //          }
  //          std::cout << std::endl;

            interfaceMap = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), interfaceMapGIDs.size(), interfaceMapGIDs.data(), 0, Comm));
          }

          RCP<Epetra_CrsMatrix> duplicateMapping = rcp(new Epetra_CrsMatrix(Copy, *interfaceMap, 2));

          // Fill the matrix
          RCP<Epetra_CrsMatrix> transDuplicateMapping = Teuchos::null;
          {
            int numRows = myDuplicates.size();
            std::vector<int> rowPtr(numRows);
            std::vector<std::vector<double> > vals(numRows);
            std::vector<std::vector<int> > colInds(numRows);

            for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
              rowPtr[rowIdx] = myDuplicates[rowIdx][0];
              for (std::size_t k = 0; k < myDuplicates[rowIdx].size(); ++k) {
                vals[rowIdx].push_back(myCoarseInterfaceGIDs[rowIdx]);
                colInds[rowIdx].push_back(myDuplicates[rowIdx][k]);
              }
            }

  //          for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
  //            std::cout << myRank << " | rowIdx = " << rowIdx << ", rowGID = " << rowPtr[rowIdx] << ":";
  //            for (std::size_t i = 0; i < vals[rowIdx].size(); ++i)
  //              std::cout << "(" << colInds[rowIdx][i] << "|" << vals[rowIdx][i] << "), ";
  //            std::cout << std::endl;
  //          }

            // local dummy insertion
            for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
              err = duplicateMapping->InsertGlobalValues(rowPtr[rowIdx], colInds[rowIdx].size(), vals[rowIdx].data(), colInds[rowIdx].data());
              TEUCHOS_ASSERT(err >= 0);
            }

            err = duplicateMapping->FillComplete();
            TEUCHOS_ASSERT(err == 0);

  //          sleep(1);
  //          Comm.Barrier();
  //          duplicateMapping->Print(std::cout);

            EpetraExt::RowMatrix_Transpose transposer;
            transDuplicateMapping = rcp(new Epetra_CrsMatrix(dynamic_cast<Epetra_CrsMatrix&>(transposer(*duplicateMapping))));

            TEUCHOS_ASSERT(!transDuplicateMapping.is_null());

            err = transDuplicateMapping->FillComplete();
            TEUCHOS_ASSERT(err == 0);

  //          sleep(1);
  //          std::cout << myRank << " | Printing the tranpose ..." << std::endl;
  //          Comm.Barrier();
  //          transDuplicateMapping->Print(std::cout);
          }

          sleep(1);

          /* Extract coarse level duplicates from transDuplicateMapping
           *
           * Note: In 2D, there will be duplicates which need to be removed.
           * One could maybe use something along the lines of:
           *    sort( vec.begin(), vec.end() );
           *    vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
           *
           */
          std::vector<std::vector<double> > myCoarseInterfaceDuplicates(transDuplicateMapping->NumMyRows());
          {
            for (int i = 0; i < transDuplicateMapping->NumMyRows(); ++i) {
              int numEntries;
              double* myVals;
              int* myInds;
              err = transDuplicateMapping->ExtractMyRowView(i, numEntries, myVals, myInds);
              TEUCHOS_ASSERT(err == 0);

              myCoarseInterfaceDuplicates[i].resize(numEntries);
              myCoarseInterfaceDuplicates[i].assign(myVals, myVals+numEntries);

//              std::cout << myRank << " | myCoarseInterfaceDuplicates[" << i << "] = ";
//              for (auto id : myCoarseInterfaceDuplicates[i])
//                std::cout << id << ", ";
//              std::cout << std::endl;
            }

          }

          ////////////////////////////////////////////////////////////////////////
          // CREATE COARSE LEVEL MAPS
          ////////////////////////////////////////////////////////////////////////
          {
            sleep(1);
            std::cout << myRank << " | Printing regRowMaps[" << l+1 << "][0] ..." << std::endl;
            Comm.Barrier();
            regRowMaps[l+1][0]->Print(std::cout);
            sleep(2);

            // create quasiRegional map
            {
              std::vector<int> myQuasiRegGIDs;

              for (int i = 0; i < regRowMaps[l+1][0]->NumMyElements(); ++i) {
                // grab current regional GID to be processed
                int currGID = regRowMaps[l+1][0]->GID(i);
                int quasiGID = currGID; // assign dummy value

                // find quasiRegional counterpart
                {
                  // Is this an interface GID?
                  bool isInterface = false;
                  for (std::size_t k = 0; k < myCoarseInterfaceDuplicates.size(); ++k) {
                    for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk) {
                      if (currGID == myCoarseInterfaceDuplicates[k][kk])
                        isInterface = true;
                    }
                  }

                  if (isInterface) {
                    for (std::size_t k = 0; k < myCoarseInterfaceDuplicates.size(); ++k) {
                      bool found = false;
                      for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk) {
                        if (currGID == myCoarseInterfaceDuplicates[k][kk])
                          found = true;
                      }

                      if (found) {
                        for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk)
                          quasiGID = std::min(quasiGID, Teuchos::as<int>(myCoarseInterfaceDuplicates[k][kk]));

  //                      std::cout << myRank << " | Interface GID " << currGID << " is replaced by quasiGID " << quasiGID << std::endl;
                        break;
                      }
                    }
                  }
                  else { // interior node --> take GID from regional map
                    quasiGID = currGID;
                  }
                }

                TEUCHOS_ASSERT(quasiGID>=0);
                myQuasiRegGIDs.push_back(quasiGID);
              }

              quasiRegRowMaps[l+1][0] = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myQuasiRegGIDs.size(), myQuasiRegGIDs.data(), 0, Comm));
              TEUCHOS_ASSERT(!quasiRegRowMaps[l+1][0].is_null());

              sleep(1);
              std::cout << myRank << " | Printing quasiRegRowMaps[" << l+1 << "][0] ..." << std::endl;
              Comm.Barrier();
              quasiRegRowMaps[l+1][0]->Print(std::cout);

            }

            // create composite map
            {
              std::vector<int> myCompGIDs;
              for (int j = 0; j < maxRegPerProc; j++) {
                for (int i = 0; i < quasiRegRowMaps[l+1][j]->NumMyElements(); ++i) {
                  const int trialGID = quasiRegRowMaps[l+1][j]->GID(i);

                  if (regRowMaps[l+1][j]->MyGID(trialGID))
                    myCompGIDs.push_back(trialGID);
                }
              }

              compMaps[l+1] = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myCompGIDs.size(), myCompGIDs.data(), 0, Comm));
              TEUCHOS_ASSERT(!compMaps[l+1].is_null());

              sleep(1);
              std::cout << myRank << " | Printing compMaps["<< l+1 << "] ..." << std::endl;
              Comm.Barrier();
              compMaps[l+1]->Print(std::cout);
            }

            // create regRowImporter
            for (int j = 0; j < maxRegPerProc; ++j) {
              regRowImporters[l+1][j] =
                  rcp(new Epetra_Import(*quasiRegRowMaps[l+1][j], *compMaps[l+1]));
              TEUCHOS_ASSERT(!regRowImporters[l+1][j].is_null());
            }
          }
        }
      }
    }
    else if (strcmp(command,"MakeRegionTransferOperators") == 0) {
      /* Populate a fine grid vector with 1 at coarse nodes and 0 at fine nodes.
       * Then, transform to regional layout and find GIDs with entry 1
       */
      int numNodes = mapComp->NumGlobalPoints() / 3 + 1;
      std::vector<double> vals(numNodes);
      std::vector<int> ind(numNodes);
      for (int i = 0; i < numNodes; ++i)
      {
        vals[i] = 1.0;
        ind[i] = 3*i;
      }
      Teuchos::RCP<Epetra_Vector> coarseGridToggle = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      coarseGridToggle->ReplaceGlobalValues(numNodes, vals.data(), ind.data());

      // compute coarse composite row map
      {
        std::vector<int> coarseGIDs;
        for (int i = 0; i < coarseGridToggle->MyLength(); ++i) {
          if ((*coarseGridToggle)[i] != 0)
            coarseGIDs.push_back(coarseGridToggle->Map().GID(i));
        }
        coarseCompRowMap = Teuchos::rcp(new Epetra_Map(-1, coarseGIDs.size(), coarseGIDs.data(), 0, Comm));
      }

      // transform composite vector to regional layout
      std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegCoarseGridToggle(maxRegPerProc);
      std::vector<Teuchos::RCP<Epetra_Vector> > regCoarseGridToggle(maxRegPerProc);
      compositeToRegional(coarseGridToggle, quasiRegCoarseGridToggle,
          regCoarseGridToggle, maxRegPerProc, rowMapPerGrp, revisedRowMapPerGrp,
          rowImportPerGrp);

      // create coarse grid row maps in region layout
      for (int j = 0; j < maxRegPerProc; j++) {
        std::vector<int> coarseRowGIDsReg;
        std::vector<int> coarseQuasiRowGIDsReg;
        for (int i = 0; i < regCoarseGridToggle[j]->Map().NumMyElements(); ++i) {
          if ((*regCoarseGridToggle[j])[i] == 1.0) {
            coarseQuasiRowGIDsReg.push_back(quasiRegCoarseGridToggle[j]->Map().GID(i));
            coarseRowGIDsReg.push_back(regCoarseGridToggle[j]->Map().GID(i));
          }
        }

        coarseQuasiRowMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(-1, coarseQuasiRowGIDsReg.size(), coarseQuasiRowGIDsReg.data(), 0, Comm));
        coarseRowMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(-1, coarseRowGIDsReg.size(), coarseRowGIDsReg.data(), 0, Comm));
      }

      // setup coarse level row importer
      for (int j = 0; j < maxRegPerProc; j++) {
        coarseRowImportPerGrp[j] = Teuchos::rcp(new Epetra_Import(*(coarseQuasiRowMapPerGrp[j]),
            *coarseCompRowMap));
      }

      // create coarse grid column map
      Teuchos::RCP<Epetra_Vector> compGIDVec = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
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
          TEUCHOS_ASSERT(false);
      }

      // transform composite vector to regional layout
      std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegGIDVec(maxRegPerProc);
      std::vector<Teuchos::RCP<Epetra_Vector> > regGIDVec(maxRegPerProc);
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

      for (int j = 0; j < maxRegPerProc; j++) {
        std::vector<int> regColGIDs;
        for (int i = 0; i < revisedRowMapPerGrp[j]->NumMyElements(); ++i) {
          const int currGID = revisedRowMapPerGrp[j]->GID(i);

          // loop over existing regColGIDs and see if current GID is already in that list
          bool found = false;
          for (std::size_t k = 0; k < regColGIDs.size(); ++k) {
            if ((*regGIDVec[j])[regGIDVec[j]->Map().LID(currGID)] == regColGIDs[k]) {
              found = true;
              break;
            }
          }

          if (not found)
            regColGIDs.push_back((*regGIDVec[j])[regGIDVec[j]->Map().LID(currGID)]);
        }

        coarseColMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(-1, regColGIDs.size(), regColGIDs.data(), 0, Comm));
      }

      for (int j = 0; j < maxRegPerProc; j++) {
        std::vector<int> regAltColGIDs;
        int *colGIDs = regionGrpMats[j]->ColMap().MyGlobalElements();
        if (regionGrpMats[j]->ColMap().NumMyElements() > 0) {

          // check if leftmost point is next to a cpt to the left in which case 1st ghost point is a cpt
          if ( (appData.lowInd[3*j]%3) == 1 ) {
           int NIOwn = regionGrpMats[j]->RowMap().NumMyElements();
           regAltColGIDs.push_back(colGIDs[NIOwn]);
          }
          int start = 3 - (appData.lowInd[3*j]%3);
          if (start == 3) start = 0;
          for (int i = start; i < appData.lDim[3*j] ; i += 3)
            regAltColGIDs.push_back(colGIDs[i]);

          // check if pt to the right of rightmost point is is cpt, i.e. last ghost is a cpt
          if ( ((appData.lowInd[3*j]+appData.lDim[3*j])%3) == 0 ) {
            int NLast = regionGrpMats[j]->ColMap().NumMyElements()-1;
            regAltColGIDs.push_back(colGIDs[NLast]);
          }
        }
        coarseAltColMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(-1, regAltColGIDs.size(), regAltColGIDs.data(), 0, Comm));
      }

      // Build the actual prolongator
      for (int j = 0; j < maxRegPerProc; j++) {
        regionGrpProlong[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *revisedRowMapPerGrp[j], *coarseColMapPerGrp[j], 1, false));
        for (int c = 0; c < regionGrpProlong[j]->NumMyCols(); ++c) {
          for (int r = 0; r < regionGrpProlong[j]->NumMyRows(); ++r) {
            if (regionGrpProlong[j]->ColMap().GID(c) == (*regGIDVec[j])[r]) {
              double myVals[1];
              int inds[1];
              myVals[0] = 1.0; //1.0/sqrt(3.0); //1.0; // use all ones of the prolongator
              inds[0] = c;
              regionGrpProlong[j]->InsertMyValues(r, 1, &*myVals, &*inds);
            }
          }
        }
        int err = regionGrpProlong[j]->FillComplete(*coarseRowMapPerGrp[j], *revisedRowMapPerGrp[j]);
        TEUCHOS_ASSERT(err == 0);
//        regionGrpProlong[j]->Print(std::cout);
      }

      for (int j = 0; j < maxRegPerProc; j++) {
        double myVals[1];
        int inds[1];
        myVals[0] = 1.0/3.0;
        regionAltGrpProlong[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *revisedRowMapPerGrp[j], *coarseAltColMapPerGrp[j], 1, false));
//        int *coarseCol = coarseAltColMapPerGrp[j]->MyGlobalElements();
        int NccSize    = coarseAltColMapPerGrp[j]->NumMyElements();
//        int *fineRow   = revisedRowMapPerGrp[j]->MyGlobalElements();
        int NfrSize    = revisedRowMapPerGrp[j]->NumMyElements();
        if (NfrSize > 0) {
          int fstart = 3 - (appData.lowInd[3*j]%3);
          if (fstart == 3) fstart = 0;
          int cstart = 0;
          if (fstart == 2) cstart = 1;

          // need to add 1st prolongator row as this is not addressed by loop below
          if (cstart == 1) {
            inds[0] = 0;
            regionAltGrpProlong[j]->InsertMyValues(0, 1, &*myVals, &*inds);
          }
          int i;
          for (i = fstart; i < appData.lDim[3*j] ; i += 3) {
            inds[0] = cstart;
            regionAltGrpProlong[j]->InsertMyValues(i, 1, &*myVals, &*inds);
            if (i > 0)         regionAltGrpProlong[j]->InsertMyValues(i-1, 1, &*myVals, &*inds);
            if (i < NfrSize-1) regionAltGrpProlong[j]->InsertMyValues(i+1, 1, &*myVals, &*inds);
            cstart++;
          }
          // last cpoint hasn't been addressed because someone else owns it
          if (cstart < NccSize) {
            inds[0] = cstart;
            regionAltGrpProlong[j]->InsertMyValues(NfrSize-1, 1, &*myVals, &*inds);
          }
        }
        regionAltGrpProlong[j]->FillComplete(*coarseRowMapPerGrp[j], *revisedRowMapPerGrp[j]); // ToDo: test this
//        regionAltGrpProlong[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"PrintRegProlongatorPerGrp") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        regionGrpProlong[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"MakeCoarseLevelOperator") == 0) {

      std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regAP(maxRegPerProc); // store intermediate result A*P
      for (int j = 0; j < maxRegPerProc; j++) {
        // Compute A*P
        {
          regAP[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *revisedRowMapPerGrp[j], 3, false));

          int err = EpetraExt::MatrixMatrix::Multiply(*regionGrpMats[j], false, *regionGrpProlong[j], false, *regAP[j], false);
          TEUCHOS_ASSERT(err == 0);
          err = regAP[j]->FillComplete();
          TEUCHOS_ASSERT(err == 0);
        }

        // Compute R*(AP) = P'*(AP)
        {
          regCoarseMatPerGrp[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *coarseRowMapPerGrp[j], 3, false));

          int err = EpetraExt::MatrixMatrix::Multiply(*regionGrpProlong[j], true, *regAP[j], false, *regCoarseMatPerGrp[j], false);
          TEUCHOS_ASSERT(err == 0);
          err = regCoarseMatPerGrp[j]->FillComplete();
          TEUCHOS_ASSERT(err == 0);
        }
      }

//      for (int j = 0; j < maxRegPerProc; j++) {
//          regAP[j]->Print(std::cout);
//      }
    }
    else if (strcmp(command,"PrintRegCoarseMatPerGrp") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        regCoarseMatPerGrp[j]->Print(std::cout);
      }
    }
    else if (strcmp(command,"RunTwoLevelMethod") == 0) {
      // initial guess for solution
      compX = Teuchos::rcp(new Epetra_Vector(*mapComp, true));

// debugging using shadow.m
//double *z;
//compX->ExtractView(&z); for (int kk = 0; kk < compX->MyLength(); kk++) z[kk] = (double) kk*kk;
//compX->Print(std::cout);
      // forcing vector
      Teuchos::RCP<Epetra_Vector> compB = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      {
        compB->ReplaceGlobalValue(compB->GlobalLength() - 1, 0, 1.0e-3);
      }
//      {
//      compB->PutScalar(1.0);
//      compB->ReplaceGlobalValue(0, 0, 0.0);
//      }
//      {
//        compB->ReplaceGlobalValue(15, 0, 1.0);
//      }
//      {
//        compB->ReplaceGlobalValue(16, 0, 1.0);
//      }

      // residual vector
      Teuchos::RCP<Epetra_Vector> compRes = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      {
        int err = AComp->Multiply(false, *compX, *compRes);
        TEUCHOS_ASSERT(err == 0);
        err = compRes->Update(1.0, *compB, -1.0);
        TEUCHOS_ASSERT(err == 0);
      }

      // transform composite vectors to regional layout
      compositeToRegional(compX, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);

      std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegB(maxRegPerProc);
      std::vector<Teuchos::RCP<Epetra_Vector> > regB(maxRegPerProc);
      compositeToRegional(compB, quasiRegB, regB, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);

      std::vector<Teuchos::RCP<Epetra_Vector> > regRes(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) { // step 1
        regRes[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
      }

      // define max iteration counts
      const int maxVCycle = 2;
      const int maxFineIter = 10;
      const int maxCoarseIter = 100;
      const double omega = 0.67;

      // Richardson iterations
      for (int cycle = 0; cycle < maxVCycle; ++cycle) {

        // a single 2-level V-Cycle
        {
          // -----------------------------------------------------------------------
          // pre-smoothing on fine level
          // -----------------------------------------------------------------------
          jacobiIterate(maxFineIter, omega, regX, regB, regionGrpMats,
              regionInterfaceScaling, maxRegPerProc, mapComp, rowMapPerGrp,
              revisedRowMapPerGrp, rowImportPerGrp);

          computeResidual(regRes, regX, regB, regionGrpMats, mapComp,
              rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

          // -----------------------------------------------------------------------
          // Transfer to coarse level
          // -----------------------------------------------------------------------
          std::vector<Teuchos::RCP<Epetra_Vector> > coarseRegX(maxRegPerProc);
          std::vector<Teuchos::RCP<Epetra_Vector> > coarseRegB(maxRegPerProc);
          for (int j = 0; j < maxRegPerProc; j++) {
            coarseRegX[j] = Teuchos::rcp(new Epetra_Vector(regionGrpProlong[j]->DomainMap(), true));
            coarseRegB[j] = Teuchos::rcp(new Epetra_Vector(regionGrpProlong[j]->DomainMap(), true));

            for (int i = 0; i < regRes[j]->MyLength(); ++i)
              (*regRes[j])[i] /= (*regionInterfaceScaling[j])[i];

            int err = regionGrpProlong[j]->Multiply(true, *regRes[j], *coarseRegB[j]);
            TEUCHOS_ASSERT(err == 0);
            TEUCHOS_ASSERT(regionGrpProlong[j]->RangeMap().PointSameAs(regRes[j]->Map()));
            TEUCHOS_ASSERT(regionGrpProlong[j]->DomainMap().PointSameAs(coarseRegB[j]->Map()));
          }

          sumInterfaceValues(coarseRegB, coarseCompRowMap, maxRegPerProc,
              coarseQuasiRowMapPerGrp, coarseRowMapPerGrp,
              coarseRowImportPerGrp);

          // -----------------------------------------------------------------------
          // Perform region-wise Jacobi on coarse level
          // -----------------------------------------------------------------------
          jacobiIterate(maxCoarseIter, omega, coarseRegX, coarseRegB, regCoarseMatPerGrp,
              coarseRegionInterfaceScaling, maxRegPerProc, coarseCompRowMap, coarseQuasiRowMapPerGrp,
              coarseRowMapPerGrp, coarseRowImportPerGrp);

          // -----------------------------------------------------------------------
          // Transfer to fine level
          // -----------------------------------------------------------------------
          std::vector<Teuchos::RCP<Epetra_Vector> > regCorrection(maxRegPerProc);
          for (int j = 0; j < maxRegPerProc; j++) {
            regCorrection[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
            int err = regionGrpProlong[j]->Multiply(false, *coarseRegX[j], *regCorrection[j]);
            TEUCHOS_ASSERT(err == 0);
            TEUCHOS_ASSERT(regionGrpProlong[j]->DomainMap().PointSameAs(coarseRegX[j]->Map()));
            TEUCHOS_ASSERT(regionGrpProlong[j]->RangeMap().PointSameAs(regCorrection[j]->Map()));
          }

          // apply coarse grid correction
          for (int j = 0; j < maxRegPerProc; j++) {
            int err = regX[j]->Update(1.0, *regCorrection[j], 1.0);
            TEUCHOS_ASSERT(err == 0);
          }

          // -----------------------------------------------------------------------
          // post-smoothing on fine level
          // -----------------------------------------------------------------------
          jacobiIterate(maxFineIter, omega, regX, regB, regionGrpMats,
              regionInterfaceScaling, maxRegPerProc, mapComp, rowMapPerGrp,
              revisedRowMapPerGrp, rowImportPerGrp);
        }

        // check for convergence
        {
          computeResidual(regRes, regX, regB, regionGrpMats, mapComp,
              rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

          compRes = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
          regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
              rowImportPerGrp, Add, true);
          double normRes = 0.0;
          compRes->Norm2(&normRes);

          if (normRes < 1.0e-12)
            break;
        }
      }

      // -----------------------------------------------------------------------
      // Print fine-level solution
      // -----------------------------------------------------------------------

      compX->Comm().Barrier();
      sleep(1);

      regionalToComposite(regX, compX, maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Zero, false);
      std::cout << "compX after V-cycle" << std::endl;
      compX->Print(std::cout);
      sleep(2);
    }
    else if (strcmp(command,"RunVCycle") == 0) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!numLevels>0, "We require numLevel > 0. Probaly, numLevel has not been set, yet.");

      /* We first use the non-level container variables to setup the fine grid problem.
       * This is ok since the initial setup just mimics the application and the outer
       * Krylov method.
       *
       * We switch to using the level container variables as soon as we enter the
       * recursive part of the algorithm.
       */

      // initial guess for solution
      compX = Teuchos::rcp(new Epetra_Vector(*mapComp, true));

// debugging using shadow.m
//double *z;
//compX->ExtractView(&z); for (int kk = 0; kk < compX->MyLength(); kk++) z[kk] = (double) kk*kk;
//compX->Print(std::cout);
      // forcing vector
      Teuchos::RCP<Epetra_Vector> compB = Teuchos::rcp(new Epetra_Vector(*mapComp, true));

      // 1D
      {
        {
          compB->ReplaceGlobalValue(compB->GlobalLength() - 1, 0, 1.0e-3);
        }
//        {
//        compB->PutScalar(1.0);
//        compB->ReplaceGlobalValue(0, 0, 0.0);
//        }
//        {
//          compB->ReplaceGlobalValue(15, 0, 1.0);
//        }
//        {
//          compB->ReplaceGlobalValue(16, 0, 1.0);
//        }
      }

      std::vector<int> dirichletGIDs;

      // 2D
      {
        const int nx = 19; // global number of nodes in x-direction
        for (int i = 0; i < nx; ++i)
          dirichletGIDs.push_back(i);
        for (int i = 0; i < nx; ++i) {
          dirichletGIDs.push_back(i*nx);
          dirichletGIDs.push_back((i+1)*nx - 1);
        }
        for (int i = 0; i < nx; ++i)
          dirichletGIDs.push_back((nx*(nx-1) + i));
      }

//      for (auto gid : dirichletGIDs)
//        std::cout << gid << ", ";
//      std::cout << std::endl;

      {
        compB->PutScalar(1.0e-3);
        for (std::size_t i = 0; i < dirichletGIDs.size(); ++i)
          compB->ReplaceGlobalValue(dirichletGIDs[i], 0, 0.0);
      }

      // residual vector
      Teuchos::RCP<Epetra_Vector> compRes = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
      {
        int err = AComp->Multiply(false, *compX, *compRes);
        TEUCHOS_ASSERT(err == 0);
        err = compRes->Update(1.0, *compB, -1.0);
        TEUCHOS_ASSERT(err == 0);
      }

      // transform composite vectors to regional layout
      compositeToRegional(compX, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);

      std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegB(maxRegPerProc);
      std::vector<Teuchos::RCP<Epetra_Vector> > regB(maxRegPerProc);
      compositeToRegional(compB, quasiRegB, regB, maxRegPerProc, rowMapPerGrp,
          revisedRowMapPerGrp, rowImportPerGrp);

      std::vector<Teuchos::RCP<Epetra_Vector> > regRes(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) { // step 1
        regRes[j] = Teuchos::rcp(new Epetra_Vector(*revisedRowMapPerGrp[j], true));
      }

      /////////////////////////////////////////////////////////////////////////
      // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
      /////////////////////////////////////////////////////////////////////////

      // define max iteration counts
      const int maxVCycle = 200;
      const int maxFineIter = 20;
      const int maxCoarseIter = 100;
      const double omega = 0.67;

      // Richardson iterations
      for (int cycle = 0; cycle < maxVCycle; ++cycle) {

        vCycle(0, numLevels, maxFineIter, maxCoarseIter, omega, maxRegPerProc, regX, regB, regMatrices,
            regProlong, compMaps, quasiRegRowMaps, regRowMaps, regRowImporters,
            regInterfaceScalings);

        ////////////////////////////////////////////////////////////////////////
        // SWITCH BACK TO NON-LEVEL VARIABLES
        ////////////////////////////////////////////////////////////////////////

        // check for convergence
        {
          computeResidual(regRes, regX, regB, regionGrpMats, mapComp,
              rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

          compRes = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
          regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
              rowImportPerGrp, Add, true);
          double normRes = 0.0;
          compRes->Norm2(&normRes);

          if (normRes < 1.0e-12)
            break;
        }
      }

      // -----------------------------------------------------------------------
      // Print fine-level solution
      // -----------------------------------------------------------------------

      compX->Comm().Barrier();
      sleep(1);

      regionalToComposite(regX, compX, maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Zero, false);

      std::cout << myRank << " | compX after V-cycle" << std::endl;
      sleep(1);
      compX->Print(std::cout);
      sleep(2);

      // Write solution to file for printing
      std::string outFileName = "compX.mm";
      Xpetra::IO<double,int,int,Xpetra::EpetraNode>::Write(outFileName, *Xpetra::toXpetra<int,Xpetra::EpetraNode>(compX));
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
    else if (strcmp(command,"PrintRegVectorInterfaceScaling") == 0) {
      sleep(myRank);
      for (int j = 0; j < maxRegPerProc; j++) {
        printf("%d: region vector regionInterfaceScaling Grp %d\n", myRank, j);
        regionInterfaceScaling[j]->Print(std::cout);
      }
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
   Teuchos::RCP<Epetra_MultiVector> regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

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
int LID2Dregion(void *ptr, int LIDcomp, int whichGrp)
{
   struct widget * myWidget = (struct widget *) ptr;

//   int        *minGIDComp  = myWidget->minGIDComp;
//   int        *maxGIDComp  = myWidget->maxGIDComp;
   int        *myRegions   = myWidget->myRegions;
   Epetra_Map *colMap      = myWidget->colMap;
   int        maxRegPerGID = myWidget->maxRegPerGID;
   int       *trueCornerx  = myWidget->trueCornerx; // global coords of region
   int       *trueCornery  = myWidget->trueCornery; // corner within entire 2D mesh
                                                    // across all regions/procs
   int       *relcornerx   = myWidget->relcornerx;  // coords of corner relative
   int       *relcornery   = myWidget->relcornery;  // to region corner
   int       *lDimx        = myWidget->lDimx;
   int       *lDimy        = myWidget->lDimy;
   Teuchos::RCP<Epetra_MultiVector> regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

   int curRegion = myRegions[whichGrp];
   if (LIDcomp >=  colMap->NumMyElements()) return(-1);

   double *jthRegions;
   const int *colGIDsComp = colMap->MyGlobalElements();

   int xGIDComp =  colGIDsComp[LIDcomp]%(myWidget->nx);
   int yGIDComp = (colGIDsComp[LIDcomp] - xGIDComp)/myWidget->nx;

   if (xGIDComp < trueCornerx[whichGrp]+relcornerx[whichGrp]) return(-1);
   if (yGIDComp < trueCornery[whichGrp]+relcornery[whichGrp]) return(-1);
   if (xGIDComp > trueCornerx[whichGrp]+relcornerx[whichGrp]+lDimx[whichGrp]-1) return(-1);
   if (yGIDComp > trueCornery[whichGrp]+relcornery[whichGrp]+lDimy[whichGrp]-1) return(-1);


   bool found = false;
   for (int j = 0; j <  maxRegPerGID; j++) {
      jthRegions = (*regionsPerGIDWithGhosts)[j];
      if  ( ((int) jthRegions[LIDcomp]) == curRegion ) {
         found = true;
         break;
      }
   }
   if (found == false) return(-1);

   return(
    (yGIDComp - relcornery[whichGrp]-trueCornery[whichGrp])*lDimx[whichGrp]+ (xGIDComp - relcornerx[whichGrp]-trueCornerx[whichGrp]));
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
void printGrpMaps(std::vector <Teuchos::RCP<Epetra_Map> > &mapPerGrp, int maxRegPerProc, char *str)
{
   for (int j = 0; j < maxRegPerProc; j++) {
      printf("%s in grp(%d):",str,j);
      Teuchos::RCP<Epetra_Map> jthMap = mapPerGrp[j];
      const int *GIDsReg = jthMap->MyGlobalElements();
      for (int i = 0; i < mapPerGrp[j]->NumMyElements(); i++) {
         printf("%d ",GIDsReg[i]);
      }
      printf("\n");
   }
}
