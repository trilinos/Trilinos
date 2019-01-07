// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <numeric>

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#define RegionsSpanProcs  1
#define MultipleRegionsPerProc  2
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <Kokkos_DefaultNode.hpp>

#include <MueLu.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

/* This is a driver that is not included in any other file.
 * So, we should be fine to create useful typedefs for Xpetra here and use them in the entire file.
 */
typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#include "Xpetra_UseShortNames.hpp"

LocalOrdinal LIDregion(void *ptr, int LIDcomp, int whichGrp);
LocalOrdinal LID2Dregion(void *ptr, int LIDcomp, int whichGrp);

extern void edgeGhosts(int ArowPtr[], int Acols[], int &nGhostFound, int ghostCompLIDs[], int edgeLength, int alongX, int ownedEdge, int interiorEdge, int start, int ownedX, int ownedY);

extern void fillCircleSquareData(int ArowPtr[], int Acols[], int ownedX, int ownedY, int Cx, int Cy, int Rx, int Ry, std::vector<int> &appData);

extern LocalOrdinal LIDregionCircleSquare(void *ptr, int compLID,int whichGrp);

// Input data is read into a generic vector.
// Use these enums to access entries in this vector.
enum InputDataIndices
{
  inpData_isStructured,
  inpData_ownedX,
  inpData_ownedY,
  inpData_regionX,
  inpData_regionY,
  inpData_cornerX,
  inpData_cornerY,
  inpData_nGhosts,
  inpData_firstLIDsOfGhosts
};

// this little widget handles application specific data
// used to implement LIDregion()
struct widget {
   int *minGIDComp;
   int *maxGIDComp;
   int *myRegions;
   Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> *colMap;
   int maxRegPerGID;
   Teuchos::RCP<MultiVector> regionsPerGIDWithGhosts;
   int *gDim, *lDim, *lowInd;
   int *trueCornerx; // global coordinates of region
   int *trueCornery; // corner within entire 2D mesh
   int *relcornerx;  // coordinates of corner relative
   int *relcornery;  // to region corner
   int *lDimx;
   int *lDimy;
   int nx;
   int myRank;
};

//! Print an object in regional layout to screen
template <class T>
void printRegionalObject(const std::string objName, ///< string to be used for screen output
    const std::vector<Teuchos::RCP<T> > regObj, ///< regional object to be printed to screen
    const int myRank, ///< rank of calling proc
    Teuchos::FancyOStream& outstream ///< output stream
    )
{
  for (int j = 0; j < (int) regObj.size(); j++) {
    outstream << myRank << ": " << objName << " " << j << std::endl;
    regObj[j]->describe(outstream, Teuchos::VERB_EXTREME);
  }
}

/*! \brief Transform composite vector to regional layout
 *
 *  Starting from a vector in composite layout, we
 *  1. import it into an auxiliary vector in the quasiRegional layout
 *  2. replace the quasiRegional map of the auxiliary vector with the regional map
 */
void compositeToRegional(Teuchos::RCP<const Vector> compVec, ///< Vector in composite layout [in]
    std::vector<Teuchos::RCP<Vector> >& quasiRegVecs, ///< Vector in quasiRegional layout [in/out]
    std::vector<Teuchos::RCP<Vector> >& regVecs, ///< Vector in regional layout [in/out]
    const int maxRegPerProc, ///< max number of regions per proc [in]
    const std::vector<Teuchos::RCP<Map> > rowMapPerGrp, ///< row maps in region layout [in]
    const std::vector<Teuchos::RCP<Map> > revisedRowMapPerGrp, ///< revised row maps in region layout [in]
    const std::vector<Teuchos::RCP<Import> > rowImportPerGrp ///< row importer in region layout [in]
    )
{
  // quasiRegional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create empty vectors and fill it by extracting data from composite vector
    quasiRegVecs[j] = VectorFactory::Build(rowMapPerGrp[j], true);
    TEUCHOS_ASSERT(!quasiRegVecs[j].is_null());
    quasiRegVecs[j]->doImport(*compVec, *(rowImportPerGrp[j]), Xpetra::INSERT);
  }

  // regional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create regVecs vector (copy from quasiRegVecs and swap the map)
    regVecs[j] = quasiRegVecs[j]; // assignment operator= does deep copy in Xpetra
    TEUCHOS_ASSERT(!regVecs[j].is_null());
    regVecs[j]->replaceMap(revisedRowMapPerGrp[j]);
  }

  return;
}

/*! \brief Transform regional vector to composite layout
 *
 *  Starting from a \c Epetra_Vector in regional layout, we
 *  1. replace the regional map with the quasiRegional map
 *  2. export it into a vector with composite layout using the \c CombineMode \c Add.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 */
void regionalToComposite(const std::vector<Teuchos::RCP<Vector> >& regVec, ///< Vector in region layout [in]
    Teuchos::RCP<Vector> compVec, ///< Vector in composite layout [in/out]
    const int maxRegPerProc, ///< max number of regions per proc
    const std::vector<Teuchos::RCP<Map> > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
    const std::vector<Teuchos::RCP<Import> > rowImportPerGrp, ///< row importer in region layout [in]
    const Xpetra::CombineMode combineMode ///< Combine mode for import/export [in]
    )
{
  /* Let's fake an ADD combine mode that also adds local values by
   * 1. exporting quasiRegional vectors to auxiliary composite vectors (1 per group)
   * 2. add all auxiliary vectors together
   */

  Teuchos::RCP<MultiVector> partialCompVec = MultiVectorFactory::Build(compVec->getMap(), maxRegPerProc, true);
  TEUCHOS_ASSERT(!partialCompVec.is_null());

  std::vector<Teuchos::RCP<Vector> > quasiRegVec(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    // copy vector and replace map
    quasiRegVec[j] = regVec[j]; // Xpetra: operator= mimics copy constructor
    TEUCHOS_ASSERT(!quasiRegVec[j].is_null());

    quasiRegVec[j]->replaceMap(rowMapPerGrp[j]);

    // ToDo (mayr.mt) Use input variable 'combineMode'
    partialCompVec->getVectorNonConst(j)->doExport(*quasiRegVec[j], *(rowImportPerGrp[j]), Xpetra::ADD);
  }

  compVec->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
  for (int j = 0; j < maxRegPerProc; j++) {
    compVec->update(Teuchos::ScalarTraits<Scalar>::one(), *partialCompVec->getVector(j), Teuchos::ScalarTraits<Scalar>::one());
  }

  return;
}

/*! \brief Transform regional matrix to composite layout
 *
 *  Starting from a \c Matrix in regional layout, we
 *  1. copy data from regional matrix into quasiRegional matrix and set all maps
 *     to be quasiRegional maps
 *  2. export it into a \c Matrix with composite layout using the given \c combineMode.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 *
 *  \return Composite matrix that is fill-completed
 */
void regionalToComposite(const std::vector<Teuchos::RCP<Matrix> >& regMat, ///< Matrix in region layout [in]
    Teuchos::RCP<Matrix> compMat, ///< Matrix in composite layout [in/out]
    const int maxRegPerProc, ///< max number of regions per proc
    const std::vector<Teuchos::RCP<Map> > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
    const std::vector<Teuchos::RCP<Map> > colMapPerGrp, ///< col maps in quasiRegion layout [in]
    const std::vector<Teuchos::RCP<Import> > rowImportPerGrp, ///< row importer in region layout [in]
    const Xpetra::CombineMode combineMode ///< Combine mode for import/export [in]
    )
{
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::size_t;

  const Scalar scalarZero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar scalarOne = Teuchos::ScalarTraits<Scalar>::one();

  // Make sure we add into an zero composite matrix
  compMat->setAllToScalar(scalarZero);

  /* Let's fake an ADD combine mode that also adds local values by
     * 1. exporting quasiRegional matrices to auxiliary composite matrices (1 per group)
     * 2. add all auxiliary matrices together
     */

  // Copy data from regMat into quasiRegMat
  std::vector<RCP<Matrix> > quasiRegMat(maxRegPerProc);
  {
    for (int j = 0; j < maxRegPerProc; j++) {
      quasiRegMat[j] = rcp(new CrsMatrixWrap(rowMapPerGrp[j], colMapPerGrp[j], 9, Xpetra::DynamicProfile));

      // Extract current quasi-region CrsMatrix
      RCP<CrsMatrix> quasiRegionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegMat[j])->getCrsMatrix();

      // Extract current region CrsMatrix
      RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regMat[j])->getCrsMatrix();

      // Pull out the data from the region CrsMatrix
      Teuchos::ArrayRCP<const size_t> rowptrRegion;
      Teuchos::ArrayRCP<const LocalOrdinal> colindRegion;
      Teuchos::ArrayRCP<const Scalar> valuesRegion;
      regionCrsMat->getAllValues(rowptrRegion, colindRegion, valuesRegion);

      // Do a deep copy of values
      // (at least we've been doing deep copies so far, maybe we could do shallow copies to save time?)
      Teuchos::ArrayRCP<size_t> rowptrQuasiRegion(rowptrRegion.size());
      Teuchos::ArrayRCP<LocalOrdinal> colindQuasiRegion(colindRegion.size());
      Teuchos::ArrayRCP<Scalar> valuesQuasiRegion(valuesRegion.size());

      quasiRegionCrsMat->allocateAllValues(valuesQuasiRegion.size(), rowptrQuasiRegion, colindQuasiRegion, valuesQuasiRegion);

      for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(rowptrQuasiRegion.size()); ++idx) {
        rowptrQuasiRegion[idx] = rowptrRegion[idx];
      }

      for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(colindQuasiRegion.size()); ++idx) {
        colindQuasiRegion[idx] = colindRegion[idx];
        valuesQuasiRegion[idx] = valuesRegion[idx];
      }

      // Set and fillComplete the quasiRegion CrsMatrix
      quasiRegionCrsMat->setAllValues(rowptrQuasiRegion, colindQuasiRegion, valuesQuasiRegion);
      quasiRegionCrsMat->expertStaticFillComplete(rowMapPerGrp[j], rowMapPerGrp[j]);
    }
  }

  // Export from quasiRegional format to composite layout
  std::vector<RCP<Matrix> > partialCompMat(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    partialCompMat[j] = MatrixFactory::Build(compMat->getRowMap(), 3, Xpetra::DynamicProfile);
    partialCompMat[j]->doExport(*(quasiRegMat[j]), *(rowImportPerGrp[j]), Xpetra::INSERT);
    partialCompMat[j]->fillComplete(compMat->getDomainMap(), compMat->getRangeMap());
  }

  // Add all partialCompMat together
  for (int j = 0; j < maxRegPerProc; j++) {
    MatrixMatrix::TwoMatrixAdd(*partialCompMat[j], false, scalarOne, *compMat, scalarOne);
  }

  compMat->fillComplete();

  return;
}

/*! \brief Sum region interface values
 *
 *  Sum values of interface GIDs using the underlying Export() routines. Technically, we perform the
 *  exchange/summation of interface data by exporting a regional vector to the composite layout and
 *  then immediately importing it back to the regional layout. The Export() involved when going to the
 *  composite layout takes care of the summation of interface values.
 */
void sumInterfaceValues(std::vector<Teuchos::RCP<Vector> >& regVec,
    Teuchos::RCP<const Map> compMap,
    const int maxRegPerProc, ///< max number of regions per proc [in]
    const std::vector<Teuchos::RCP<Map> > rowMapPerGrp,///< row maps in region layout [in]
    const std::vector<Teuchos::RCP<Map> > revisedRowMapPerGrp,///< revised row maps in region layout [in]
    const std::vector<Teuchos::RCP<Import> > rowImportPerGrp ///< row importer in region layout [in])
    )
{
  Teuchos::RCP<Vector> compVec = VectorFactory::Build(compMap, true);
  TEUCHOS_ASSERT(!compVec.is_null());

  std::vector<Teuchos::RCP<Vector> > quasiRegVec(maxRegPerProc);
  regionalToComposite(regVec, compVec, maxRegPerProc, rowMapPerGrp,
      rowImportPerGrp, Xpetra::ADD);

  compositeToRegional(compVec, quasiRegVec, regVec, maxRegPerProc,
      rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

  return;
}

//! Create an empty vector in the regional layout
void createRegionalVector(std::vector<Teuchos::RCP<Vector> >& regVecs, ///< regional vector to be filled
    const int maxRegPerProc, ///< max number of regions per process
    const std::vector<Teuchos::RCP<Map> > revisedRowMapPerGrp ///< regional map
    )
{
  regVecs.resize(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++)
    regVecs[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);

  return;
}

/*! \brief Compute the residual \f$r = b - Ax\f$
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute y = A*x in regional layout.
 *  2. Sum interface values of y to account for duplication of interface DOFs.
 *  3. Compute r = b - y
 */
std::vector<Teuchos::RCP<Vector> > computeResidual(
    std::vector<Teuchos::RCP<Vector> >& regRes, ///< residual (to be evaluated)
    const std::vector<Teuchos::RCP<Vector> > regX, ///< left-hand side (solution)
    const std::vector<Teuchos::RCP<Vector> > regB, ///< right-hand side (forcing term)
    const std::vector<Teuchos::RCP<Matrix> > regionGrpMats,
    Teuchos::RCP<const Map> mapComp, ///< composite map, computed by removing GIDs > numDofs in revisedRowMapPerGrp
    const std::vector<Teuchos::RCP<Map> > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
    const std::vector<Teuchos::RCP<Map> > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
    const std::vector<Teuchos::RCP<Import> > rowImportPerGrp ///< row importer in region layout [in]
    )
{
  const int maxRegPerProc = regX.size();

  /* Update the residual vector
   * 1. Compute tmp = A * regX in each region
   * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
   * 3. Compute r = B - tmp
   */
  for (int j = 0; j < maxRegPerProc; j++) { // step 1
    regionGrpMats[j]->apply(*regX[j], *regRes[j]);
//    TEUCHOS_ASSERT(regionGrpMats[j]->getDomainMap()->isSameAs(*regX[j]->getMap()));
//    TEUCHOS_ASSERT(regionGrpMats[j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
  }

  sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
      revisedRowMapPerGrp, rowImportPerGrp);

  for (int j = 0; j < maxRegPerProc; j++) { // step 3
    regRes[j]->update(1.0, *regB[j], -1.0);
//    TEUCHOS_ASSERT(regRes[j]->getMap()->isSameAs(*regB[j]->getMap()));
  }

  return regRes;
}

/*! \brief Do Jacobi smoothing
 *
 *  Perform Jacobi smoothing in the region layout using the true diagonal value
 *  recovered from the splitted matrix.
 */
void jacobiIterate(const int maxIter,
    const double omega,
    std::vector<Teuchos::RCP<Vector> >& regX, // left-hand side (or solution)
    const std::vector<Teuchos::RCP<Vector> > regB, // right-hand side (or residual)
    const std::vector<Teuchos::RCP<Matrix> > regionGrpMats, // matrices in true region layout
    const std::vector<Teuchos::RCP<Vector> > regionInterfaceScaling, // recreate on coarse grid by import Add on region vector of ones
    const int maxRegPerProc, ///< max number of regions per proc [in]
    Teuchos::RCP<const Map> mapComp, ///< composite map
    const std::vector<Teuchos::RCP<Map> > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
    const std::vector<Teuchos::RCP<Map> > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
    const std::vector<Teuchos::RCP<Import> > rowImportPerGrp ///< row importer in region layout [in]
    )
{
  using Teuchos::RCP;

  const Scalar scalarZero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar scalarOne = Teuchos::ScalarTraits<Scalar>::one();

  std::vector<RCP<Vector> > regRes(maxRegPerProc);
  createRegionalVector(regRes, maxRegPerProc, revisedRowMapPerGrp);

  // extract diagonal from region matrices, recover true diagonal values, invert diagonal
  std::vector<RCP<Vector> > diag(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    // extract inverse of diagonal from matrix
    diag[j] = VectorFactory::Build(regionGrpMats[j]->getRowMap(), true);
    regionGrpMats[j]->getLocalDiagCopy(*diag[j]);
    diag[j]->elementWiseMultiply(scalarOne, *diag[j], *regionInterfaceScaling[j], scalarZero); // ToDo Does it work to pass in diag[j], but also return into the same variable?
    diag[j]->reciprocal(*diag[j]);
  }

  for (int iter = 0; iter < maxIter; ++iter) {

    /* Update the residual vector
     * 1. Compute tmp = A * regX in each region
     * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
     * 3. Compute r = B - tmp
     */
    for (int j = 0; j < maxRegPerProc; j++) { // step 1

//      Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//      regRes[j]->getMap()->describe(*fos, Teuchos::VERB_EXTREME);
//      regionGrpMats[j]->getRangeMap()->describe(*fos, Teuchos::VERB_EXTREME);

//      TEUCHOS_ASSERT(regionGrpMats[j]->getDomainMap()->isSameAs(*regX[j]->getMap()));
//      TEUCHOS_ASSERT(regionGrpMats[j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));

      regionGrpMats[j]->apply(*regX[j], *regRes[j]);
    }

    sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
        revisedRowMapPerGrp, rowImportPerGrp);

    for (int j = 0; j < maxRegPerProc; j++) { // step 3
      regRes[j]->update(1.0, *regB[j], -1.0);
    }

    // check for convergence
    {
      RCP<Vector> compRes = VectorFactory::Build(mapComp, true);
      regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
          rowImportPerGrp, Xpetra::ADD);
      Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();

      if (normRes < 1.0e-12)
        return;
    }

    for (int j = 0; j < maxRegPerProc; j++) {
      // update solution according to Jacobi's method
      regX[j]->elementWiseMultiply(omega, *diag[j], *regRes[j], scalarOne);
    }
  }

  return;
}

/*! \brief Find common regions of two nodes
 *
 */
std::vector<int> findCommonRegions(const GlobalOrdinal nodeA, ///< GID of first node
    const GlobalOrdinal nodeB, ///< GID of second node
    const MultiVector& nodesToRegions ///< mapping of nodes to regions
    )
{
  // extract node-to-regions mapping for both nodes A and B
  Teuchos::Array<int> regionsA;
  Teuchos::Array<int> regionsB;
  {
    Teuchos::RCP<const Map> map = nodesToRegions.getMap();
    for (std::size_t i = 0; i < nodesToRegions.getNumVectors(); ++i) {
      regionsA.push_back(nodesToRegions.getData(i)[map->getLocalElement(nodeA)]);
      regionsB.push_back(nodesToRegions.getData(i)[map->getLocalElement(nodeB)]);
    }
  }

//  // Print list of regions for both nodes
//  {
//    int myRank = nodesToRegions.getMap()->getComm()->getRank();
//    Teuchos::RCP<const Map> map = nodesToRegions.getMap();
//    std::cout << myRank << ": nodeA = " << map->getGlobalElement(nodeA) << ": ";
//    for (std::size_t i = 0; i < nodesToRegions.getNumVectors(); ++i)
//      std::cout << ", " << regionsA[i];
//    std::cout << std::endl;
//    std::cout << myRank << ": nodeB = " << map->getGlobalElement(nodeB) << ": ";
//      for (std::size_t i = 0; i < nodesToRegions.getNumVectors(); ++i)
//        std::cout << ", " << regionsB[i];
//      std::cout << std::endl;
//  }

  // identify common regions
  std::vector<int> commonRegions(nodesToRegions.getNumVectors());
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
//    int myRank = nodesToRegions.getMap()->getComm()->getRank();
//    std::cout << myRank << ": " << nodeA << "/" << nodeB << ": ";
//    for (int i = 0; i < finalCommonRegions.size(); ++i)
//      std::cout << ", " << finalCommonRegions[i];
//    std::cout << std::endl;
//  }

  return finalCommonRegions;
}

/*! \brief Create coarse level maps with continuous GIDs
 *
 *  The direct solver requires maps with continuous GIDs. Starting from the
 *  coarse level composite maps with discontinuous GIDs, we create a new row map
 *  and a matching column map.
 *
 *  Range and Domain map happen to correspond to the Row map, so we don't have
 *  to deal with them in particular.
 */
void createContinuousCoarseLevelMaps(Teuchos::RCP<const Map> rowMap, ///< row map
    Teuchos::RCP<const Map> colMap, ///< column map
    Teuchos::RCP<Map>& contRowMap, ///< row map with continuous GIDs
    Teuchos::RCP<Map>& contColMap ///< column map with continuous GIDs
    )
{
  // Create row map with continuous GIDs
  contRowMap = MapFactory::Build(rowMap->lib(), Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
      rowMap->getNodeNumElements(), Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), rowMap->getComm());

  /* Create column map based on row map with continuous GIDs
   *
   * We use an Importer to create an auxiliary vector in colMap format containing
   * the GIDs of the contRowMap as its entries. By looping over its LIDs, we can
   * then form the contColMap.
   */
  Teuchos::RCP<Import> rowColImport = ImportFactory::Build(rowMap, colMap);
  Teuchos::RCP<Vector> colGIDVec = VectorFactory::Build(rowMap, true);
  Teuchos::ArrayRCP<Scalar> colGIDVecData = colGIDVec->getDataNonConst(0);
  Teuchos::RCP<Vector> contColGIDVec = VectorFactory::Build(colMap, true);
  contColGIDVec->doImport(*colGIDVec, *rowColImport, Xpetra::INSERT);

  Teuchos::ArrayRCP<const Scalar> constColGIDVecData = colGIDVec->getData(0);
  std::vector<GlobalOrdinal> contColGIDs;
  for (size_t i = 0; i < contColGIDVec->getLocalLength(); ++i) {
    contColGIDs.push_back(constColGIDVecData[i]);
  }
  contColMap = MapFactory::Build(rowMap->lib(), Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
      contColGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), rowMap->getComm());

  return;
}

//! Recursive V-cycle in region fashion
void vCycle(const int l, ///< ID of current level
    const int numLevels, ///< Total number of levels
    const int maxFineIter, ///< max. sweeps on fine and intermediate levels
    const int maxCoarseIter, ///< max. sweeps on coarse level
    const double omega, ///< damping parameter for Jacobi smoother
    const int maxRegPerProc, ///< Max number of regions per process
    std::vector<Teuchos::RCP<Vector> >& fineRegX, ///< solution
    std::vector<Teuchos::RCP<Vector> > fineRegB, ///< right hand side
    Teuchos::Array<std::vector<Teuchos::RCP<Matrix> > > regMatrices, ///< Matrices in region layout
    Teuchos::Array<std::vector<Teuchos::RCP<Matrix> > > regProlong, ///< Prolongators in region layout
    Teuchos::Array<Teuchos::RCP<Map> > compRowMaps, ///< composite maps
    Teuchos::Array<std::vector<Teuchos::RCP<Map> > > quasiRegRowMaps, ///< quasiRegional row maps
    Teuchos::Array<std::vector<Teuchos::RCP<Map> > > regRowMaps, ///< regional row maps
    Teuchos::Array<std::vector<Teuchos::RCP<Import> > > regRowImporters, ///< regional row importers
    Teuchos::Array<std::vector<Teuchos::RCP<Vector> > > regInterfaceScalings, ///< regional interface scaling factors
    Teuchos::RCP<Matrix> coarseCompMat ///< Coarsest level composite operator
    )
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  const Scalar scalarZero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar scalarOne = Teuchos::ScalarTraits<Scalar>::one();

  if (l < numLevels - 1) { // fine or intermediate levels

//    std::cout << "level: " << l << std::endl;

    // pre-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
        regInterfaceScalings[l], maxRegPerProc, compRowMaps[l],
        quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    std::vector<RCP<Vector> > regRes(maxRegPerProc);
    createRegionalVector(regRes, maxRegPerProc, regRowMaps[l]);
    computeResidual(regRes, fineRegX, fineRegB, regMatrices[l], compRowMaps[l],
        quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    // Transfer to coarse level
    std::vector<RCP<Vector> > coarseRegX(maxRegPerProc);
    std::vector<RCP<Vector> > coarseRegB(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      coarseRegX[j] = VectorFactory::Build(regRowMaps[l+1][j], true);
      coarseRegB[j] = VectorFactory::Build(regRowMaps[l+1][j], true);

      regRes[j]->elementWiseMultiply(scalarOne, *regRes[j], *((regInterfaceScalings[l])[j]), scalarZero);

      regProlong[l+1][j]->apply(*regRes[j], *coarseRegB[j], Teuchos::TRANS);
      TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
      TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegB[j]->getMap()));
    }
    sumInterfaceValues(coarseRegB, compRowMaps[l+1], maxRegPerProc,
        quasiRegRowMaps[l+1], regRowMaps[l+1], regRowImporters[l+1]);

    // Call V-cycle recursively
    vCycle(l+1, numLevels, maxFineIter, maxCoarseIter, omega, maxRegPerProc,
        coarseRegX, coarseRegB, regMatrices, regProlong, compRowMaps,
        quasiRegRowMaps, regRowMaps, regRowImporters, regInterfaceScalings, coarseCompMat);

    // Transfer coarse level correction to fine level
    std::vector<RCP<Vector> > regCorrection(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      regCorrection[j] = VectorFactory::Build(regRowMaps[l][j], true);
      regProlong[l+1][j]->apply(*coarseRegX[j], *regCorrection[j]);
      TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegX[j]->getMap()));
      TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regCorrection[j]->getMap()));
    }

    // apply coarse grid correction
    for (int j = 0; j < maxRegPerProc; j++) {
      fineRegX[j]->update(scalarOne, *regCorrection[j], scalarOne);
    }

//    std::cout << "level: " << l << std::endl;

    // post-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
        regInterfaceScalings[l], maxRegPerProc, compRowMaps[l], quasiRegRowMaps[l],
        regRowMaps[l], regRowImporters[l]);
  }

  else {

    // create row and column maps with continuous GIDs
    RCP<Map> contigRowMap = Teuchos::null; //= MapFactory::Build(coarseCompMat->getRowMap());
    RCP<Map> contigColMap = Teuchos::null; //MapFactory::Build(coarseCompMat->getColMap());
    createContinuousCoarseLevelMaps(coarseCompMat->getRowMap(),
        coarseCompMat->getColMap(), contigRowMap, contigColMap);

    // Store non contiguous maps for later
    RCP<const Map> noncontigRowMap = coarseCompMat->getRowMap();
    RCP<const Map> noncontigColMap = coarseCompMat->getColMap();

    // Create composite error vector (zero initial guess)
    RCP<Vector> compX = VectorFactory::Build(coarseCompMat->getRowMap(), true);

    // Create composite right-hand side vector
    RCP<Vector> compRhs = VectorFactory::Build(coarseCompMat->getRowMap(), true);
    {
      for (int j = 0; j < maxRegPerProc; j++) {
        RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(regInterfaceScalings[l][j]->getMap());
        inverseInterfaceScaling->reciprocal(*regInterfaceScalings[l][j]);
        fineRegB[j]->elementWiseMultiply(scalarOne, *fineRegB[j], *inverseInterfaceScaling, scalarZero);
      }

      regionalToComposite(fineRegB, compRhs, maxRegPerProc, quasiRegRowMaps[l],
        regRowImporters[l], Xpetra::ADD);
    }

    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Coarse level solver not migrated to Xpetra, yet.");

//    {
//      Level level;
//      RCP<FactoryManager> factoryHandler = rcp(new FactoryManager());
//      factoryHandler->SetKokkosRefactor(false);
//      level.SetFactoryManager(factoryHandler);
//      level.SetLevelID(0);
//      level.Set("A", coarseCompMat);
//      level.setlib(compX->getMap()->UnderlyingLib());
//    }


    {
      /*

       1. Create DirectSolver by calling its constructor
       2. Call DirectSolver::Copy() to obtain Amesos/Amesos2 object wrapped into a SmootherPrototype
       3. Call Setup() and Apply() on the SmootherPrototype

       */
    }

//    {
//      using DirectSolver = MueLu::DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using FactoryManager = MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using Hierarchy = MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using SmootherPrototype = MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using SmootherFactory = MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//
//      RCP<Hierarchy> H = rcp(new Hierarchy(coarseCompMat));
//
//      RCP<SmootherPrototype> coarseProto = rcp(new DirectSolver("Klu"));
//      RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(coarseProto));
//
//      FactoryManager M;
//      M.SetFactory("CoarseSolver", coarseSolveFact);
//
//      H->Setup(M, 0, 1);
//      H->Iterate(*compRhs, *compX, 1);
//    }

    // Replace non-continuos maps by continuous maps
//    compX->replaceMap(contigRowMap);
//    compRhs->replaceMap(contigRowMap);



//    coarseCompMat->replaceRowMap(*contigRowMap);

    //    TEUCHOS_ASSERT(err==0);
//    err = coarseCompMat->ReplaceColMap(*contigColMap);
//    TEUCHOS_ASSERT(err==0);
//    err = coarseCompMat->ExpertStaticFillComplete(*contigRowMap, *contigRowMap);
//    TEUCHOS_ASSERT(err==0);
//
//    // create a linear problem object
//    Epetra_LinearProblem problem(coarseCompMat.get(), &(*compX), &(*compRhs));
//
//    // Direct solver
//    {
//      Teuchos::ParameterList pList;
//      pList.set("PrintTiming",true);
//      pList.set("PrintStatus",true);
//      pList.set("MaxProcs", coarseCompMat->Comm().NumProc());
//
//      Amesos Factory;
//      Amesos_BaseSolver* solver = Factory.Create("Amesos_Umfpack", problem);
//      TEUCHOS_ASSERT(solver!=NULL);
//
//      solver->SetParameters(pList);
//      solver->SetUseTranspose(false);
//
//      solver->SymbolicFactorization();
//      solver->NumericFactorization();
//      solver->Solve();
//    }
//
//    // Replace maps with the original non-continuous maps
//    compX->ReplaceMap(*noncontigColMap);
//    compRhs->ReplaceMap(*noncontigRowMap);
  }

  return;
}

Teuchos::Array<LocalOrdinal> setLocalNodesPerDim(const std::string& problemType,
    const bool doing1D, const Map& rowMap, const LocalOrdinal dimX, const LocalOrdinal dimY, const LocalOrdinal dimZ = 1)
{
  // Number of nodes per x/y/z-direction per processor
  Teuchos::Array<LocalOrdinal> lNodesPerDim(3);

  if (problemType == "structured") {
    if (doing1D) { // One-dimensional problems
      lNodesPerDim[0] = rowMap.getNodeNumElements();
      lNodesPerDim[1] = 1;
      lNodesPerDim[2] = 1;
    }
    else { // Two-dimensional problems

      // caseFifteen
//      {
//        lNodesPerDim[0] = 4;
//        lNodesPerDim[1] = 4;
//        lNodesPerDim[2] = 1;
//      }
//
//      // caseSixteen
//      {
//        lNodesPerDim[0] = 7;
//        lNodesPerDim[1] = 7;
//        lNodesPerDim[2] = 1;
//      }
//
//      // caseSeventeen
//      {
//        lNodesPerDim[0] = 31;
//        lNodesPerDim[1] = 31;
//        lNodesPerDim[2] = 1;
//      }
//
//      // caseEightteen / caseNineteen
//      {
//        lNodesPerDim[0] = 16;
//        lNodesPerDim[1] = 16;
//        lNodesPerDim[2] = 1;
//      }

      // caseTwenty
      {
        lNodesPerDim[0] = 10;
        lNodesPerDim[1] = 10;
        lNodesPerDim[2] = 1;
      }
    }
  }
  else {

    // Circle-in-a-disk example on 5 processors
    {
      if (rowMap.getComm()->getRank() == 4)
      {
        // This is the unstructured region, so we set dummy values
        lNodesPerDim [0] = -1;
        lNodesPerDim [1] = -1;
        lNodesPerDim [2] = -1;
      }
      else
      {
        lNodesPerDim[0] = dimX;
        lNodesPerDim[1] = dimY;
        lNodesPerDim[2] = dimZ;
      }
    }
  }

  return lNodesPerDim;
}

// Select the type of aggregation in each region
std::string setAggregationTypePerRegion (const std::string& problemType, const int myRank)
{
  std::string aggregationType;

  if (problemType == "structured")
  {
    aggregationType = "structured";
  }
  else
  {
    // Circle-in-a-disk example on 5 processors
    {
      if (myRank == 4)
        aggregationType = "uncoupled";
      else
        aggregationType = "structured";
    }
  }

  return aggregationType;
}

/* To run the region MG solver, first run the Matlab program 'createInput.m'
 * to write a bunch of files with region information to the disk.
 * Then start this executable with the appropriate number of MPI ranks.
 * See unpublished latex document
 * ``Using Trilinos Capabilities to Facilitate Composite-Regional Transitions''
 * for further information.
 */

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm CommEpetra(MPI_COMM_WORLD);

  MPI_Group world_group;
  MPI_Comm_group(CommEpetra.Comm(),&world_group);
#else
  Epetra_SerialComm CommEpetra;
#endif

  // wrap communicator into Teuchos::Comm
  Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();

  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  int  myRank;
  FILE *fp;
  char command[40];
  bool doing1D = false;
  int numDimensions = 0;
  int  globalNx, globalNy;
  std::string xmlFileName;
  std::string problemType;
  std::string regionDataDirectory;

  myRank = Comm->getRank();

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;

  Comm->barrier();

  // read xml filename from command line
  Teuchos::CommandLineProcessor clp;
  {
    // define a help message
    clp.setDocString("Driver for region multigrid\n\nUse Matlab script 'createInput.m' to create necessary input data on the hard dist.\n\nProvide filename of MueLu xml configuration via '--xml=...'.");

    // define command line arguments
    clp.setOption("xml", &xmlFileName, "filename of xml-file with MueLu configuration", true);
    clp.setOption("probType", &problemType, "Problem type [structured, hybrid]", true);
    clp.setOption("regDataDir", &regionDataDirectory, "directory with all region information/files", true);

    // force user to specify all options
    clp.recogniseAllOptions(true);

    /* We now parse the command line where argc and argv are passed to
     * the parse method.
     */
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = clp.parse(argc, argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      return 0;
    }
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return 1; // Error!
    }

    // Check for valid command line arguments
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(problemType == "structured" || problemType == "hybrid"),
        "Unknown problem type. Use either 'structured' or 'hybrid'.\n");
  }

  Comm->barrier();

  // read maxRegPerGID (maximum # of regions shared by any one node) and
  //      maxRegPerProc (maximum # of partial regions owned by any proc)
  //      whichCase (MultipleRegionsPerProc: No regions spans across multiple
  //                                procs but a proc might own multiple regions
  //                 RegionsSpanProcs: processors own only a piece of 1 region
  //                                but regions may span across many procs

  // Provide some feedback to the user
  if (myRank == 0) {

    std::cout << "User input:" << std::endl
        << "  xml-file with MueLu configuration: " << xmlFileName << std::endl
        << "  problem type: " << problemType << std::endl
        << "  Path to directory with region data: " << regionDataDirectory << std::endl;
  }
  Comm->barrier();

  // Define some basic scalar values
  Scalar scalarZero = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar scalarOne = Teuchos::ScalarTraits<Scalar>::one();

  int maxRegPerGID = 0;
  int maxRegPerProc = 0;
  int whichCase = 0;
  int iii = 0;

  // Read region info from file
  {
    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegionInfo_" << myRank;
    while ((fp = fopen(fileNameSS.str().c_str(),"r") ) == NULL) sleep(1);

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
  }

  Comm->barrier();

  // check for 1D or 2D problem
  if (globalNy == 1)
    doing1D = true;
  else
    doing1D = false;

  if (doing1D)
    numDimensions = 1;
  else
    numDimensions = 2;

  // ******************************************************************
  // Application Specific Data for LIDregion()
  // ******************************************************************
  struct widget appData;                            // ****************
  std::vector<GlobalOrdinal>  minGIDComp(maxRegPerProc);      // ****************
  std::vector<GlobalOrdinal>  maxGIDComp(maxRegPerProc);      // ****************
  // ******************************************************************
  // ******************************************************************


  std::vector<LocalOrdinal> myRegions; // regions that myRank owns
  Teuchos::RCP<CrsMatrixWrap> AComp = Teuchos::null; // composite form of matrix
  Teuchos::RCP<Matrix> ACompSplit = Teuchos::null; // composite form of matrix
  Teuchos::RCP<Map> mapComp = Teuchos::null; // composite map used to build AComp

  // regionsPerGID[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix row Map.
  Teuchos::RCP<MultiVector> regionsPerGID = Teuchos::null;

  // regionsPerGIDWithGhosts[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix col Map.
  Teuchos::RCP<MultiVector> regionsPerGIDWithGhosts = Teuchos::null;

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

  std::vector <Teuchos::RCP<Map> > rowMapPerGrp(maxRegPerProc); // row map associated with myRank's ith region in composite layout
  std::vector <Teuchos::RCP<Map> > colMapPerGrp(maxRegPerProc); // column map associated with myRank's ith region in composite layout
  std::vector <Teuchos::RCP<Map> > revisedRowMapPerGrp(maxRegPerProc); // revised row map associated with myRank's ith region for regional layout
  std::vector <Teuchos::RCP<Map> > revisedColMapPerGrp(maxRegPerProc); // revised column map associated with myRank's ith region for regional layout

  std::vector<Teuchos::RCP<Import> > rowImportPerGrp(maxRegPerProc); // row importers per group
  std::vector<Teuchos::RCP<Import> > colImportPerGrp(maxRegPerProc); // column importers per group
  std::vector<Teuchos::RCP<Export> > rowExportPerGrp(maxRegPerProc); // row exporters per group

  std::vector<Teuchos::RCP<Matrix> > quasiRegionGrpMats(maxRegPerProc); // region-wise matrices with quasiRegion maps (= composite GIDs)
  std::vector<Teuchos::RCP<Matrix> > regionGrpMats(maxRegPerProc); // region-wise matrices in true region layout with unique GIDs for replicated interface DOFs

  Teuchos::RCP<Epetra_Map> coarseCompRowMap; // composite row map on the coarse grid
  std::vector<Teuchos::RCP<Map> > coarseRowMapPerGrp(maxRegPerProc); // region-wise row map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Map> > coarseQuasiRowMapPerGrp(maxRegPerProc); // region-wise row map in quasiRegion layout with original GIDs from fine level
  std::vector<Teuchos::RCP<Map> > coarseColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Map> > coarseAltColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Import> > coarseRowImportPerGrp(maxRegPerProc); // coarse level row importer per group
  std::vector<Teuchos::RCP<Matrix> > regionGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Matrix> > regionAltGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Matrix> > regCoarseMatPerGrp(maxRegPerProc); // coarse level operator 'RAP' in region layout

  std::vector<Teuchos::RCP<Vector> > regionInterfaceScaling(maxRegPerProc);
  std::vector<Teuchos::RCP<Vector> > coarseRegionInterfaceScaling(maxRegPerProc);

  std::vector<Teuchos::RCP<Vector> > regNspViolation(maxRegPerProc); // violation of nullspace property in region layout

  Teuchos::RCP<Vector> compX = Teuchos::null; // initial guess for truly composite calculations
  Teuchos::RCP<Vector> compY = Teuchos::null; // result vector for truly composite calculations
  Teuchos::RCP<Vector> regYComp = Teuchos::null; // result vector in composite layout, but computed via regional operations
  Teuchos::RCP<Vector> nspViolation = Teuchos::null; // violation of nullspace property in composite layout

  std::vector<Teuchos::RCP<Vector> > quasiRegX(maxRegPerProc); // initial guess associated with myRank's ith region in quasiRegional layout
  std::vector<Teuchos::RCP<Vector> > quasiRegY(maxRegPerProc); // result vector associated with myRank's ith region in quasiRegional layout
  std::vector<Teuchos::RCP<Vector> > regX(maxRegPerProc); // initial guess associated with myRank's ith region in regional layout
  std::vector<Teuchos::RCP<Vector> > regY(maxRegPerProc); // result vector associated with myRank's ith region in regional layout

  std::vector<LocalOrdinal> intIDs; // LIDs of interface DOFs
  std::vector<std::vector<LocalOrdinal> > regIntIDs(maxRegPerProc); // LIDs of interface DOFs

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
  Array<RCP<Map> > compRowMaps; // composite row maps on each level
  Array<RCP<Map> > compColMaps; // composite columns maps on each level
  Array<std::vector<RCP<Map> > > regRowMaps; // regional row maps on each level
  Array<std::vector<RCP<Map> > > quasiRegRowMaps; // quasiRegional row maps on each level
  Array<std::vector<RCP<Map> > > regColMaps; // regional column maps on each level
  Array<std::vector<RCP<Map> > > quasiRegColMaps; // quasiRegional column maps on each level
  Array<std::vector<RCP<Matrix> > > regMatrices; // regional matrices on each level
  Array<std::vector<RCP<Matrix> > > regProlong; // regional prolongators on each level
  Array<std::vector<RCP<Import> > > regRowImporters; // regional row importers on each level
  Array<std::vector<RCP<Vector> > > regInterfaceScalings; // regional interface scaling factors on each level

  Array<std::vector<std::vector<LocalOrdinal> > > interfaceLIDs; // local IDs of interface nodes on each level in each group
  Array<std::vector<std::vector<std::vector<GlobalOrdinal> > > > interfaceGIDPairs; // pairs of GIDs of interface nodes on each level in each group

  Teuchos::RCP<Matrix> coarseCompOp = Teuchos::null;

  /* The actual computations start here. It's a sequence of operations to
   * - read region data from files
   * - create composite and region matrices
   * - create a region multigrid hierarchy
   * - run a region-type V-cycle
   */

  Comm->barrier();

  // Load composite map
  {
    std::cout << myRank << " | Loading composite map ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myCompositeMap_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    if (fp == NULL)
      std::cout << std::endl << ">>> Check number of MPI ranks!" << std::endl << std::endl;
    TEUCHOS_ASSERT(fp!=NULL);

    Teuchos::Array<GlobalOrdinal> fileData; // composite GIDs
    GlobalOrdinal i;
    while (fscanf(fp, "%d", &i) != EOF)
      fileData.push_back(i);

    Teuchos::ArrayView<GlobalOrdinal> fileDataView = fileData;

    mapComp = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), fileDataView,
        Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
//    mapComp->describe(*fos, Teuchos::VERB_EXTREME);
  }

  Comm->barrier();

  // Read matrix from file
  {
    std::cout << myRank << " | Reading matrix from file ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/Amat.mm";

    AComp = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(IO::Read(fileNameSS.str(), lib, Comm, false));
//    AComp->describe(*fos, Teuchos::VERB_EXTREME);
  }

  Comm->barrier();

  // Load and communicate region assignments
  {
    std::cout << myRank << " | Loading and communicating region assignments ..." << std::endl;

    regionsPerGID = MultiVectorFactory::Build(AComp->getRowMap(), maxRegPerGID, true);

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegionAssignment_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    int k;
    RCP<Vector> jthRegions; // pointer to jth column in regionsPerGID
    for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(mapComp->getNodeNumElements()); i++) {
      for (int j = 0; j < maxRegPerGID; j++) {
        jthRegions = regionsPerGID->getVectorNonConst(j);

        if (fscanf(fp, "%d", &k) == EOF) {
          fprintf(stderr, "not enough region assignments\n");
          exit(1);
        }
        jthRegions->replaceLocalValue(i, (Scalar) k);

        // identify interface DOFs. An interface DOF is assigned to at least two regions.
        if (j > 0 and k != -1)
          intIDs.push_back(i);
      }
    }

    // make extended Region Assignments
    RCP<Import> Importer = ImportFactory::Build(mapComp, AComp->getColMap());
    regionsPerGIDWithGhosts = MultiVectorFactory::Build(AComp->getColMap(), maxRegPerGID, true);
    regionsPerGIDWithGhosts->doImport(*regionsPerGID, *Importer, Xpetra::INSERT);
//    regionsPerGIDWithGhosts->describe(*fos, Teuchos::VERB_EXTREME);
  }

  Comm->barrier();

  // Load Regions
  {
    std::cout << myRank << " | Loading regions ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegions_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    while (fgets(command, 80, fp) != NULL) {
      int i;
      sscanf(command, "%d", &i);
      myRegions.push_back(i);
    }
  }

  Comm->barrier();

  std::vector<LocalOrdinal> genericVector;
  // Load AppData for LID region
  {
    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myAppData_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    int retval; // dummy return value for fscanf() to suppress compiler warnings

    if (problemType == "structured") {

      // Fill this with dummy entries. Will not be used in fully structured problems.
      genericVector.resize(3);
      genericVector[inpData_isStructured] = 0;
      genericVector[inpData_ownedX] = -1;
      genericVector[inpData_ownedY] = -1;

      if (doing1D) {
        int minGID, maxGID;
        for (int i = 0; i < (int) myRegions.size(); i++) {
          retval = fscanf(fp,"%d%d",&minGID,&maxGID);
          if(retval == 0) {std::cout << "Something probably went wrong while reading minGID and maxGID from file!" << std::endl;}
          minGIDComp[i] = minGID;
          maxGIDComp[i] = maxGID;
        }
        appData.gDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lowInd= (int *) malloc(sizeof(int)*3*myRegions.size());
        for (int i = 0; i < (int) myRegions.size(); i++) {
          retval = fscanf(fp,"%d%d%d",&(appData.gDim[3*i]),&(appData.gDim[3*i+1]),&(appData.gDim[3*i+2]));
          retval = fscanf(fp,"%d%d%d",&(appData.lDim[3*i]),&(appData.lDim[3*i+1]),&(appData.lDim[3*i+2]));
          retval = fscanf(fp,"%d%d%d",&(appData.lowInd[3*i]),&(appData.lowInd[3*i+1]),&(appData.lowInd[3*i+2]));
        }
        appData.minGIDComp = minGIDComp.data();
        appData.maxGIDComp = maxGIDComp.data();
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
          retval = fscanf(fp,"%d%d%d",&(appData.gDim[3*i]),&(appData.gDim[3*i+1]),&(appData.gDim[3*i+2]));
          retval = fscanf(fp,"%d%d%d",&(appData.lDimx[i]),&(appData.lDimy[i]),&garbage);
          retval = fscanf(fp,"%d%d%d",&(appData.relcornerx[i]),&(appData.relcornery[i]),&garbage);
          retval = fscanf(fp,"%d%d%d",&(appData.trueCornerx[i]),&(appData.trueCornery[i]),&garbage);
        }
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Hybrid case not migrated to Xpetra, yet.");
//      int isStructured;
//      retval = fscanf(fp,"%d",&isStructured);
//      if (isStructured == 0) {
//         genericVector.resize(3);
//         genericVector[inpData_isStructured] = 0;
//         retval = fscanf(fp,"%d",&(genericVector[inpData_ownedX]));
//         genericVector[inpData_ownedY] = 1;
//      }
//      else {
//         int Cx, Cy, Rx, Ry;
//         int ownedX,ownedY;
//         retval = fscanf(fp,"%d%d",&Rx,&Ry);
//         retval = fscanf(fp,"%d%d",&ownedX,&ownedY);
//         retval = fscanf(fp,"%d%d",&Cx,&Cy);
//         int *rowptr, *cols;  double *values;
//         AComp->ExtractCrsDataPointers(rowptr, cols, values);
//         fillCircleSquareData(rowptr, cols, ownedX,ownedY, Cx, Cy, Rx, Ry,genericVector);
//      }
//      fclose(fp);
    }

    appData.maxRegPerGID = maxRegPerGID;
    appData.myRegions = myRegions.data();
    appData.colMap = (Map*) &*(AComp->getColMap());
    appData.regionsPerGIDWithGhosts = regionsPerGIDWithGhosts;
    appData.myRank = myRank;
  }

  Comm->barrier();

  // Make group region row maps
  {
    std::cout << myRank << " | Creating region group row maps ..." << std::endl;

    Teuchos::Array<GlobalOrdinal> rowGIDsReg;
    const Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = AComp->getColMap()->getNodeElementList();

    for (int k = 0; k < (int) myRegions.size(); k++) {
      rowGIDsReg.resize(0);
      std::vector<int> tempRegIDs(AComp->getColMap()->getNodeNumElements());
      for (int i = 0; i < (int) AComp->getColMap()->getNodeNumElements(); i++) {

        if (problemType == "structured") {
          if (doing1D)
            tempRegIDs[i] = LIDregion(&appData, i, k);
          else
            tempRegIDs[i] = LID2Dregion(&appData, i, k);
        }
        else {
          tempRegIDs[i] = LIDregionCircleSquare(genericVector.data(), i, k);
//          printf("%d: LIDRegion(Composite LID=%d or GID=%d) = %d\n",myRank,i,colGIDsComp[i],tempRegIDs[i]); fflush(stdout);
        }
      }

      std::vector<int> idx(tempRegIDs.size());
      std::iota(idx.begin(), idx.end(), 0);
      sort(idx.begin(), idx.end(), [tempRegIDs](int i1,int i2) { return tempRegIDs[i1] < tempRegIDs[i2];});

      for (int i = 0; i < Teuchos::as<int>(AComp->getColMap()->getNodeNumElements()); i++) {
        if (tempRegIDs[idx[i]] != -1)
          rowGIDsReg.push_back(colGIDsComp[idx[i]]);
      }
      rowMapPerGrp[k] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), rowGIDsReg,
          Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
    }

    for (int k=(int) myRegions.size(); k < maxRegPerProc; k++) {
      rowMapPerGrp[k] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
          Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
    }
  }

  Comm->barrier();

  // Make group region column maps
  {
    std::cout << myRank << " | Creating region group column maps ..." << std::endl;

    if (whichCase == MultipleRegionsPerProc) {
      // clone rowMap
      for (int j=0; j < maxRegPerProc; j++) {
        if (j < (int) myRegions.size()) {
          colMapPerGrp[j] = MapFactory::Build(lib, Teuchos::OrdinalTraits<int>::invalid(),
              rowMapPerGrp[j]->getNodeElementList(), Teuchos::OrdinalTraits<int>::zero(), Comm);
        }
        else
        {
          colMapPerGrp[j] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
              Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
        }
      }
    }
    else if (whichCase == RegionsSpanProcs) {//so maxRegPerProc = 1
      Teuchos::Array<GlobalOrdinal> colIDsReg;

      // copy the rowmap
      Teuchos::ArrayView<const GlobalOrdinal> rowGIDsReg = rowMapPerGrp[0]->getNodeElementList();
      for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(rowMapPerGrp[0]->getNodeNumElements()); i++)
        colIDsReg.push_back(rowGIDsReg[i]);

      // append additional ghosts who are in my region and
      // for whom I have a LID
      LocalOrdinal LID;
      Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp =  AComp->getColMap()->getNodeElementList();
      for (std::size_t i = 0; i < AComp->getColMap()->getNodeNumElements(); i++) {
      if (problemType == "structured") {
        if (doing1D)
          LID = LIDregion(&appData, i, 0);
        else
          LID = LID2Dregion(&appData, i, 0);
      }
      else {
        LID = LIDregionCircleSquare(genericVector.data(), i, 0);
      }
        if (LID == -1) {
          for (int j = 0; j < maxRegPerGID; j++) {
            Teuchos::ArrayRCP<const double> jthRegions = regionsPerGIDWithGhosts->getData(j);
            if  ( ((int) jthRegions[i]) == myRegions[0]) {
              colIDsReg.push_back(colGIDsComp[i]);
              break;
            }
          }
        }
      }
      if ((int) myRegions.size() > 0) {
       colMapPerGrp[0] = MapFactory::Build(lib, Teuchos::OrdinalTraits<int>::invalid(), colIDsReg,
           Teuchos::OrdinalTraits<int>::zero(), Comm);
      }
      else
      {
        colMapPerGrp[0] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
          Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
      }
    }
    else { fprintf(stderr,"whichCase not set properly\n"); exit(1); }
  }

  Comm->barrier();

  // Make extended group region maps
  {
    std::cout << myRank << " | Creating extended group region maps ..." << std::endl;

    LocalOrdinal nLocal = AComp->getRowMap()->getNodeNumElements();
    LocalOrdinal nExtended = AComp->getColMap()->getNodeNumElements();
    LocalOrdinal nTotal = 0;
    Teuchos::reduceAll<LocalOrdinal, LocalOrdinal>(*Comm, Teuchos::REDUCE_SUM, nLocal, Teuchos::outArg(nTotal));

    // first step of NewGID calculation just counts the number of NewGIDs
    // and sets firstNewGID[k] such that it is equal to the number of
    // NewGIDs that have already been counted in (0:k-1) that myRank
    // is responsible for.

    RCP<Vector> firstNewGID = VectorFactory::Build(mapComp, true);
    for (int k = 0; k < nLocal-1; k++) {
      firstNewGID->replaceLocalValue(k+1, (firstNewGID->getData(0))[k]-1);
      for (int j = 0; j < maxRegPerGID; j++) {
        Teuchos::ArrayRCP<const Scalar> jthRegions = regionsPerGIDWithGhosts->getData(j);
        if (jthRegions[k] != -scalarOne) firstNewGID->sumIntoLocalValue(k+1, scalarOne);
      }
    }
    // So firstNewGID[nLocal-1] is number of NewGIDs up to nLocal-2
    // To account for newGIDs associated with nLocal-1, we'll just
    // use an upper bound (to avoid the little loop above).
    // By adding maxRegPerGID-1 we account for the maximum
    // number of possible newGIDs due to last composite id
    Teuchos::ArrayRCP<const Scalar> firstNewGIDData = firstNewGID->getData(0);
    GlobalOrdinal upperBndNumNewGIDs = Teuchos::as<GlobalOrdinal>(firstNewGIDData[nLocal-1]) + maxRegPerGID-1;
    GlobalOrdinal upperBndNumNewGIDsAllProcs;
    Teuchos::reduceAll<GlobalOrdinal,GlobalOrdinal>(*Comm, Teuchos::REDUCE_MAX, upperBndNumNewGIDs, Teuchos::outArg(upperBndNumNewGIDsAllProcs));

//    std::cout << "upperBndNumNewGIDsAllProcs: " << upperBndNumNewGIDsAllProcs << std::endl;

    // Now that we have an upper bound on the maximum number of
    // NewGIDs over all procs, we sweep through firstNewGID again
    // to assign ids to the first NewGID associated with each row of
    // regionsPerGIDWithGhosts (by adding an offset)
    ArrayRCP<Scalar> firstNewGIDDataNonConst = firstNewGID->getDataNonConst(0);
    for (LocalOrdinal k = 0; k < nLocal; k++) {
//      const GlobalOrdinal tmpGID = firstNewGIDDataNonConst[k];
//      firstNewGIDDataNonConst[k] = tmpGID + upperBndNumNewGIDsAllProcs*myRank+nTotal;
      firstNewGID->sumIntoLocalValue(k, upperBndNumNewGIDsAllProcs*myRank+nTotal);
    }

    RCP<Import> Importer = ImportFactory::Build(mapComp, AComp->getColMap());
    RCP<Vector> firstNewGIDWithGhost = VectorFactory::Build(Importer->getTargetMap()); //AComp->getColMap());
    firstNewGIDWithGhost->doImport(*firstNewGID, *Importer, Xpetra::INSERT);

    Teuchos::Array<GlobalOrdinal> revisedGIDs;

    for (std::size_t k = 0; k < myRegions.size(); k++) {
      revisedGIDs.resize(0);
      int curRegion = myRegions[k];
      Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = AComp->getColMap()->getNodeElementList();
      Teuchos::Array<LocalOrdinal> tempRegIDs(nExtended);

      // must put revisedGIDs in application-provided order by
      // invoking LIDregion() and sorting
      for (LocalOrdinal i = 0; i < nExtended; i++) {
        if (problemType == "structured") {
          if (doing1D)
            tempRegIDs[i] = LIDregion(&appData, i, k);
          else
            tempRegIDs[i] = LID2Dregion(&appData, i, k);
        }
        else {
          tempRegIDs[i] = LIDregionCircleSquare(genericVector.data(), i, k);
        }
      }
      Teuchos::Array<LocalOrdinal> idx(tempRegIDs.size());
      std::iota(idx.begin(),idx.end(),0);
      sort(idx.begin(),idx.end(),[tempRegIDs](int i1,int i2){return tempRegIDs[i1] < tempRegIDs[i2];});

      // Now sweep through regionsPerGIDWithGhosts looking for those
      // associated with curRegion and record the revisedGID
      int j;
      for (LocalOrdinal i = 0; i < nExtended; i++) {

        if (tempRegIDs[idx[i]] != -1) {// if a valid LID for this region
          for (j = 0; j < maxRegPerGID; j++) {
            Teuchos::ArrayRCP<const Scalar> jthRegions = regionsPerGIDWithGhosts->getData(j);
            if (((int) jthRegions[idx[i]]) == curRegion) break;
          }

          // (*regionsPerGIDWithGhosts)[0] entries keep original GID
          // while others use firstNewGID to determine NewGID

          if (j == 0) revisedGIDs.push_back(colGIDsComp[idx[i]]);
          else if (j < maxRegPerGID) {
             revisedGIDs.push_back((GlobalOrdinal) firstNewGIDWithGhost->getData(0)[idx[i]] + j - 1);
          }

          // add entry to listDulicatedGIDs
        }
      }

      revisedRowMapPerGrp[k] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), revisedGIDs,
          Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),Comm);

      // now append more stuff to handle ghosts ... needed for
      // revised version of column map

      for (LocalOrdinal i = 0; i < nExtended; i++) {
        if (tempRegIDs[i] == -1) {// only in revised col map
                   // note: since sorting not used when making the
                   // original regional colmap, we can't use
                   // it here either ... so no idx[]'s.
          for (j = 0; j < maxRegPerGID; j++) {
            Teuchos::ArrayRCP<const Scalar> jthRegions = regionsPerGIDWithGhosts->getData(j);
            if  ( ((int) jthRegions[i]) == curRegion) break;
          }
          // (*regionsPerGIDWithGhosts)[0] entries keep original GID
          // while others use firstNewGID to determine NewGID

          if (j == 0) revisedGIDs.push_back(colGIDsComp[i]);
          else if (j < maxRegPerGID) {
            revisedGIDs.push_back((int) firstNewGIDWithGhost->getData(0)[i] + j - 1);
          }
        }
      }
      revisedColMapPerGrp[k] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
          revisedGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
    }
    for (std::size_t k = myRegions.size(); k < Teuchos::as<std::size_t>(maxRegPerProc); k++) {
      revisedRowMapPerGrp[k] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
          Teuchos::OrdinalTraits<GlobalOrdinal>::zero() ,Comm);
      revisedColMapPerGrp[k] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
          Teuchos::OrdinalTraits<GlobalOrdinal>::zero() ,Comm);
    }

    // Setup importers
    for (int j = 0; j < maxRegPerProc; j++) {
      rowImportPerGrp[j] = ImportFactory::Build(mapComp, rowMapPerGrp[j]);
      colImportPerGrp[j] = ImportFactory::Build(mapComp, colMapPerGrp[j]);
    }
  }

  Comm->barrier();

  // Make quasiRegion matrices
  {
    std::cout << myRank << " | Forming quasiRegion matrices ..." << std::endl;

    /* We use the edge-based splitting, i.e. we first modify off-diagonal
     * entries in the composite matrix, then decompose it into region matrices
     * and finally take care of diagonal entries by enforcing the nullspace
     * preservation constraint.
     */

    // copy and modify the composite matrix
    ACompSplit = MatrixFactory::BuildCopy(AComp);

    for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(ACompSplit->getNodeNumRows()); row++) { // loop over local rows of composite matrix
      GlobalOrdinal rowGID = ACompSplit->getRowMap()->getGlobalElement(row);
      std::size_t numEntries = ACompSplit->getNumEntriesInLocalRow(row); // number of entries in this row
      Teuchos::Array<Scalar> vals(numEntries); // non-zeros in this row
      Teuchos::Array<LocalOrdinal> inds(numEntries); // local column indices
      ACompSplit->getLocalRowCopy(row, inds, vals, numEntries);

      for (std::size_t c = 0; c < Teuchos::as<std::size_t>(inds.size()); ++c) { // loop over all entries in this row
        LocalOrdinal col = inds[c];
        GlobalOrdinal colGID = ACompSplit->getColMap()->getGlobalElement(col);
        std::vector<int> commonRegions;
        if (rowGID != colGID) { // Skip the diagonal entry. It will be processed later.
          commonRegions = findCommonRegions(rowGID, colGID, *regionsPerGIDWithGhosts);
        }

        if (commonRegions.size() > 1) {
          vals[c] *= 0.5;
        }
      }

      ACompSplit->resumeFill();
      ACompSplit->replaceLocalValues(row, inds, vals);
    }

    // Import data from ACompSplit into the quasiRegion matrices
    for (int j = 0; j < maxRegPerProc; j++) {
      quasiRegionGrpMats[j] = MatrixFactory::Build(rowMapPerGrp[j], colMapPerGrp[j],
          Teuchos::OrdinalTraits<int>::invalid(), Xpetra::DynamicProfile);
      quasiRegionGrpMats[j]->doImport(*ACompSplit, *(rowImportPerGrp[j]), Xpetra::INSERT);
      quasiRegionGrpMats[j]->fillComplete();
    }
  }

  Comm->barrier();

  // Make region matrices
  {
    std::cout << myRank << " | Forming region matrices ..." << std::endl;

    // Copy data from quasiRegionGrpMats, but into new map layout
    {
      for (int j = 0; j < maxRegPerProc; j++) {
        regionGrpMats[j] = rcp(new CrsMatrixWrap(revisedRowMapPerGrp[j], revisedColMapPerGrp[j], 9, Xpetra::DynamicProfile));

        // Extract current region CrsMatrix
        RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionGrpMats[j])->getCrsMatrix();

        // Extract current quasi-region CrsMatrix
        RCP<CrsMatrix> quasiRegionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegionGrpMats[j])->getCrsMatrix();

        // Pull out the data from the quasi-region CrsMatrix
        Teuchos::ArrayRCP<const size_t> rowptrQuasiRegion;
        Teuchos::ArrayRCP<const LocalOrdinal> colindQuasiRegion;
        Teuchos::ArrayRCP<const Scalar> valuesQuasiRegion;
        quasiRegionCrsMat->getAllValues(rowptrQuasiRegion, colindQuasiRegion, valuesQuasiRegion);

        // Do a deep copy of values
        // (at least we've been doing deep copies so far, maybe we could do shallow copies to save time?)
        Teuchos::ArrayRCP<size_t> rowptrRegion(rowptrQuasiRegion.size());
        Teuchos::ArrayRCP<LocalOrdinal> colindRegion(colindQuasiRegion.size());
        Teuchos::ArrayRCP<Scalar> valuesRegion(valuesQuasiRegion.size());

        regionCrsMat->allocateAllValues(valuesRegion.size(), rowptrRegion, colindRegion, valuesRegion);

        for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(rowptrRegion.size()); ++idx) {
          rowptrRegion[idx] = rowptrQuasiRegion[idx];
        }

        for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(colindRegion.size()); ++idx) {
          colindRegion[idx] = colindQuasiRegion[idx];
          valuesRegion[idx] = valuesQuasiRegion[idx];
        }

        // Set and fillComplete the region CrsMatrix
        regionCrsMat->setAllValues(rowptrRegion, colindRegion, valuesRegion);
        regionCrsMat->expertStaticFillComplete(revisedRowMapPerGrp[j], revisedRowMapPerGrp[j]);
      }
    }

    // enforce nullspace constraint
    {
      // compute violation of nullspace property close to DBCs
      RCP<Vector> nspVec = VectorFactory::Build(mapComp);
      nspVec->putScalar(scalarOne);
      nspViolation = VectorFactory::Build(mapComp, true);
      AComp->apply(*nspVec, *nspViolation);

      // move to regional layout
      std::vector<RCP<Vector> > quasiRegNspViolation(maxRegPerProc);
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
        std::vector<RCP<Vector> > interfaceScaling(maxRegPerProc);
        for (int j = 0; j < maxRegPerProc; j++) {
          interfaceScaling[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
          interfaceScaling[j]->putScalar(scalarOne);
        }

        // transform to composite layout while adding interface values via the Export() combine mode
        RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(mapComp, true);
        regionalToComposite(interfaceScaling, compInterfaceScalingSum,
            maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Xpetra::ADD);

        /* transform composite layout back to regional layout. Now, GIDs associated
         * with region interface should carry a scaling factor (!= 1).
         */
        std::vector<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc);
        compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
            interfaceScaling, maxRegPerProc, rowMapPerGrp,
            revisedRowMapPerGrp, rowImportPerGrp);

        // modify its interface entries
        for (int j = 0; j < maxRegPerProc; j++) {
          RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(interfaceScaling[j]->getMap(), true);
          inverseInterfaceScaling->reciprocal(*interfaceScaling[j]);
          regNspViolation[j]->elementWiseMultiply(scalarOne, *regNspViolation[j], *inverseInterfaceScaling, scalarZero);
        }
      }
    }

    std::vector<RCP<Vector> > regNsp(maxRegPerProc);
    std::vector<RCP<Vector> > regCorrection(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      regNsp[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
      regNsp[j]->putScalar(scalarOne);

      regCorrection[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
      regionGrpMats[j]->apply(*regNsp[j], *regCorrection[j]);
    }

    RCP<Vector> regDiag = Teuchos::null;
    for (int j = 0; j < maxRegPerProc; j++) {
      regDiag = VectorFactory::Build(revisedRowMapPerGrp[j], true);
      regionGrpMats[j]->getLocalDiagCopy(*regDiag);
      regDiag->update(-scalarOne, *regCorrection[j], scalarOne, *regNspViolation[j], scalarOne);

      // Extract current region matrix in as CrsMatrix
      RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionGrpMats[j])->getCrsMatrix();
      regionCrsMat->replaceDiag(*regDiag);
    }
  }

  Comm->barrier();

//  /* Form composite operator on fine level for debug purposes */
//  {
//    std::cout << myRank << " | Forming composite operator on fine level for debug purposes ..." << std::endl;
//
////    sleep(1);
////    std::cout << myRank << " | Printing AComp ..." << std::endl;
////    Comm->barrier();
////    AComp->describe(*fos, Teuchos::VERB_EXTREME);
//
////    RCP<Matrix> recombinedAComp = MatrixFactory::Build(AComp->getCrsGraph());
//    RCP<Matrix> recombinedAComp = MatrixFactory::Build(mapComp, 3, Xpetra::DynamicProfile);
//    recombinedAComp->setAllToScalar(scalarZero);
//    regionalToComposite(regionGrpMats, recombinedAComp, maxRegPerProc, rowMapPerGrp, colMapPerGrp, rowImportPerGrp, Xpetra::ADD);
//
////    sleep(1);
////    std::cout << myRank << " | Printing recombinedAComp ..." << std::endl;
////    Comm->barrier();
////    recombinedAComp->describe(*fos, Teuchos::VERB_EXTREME);
//
//        // Warning: I'm not sure whether the diffOp is computed correctly.
////    RCP<Matrix> diffOp = Teuchos::null; //MatrixFactory::Build(AComp->getCrsGraph());
////    MatrixMatrix::TwoMatrixAdd(*AComp, false, scalarOne, *recombinedAComp, false, -scalarOne, diffOp, *fos);
////
////    sleep(1);
////    std::cout << myRank << " | Printing diffOp ..." << std::endl;
////    Comm->barrier();
////    diffOp->describe(*fos, Teuchos::VERB_EXTREME);
////    Teuchos::ScalarTraits<Scalar>::magnitudeType frobeniusNorm = diffOp->getFrobeniusNorm();
////    Comm->barrier();
////    std::cout << myRank << " | Forbenius norm = " << frobeniusNorm << std::endl;
//
//  }

  Comm->barrier();

  // Create multigrid hierarchy
  {
    typedef MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> Hierarchy;
    typedef MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> Utilities;

    std::cout << myRank << " | Setting up MueLu hierarchies ..." << std::endl;

    // A hierarchy for each group
    std::vector<RCP<Hierarchy> > regGrpHierarchy(maxRegPerProc);

    for (int j = 0; j < maxRegPerProc; j++) {

      /* Set number of nodes per processor per dimension
       *
       * We don't use the number of owned nodes provided on input.
       * Use the region dimensions instead. This is the right thing to do
       * since duplication of interface nodes has added duplicated nodes to those regions
       * where inpData_ownedX/inpData_ownedY and inpData_regionX/inpData_regionY have been different on input.
       */
      Array<int> lNodesPerDim = setLocalNodesPerDim(problemType, doing1D,
          *revisedRowMapPerGrp[j], genericVector[inpData_regionX], genericVector[inpData_regionY]);

      // Set aggregation type for each region
      std::string aggregationRegionType = setAggregationTypePerRegion(problemType, myRank);

      // create nullspace vector
      RCP<MultiVector> nullspace = MultiVectorFactory::Build(revisedRowMapPerGrp[j], 1);
      nullspace->putScalar(scalarOne);

      // create dummy coordinates vector
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType real_type;
      typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;
      RCP<RealValuedMultiVector> coordinates =
          Xpetra::MultiVectorFactory<real_type, LocalOrdinal, GlobalOrdinal, Node>::Build(revisedRowMapPerGrp[j], 3);
      coordinates->putScalar(scalarOne);

      // Read MueLu parameter list form xml file
      RCP<ParameterList> mueluParams = Teuchos::getParametersFromXmlFile(xmlFileName);

      // Insert region-specific data into parameter list
      const std::string userName = "user data";
      Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);
      userParamList.set<int>("int numDimensions", numDimensions);
      userParamList.set<Array<int> >("Array<LO> lNodesPerDim", lNodesPerDim);
      userParamList.set<std::string>("string aggregationRegionType", aggregationRegionType);
      userParamList.set<RCP<Epetra_MultiVector> >("Coordinates", coordinates);
      userParamList.set<RCP<Epetra_MultiVector> >("Nullspace", nullspace);

      // Setup hierarchy
      regGrpHierarchy[j] = MueLu::CreateXpetraPreconditioner(regionGrpMats[j], *mueluParams, coordinates, nullspace);
    }

    // resize Arrays and vectors
    {
      // resize level containers
      numLevels = regGrpHierarchy[0]->GetNumLevels();
      compRowMaps.resize(numLevels);
      compColMaps.resize(numLevels);
      regRowMaps.resize(numLevels);
      regColMaps.resize(numLevels);
      quasiRegRowMaps.resize(numLevels);
      quasiRegColMaps.resize(numLevels);
      regMatrices.resize(numLevels);
      regProlong.resize(numLevels);
      regRowImporters.resize(numLevels);
      regInterfaceScalings.resize(numLevels);

      // resize group containers on each level
      for (int l = 0; l < numLevels; ++l) {
        regRowMaps[l].resize(maxRegPerProc);
        regColMaps[l].resize(maxRegPerProc);
        quasiRegRowMaps[l].resize(maxRegPerProc);
        quasiRegColMaps[l].resize(maxRegPerProc);
        regMatrices[l].resize(maxRegPerProc);
        regProlong[l].resize(maxRegPerProc);
        regRowImporters[l].resize(maxRegPerProc);
        regInterfaceScalings[l].resize(maxRegPerProc);
      }
    }

    // Fill fine level with our data
    {
      compRowMaps[0] = mapComp;
      quasiRegRowMaps[0] = rowMapPerGrp;
      quasiRegColMaps[0] = colMapPerGrp;
      regRowMaps[0] = revisedRowMapPerGrp;
      regColMaps[0] = revisedColMapPerGrp;
      regRowImporters[0] = rowImportPerGrp;
      regMatrices[0] = regionGrpMats;

      /* MueLu stores prolongator on coarse level, so there is no prolongator
       * on the fine level. To have level containers of the same size, let's
       * just put in dummy data
       */
      std::vector<RCP<Matrix> > fineLevelProlong(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; ++j)
        fineLevelProlong[j] = Teuchos::null;
      regProlong[0] = fineLevelProlong;
    }

    /* Get coarse level matrices and prolongators from MueLu hierarchy
     * Note: fine level has been dealt with previously, so we start at level 1 here.
     */
    for (int l = 1; l < numLevels; ++l) { // Note: we start at level 1 (which is the first coarse level)
      for (int j = 0; j < maxRegPerProc; ++j) {
        RCP<MueLu::Level> level = regGrpHierarchy[j]->GetLevel(l);

        regProlong[l][j] = level->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
        regMatrices[l][j] = level->Get<RCP<Matrix> >("A", MueLu::NoFactory::get());

        regRowMaps[l][j] = Teuchos::rcp_const_cast<Map>(regMatrices[l][j]->getRowMap()); // ToDo (mayr.mt) Should we rather copy?
        regColMaps[l][j] = Teuchos::rcp_const_cast<Map>(regMatrices[l][j]->getColMap()); // ToDo (mayr.mt) Should we rather copy?
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
       * We need three row maps on the next coarser level:
       * - regional row map: just extract the RowMap() from the prolongator
       * - quasiRegional row map: needs to be formed manually
       * - composite row map: needs to be formed manually
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
       *
       * To form a composite operator on the coarsest level, we also need the
       * column map on the coarsest level. Hence, we have to recursively
       * re-construct all column maps on each level, namely
       * - regional col map: just extract the ColMap() from the prolongator
       * - quasiRegional col map: needs to be formed manually
       * - composite col map: needs to be formed manually
       *
       * To form column maps, lots of information that we used to generate row
       * maps can be reused. The central piece of information is the mapping
       * of duplicated interface GIDs.
       */

      int err = 0;

      ////////////////////////////////////////////////////////////////////////
      // IDENTIFY DUPLICATED GIDS ON FINE LEVEL
      ////////////////////////////////////////////////////////////////////////

      RCP<Matrix> summedDuplicateMatrix = Teuchos::null;

      // Find lists of duplicated GIDs locally
      interfaceLIDs.resize(numLevels);
      interfaceGIDPairs.resize(numLevels);
      for (int l = 0; l < numLevels; ++l) {
        interfaceLIDs[l].resize(maxRegPerProc);
        interfaceGIDPairs[l].resize(maxRegPerProc);
      }

      for (int l = 0; l < numLevels - 1; ++l) {
//        Comm->barrier();
//        if (Comm->getRank() == 0) {
//          std::cout << std::endl << std::endl
//              << "Processing GID pairs on level " << l
//              << std::endl << std::endl;
//        }
//        Comm->barrier();

//        sleep(1);
//        std::cout << "Prolongator" << std::endl;
//        regProlong[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);

        // create list of LIDs per group
        for (int j = 0; j < maxRegPerProc; j++) {
          for (size_t i = 0; i < regRowMaps[l][j]->getNodeNumElements(); ++i) {
            if (regRowMaps[l][j]->getGlobalElement(i) != quasiRegRowMaps[l][j]->getGlobalElement(i)) {
              // This is an interface node
              interfaceLIDs[l][j].push_back(i);
            }
          }
        }

//        for (int j = 0; j < maxRegPerProc; j++) {
//          std::cout << myRank << " | group = " << j << " | LIDs = ";
//          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
//            std::cout << interfaceLIDs[l][j][i] << ", ";
//          std::cout << std::endl;
//        }

        for (int j = 0; j < maxRegPerProc; j++)
          interfaceGIDPairs[l][j].resize(interfaceLIDs[l][j].size());

        // create list of LIDPairs per group
        for (int j = 0; j < maxRegPerProc; j++)
          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
            interfaceGIDPairs[l][j][i].push_back(quasiRegRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]));
        for (int j = 0; j < maxRegPerProc; j++)
          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
            if (regRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]) != quasiRegRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]))
              interfaceGIDPairs[l][j][i].push_back(regRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]));

//        std::cout << myRank << " | Print interfaceGIDPairs:" << std::endl;
//        for (int j = 0; j < maxRegPerProc; j++) {
//          std::cout << myRank << " | " << "Level " << l <<" | Group " << j << ":" << std::endl;
//          for (std::size_t i = 0; i < interfaceGIDPairs[l][j].size(); ++i) {
//            std::cout << "   " << myRank << " | ";
//            for (std::size_t k = 0; k < interfaceGIDPairs[l][j][i].size(); ++k)
//              std::cout << interfaceGIDPairs[l][j][i][k] << ", ";
//            std::cout << std::endl;
//          }
//        }

        std::vector<RCP<Vector> > regDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
        std::vector<RCP<Vector> > quasiRegDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
        for (int j = 0; j < maxRegPerProc; ++j) {
          regDupGIDVec[j] = VectorFactory::Build(regRowMaps[l][j], true);
          quasiRegDupGIDVec[j] = VectorFactory::Build(quasiRegRowMaps[l][j], true);

          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
            regDupGIDVec[j]->replaceLocalValue(interfaceLIDs[l][j][i], scalarOne);
        }

        RCP<Vector> compDupGIDVec = VectorFactory::Build(compRowMaps[l], true);
        regionalToComposite(regDupGIDVec, compDupGIDVec, maxRegPerProc,
            quasiRegRowMaps[l], regRowImporters[l], Xpetra::ADD);

        compositeToRegional(compDupGIDVec, quasiRegDupGIDVec, regDupGIDVec,
            maxRegPerProc, quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

        // create row/range/domain map for fine level duplicates mapping matrix
        RCP<Map> duplicateMap = Teuchos::null;
        {
          Teuchos::Array<GlobalOrdinal> myIntGIDs;
          for (int j = 0; j < maxRegPerProc; j++) {
            for (std::size_t i = 0; i < interfaceGIDPairs[l][j].size(); ++i) {
              for (std::size_t k = 0; k < interfaceGIDPairs[l][j][i].size(); ++k) {
                if (regRowMaps[l][j]->isNodeGlobalElement(interfaceGIDPairs[l][j][i][k]))
                  myIntGIDs.push_back(interfaceGIDPairs[l][j][i][k]);
              }
            }
          }

          duplicateMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
              myIntGIDs,Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);

//          sleep(1);
//          Comm->barrier();
//          std::cout << myRank << " | duplicateMap:" << std::endl;
//          duplicateMap->describe(*fos, Teuchos::VERB_HIGH);
        }

        // create row/range/domain map for the transpose of the fine level duplicates mapping matrix
        RCP<Map> fullDuplicateMap = Teuchos::null;
        {
          Array<GlobalOrdinal> myIntGIDs;
          ArrayRCP<const Scalar> regDupGIDVecData = regDupGIDVec[0]->getData(0);
          for (size_t i = 0; i < regRowMaps[l][0]->getNodeNumElements(); ++i) {
            if (regDupGIDVecData[i] != 0)
              myIntGIDs.push_back(Teuchos::as<GlobalOrdinal>(regRowMaps[l][0]->getGlobalElement(i)));
          }

//          std::cout << myRank << " | myIntGIDs = ";
//          for (auto gid : myIntGIDs)
//            std::cout << gid << ", ";
//          std::cout << std::endl;

          fullDuplicateMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
              myIntGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);

//          sleep(1);
//          Comm->barrier();
//          std::cout << myRank << " | fullDuplicateMap:" << std::endl;
//          fullDuplicateMap->describe(*fos, TEUCHOS::VERB_HIGH);
        }

        // Create and fill matrix
        RCP<Matrix> duplicateMatrix = MatrixFactory::Build(fullDuplicateMap, 2, Xpetra::DynamicProfile);
        {
          // Fill diagonal
          {
            Array<Scalar> vals(1);
            Array<GlobalOrdinal> cols(1);
            for (size_t i = 0; i < fullDuplicateMap->getNodeNumElements(); ++i) {
              GlobalOrdinal globalRow = fullDuplicateMap->getGlobalElement(i);

              vals[0] = static_cast<Scalar>(globalRow);
              cols[0] = globalRow;

              duplicateMatrix->insertGlobalValues(globalRow, cols, vals);
            }
          }

          // Fill off-diagonals (if known)
          Array<Scalar> vals(2);
          Array<GlobalOrdinal> cols(2);
          for (size_t i = 0; i < duplicateMap->getNodeNumElements(); ++i) {
            GlobalOrdinal globalRow = duplicateMap->getGlobalElement(i);

            for (std::size_t k = 0; k < interfaceGIDPairs[l][0][i].size(); ++k) {
              vals[k] = static_cast<Scalar>(globalRow);
              cols[k] = interfaceGIDPairs[l][0][i][k];
            }

            duplicateMatrix->insertGlobalValues(globalRow, cols, vals);
          }
        }

        duplicateMatrix->fillComplete();

//        sleep(1);
//        Comm->barrier();
//        std::cout << myRank << " | duplicateMatrix:" << std::endl;
//        duplicateMatrix->describe(*fos, Teuchos::VERB_HIGH);
//
//        sleep(1);
//        Comm->barrier();
//        std::cout << myRank << " | duplicateMatrix->RangeMap():" << std::endl;
//        duplicateMatrix->getRangeMap()->describe(*fos, Teuchos::VERB_HIGH);
//
//        sleep(1);
//        Comm->barrier();
//        std::cout << myRank << " | duplicateMatrix->DomainMap():" << std::endl;
//        duplicateMatrix->getDomainMap()->describe(*fos, Teuchos::VERB_HIGH);
//
//        sleep(1);
//        Comm->barrier();
//        std::cout << myRank << " | duplicateMatrix->ColMap():" << std::endl;
//        duplicateMatrix->getColMap()->describe(*fos, Teuchos::VERB_HIGH);

        {
          RCP<Matrix> transDuplicateMatrix = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Transpose(*duplicateMatrix);
          TEUCHOS_ASSERT(!transDuplicateMatrix.is_null());

          transDuplicateMatrix->fillComplete();
          TEUCHOS_ASSERT(transDuplicateMatrix->isFillComplete());

//          sleep(1);
//          Comm->barrier();
//          std::cout << myRank << " | transDuplicateMatrix:" << std::endl;
//          transDuplicateMatrix->describe(*fos, Teuchos::VERB_HIGH);

          summedDuplicateMatrix = MatrixFactory::Build(fullDuplicateMap, 2, Xpetra::DynamicProfile);
          MatrixMatrix::TwoMatrixAdd(*transDuplicateMatrix, false, scalarOne, *duplicateMatrix, false, scalarOne, summedDuplicateMatrix, *fos);
        }

        summedDuplicateMatrix->fillComplete(fullDuplicateMap, fullDuplicateMap);
//        summedDuplicateMatrix->fillComplete();

//        sleep(1);
//        Comm->barrier();
//        std::cout << myRank << " | summedDuplicateMatrix:" << std::endl;
//        summedDuplicateMatrix->describe(*fos, Teuchos::VERB_HIGH);

        std::vector<std::vector<GlobalOrdinal> > myDuplicates; // pairs of duplicates with locally owned GID listed first.
        myDuplicates.resize(summedDuplicateMatrix->getNodeNumRows());
        for (size_t lRow = 0; lRow < summedDuplicateMatrix->getNodeNumRows(); ++lRow) {

          ArrayView<const LocalOrdinal> inds;
          ArrayView<const Scalar> vals;

          summedDuplicateMatrix->getLocalRowView(lRow, inds, vals);

          std::vector<GlobalOrdinal> gidsToSort;
          for (ArrayView<const LocalOrdinal>::size_type k = 0; k < inds.size(); ++k) {
//            std::cout << myRank << " | inds[" << k << "] = " << inds[k] << std::endl;
            gidsToSort.push_back(summedDuplicateMatrix->getColMap()->getGlobalElement(inds[k]));
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

            myDuplicates[lRow] = gidsToSort;
          }
        }

//        sleep(myRank);
//        std::cout << std::endl << std::endl << std::endl << std::endl;
//        for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
//          std::cout << myRank << " | myDuplicates:" << std::endl << "  ";
//          for (std::size_t k = 0; k < myDuplicates[i].size(); ++k)
//            std::cout << myDuplicates[i][k] << ", ";
//          std::cout << std::endl;
//        }
//        std::cout << std::endl << std::endl << std::endl << std::endl;

        ////////////////////////////////////////////////////////////////////////
        // IDENTIFY DUPLICATED GIDS ON COARSE LEVEL
        ////////////////////////////////////////////////////////////////////////

        std::vector<std::vector<GlobalOrdinal> > myCoarseDuplicates; // pairs of duplicates with locally owned GID listed first on coarse level.
        myCoarseDuplicates.resize(myDuplicates.size());

        std::vector<GlobalOrdinal> myCoarseInterfaceGIDs;
        for (size_t i = 0; i < myDuplicates.size(); ++i) {
          const GlobalOrdinal rowGID = myDuplicates[i][0];
          const LocalOrdinal rowLID = regProlong[l+1][0]->getRowMap()->getLocalElement(rowGID);

          ArrayView<const Scalar> vals;
          ArrayView<const LocalOrdinal> inds;
          regProlong[l+1][0]->getLocalRowView(rowLID, inds, vals);
//          std::cout << myRank << " | ExtractMyRowView err = " << err << std::endl;
          TEUCHOS_ASSERT(inds.size() == 1); // tentative P: only one entry per row
          TEUCHOS_ASSERT(vals.size() == 1); // tentative P: only one entry per row

          myCoarseInterfaceGIDs.push_back(regProlong[l+1][0]->getColMap()->getGlobalElement(inds[0]));
        }

//        std::cout << "myCoarseInterfaceGIDs on proc " << myRank << ": ";
//        for (auto gid : myCoarseInterfaceGIDs) {
//          std::cout << gid << ", ";
//        }
//        std::cout << std::endl;

        // Build row/range/domain map of duplicate mapping matrix
        RCP<Map> interfaceMap = Teuchos::null;
        {
          std::vector<GlobalOrdinal> interfaceMapGIDs(myDuplicates.size());
          for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
            interfaceMapGIDs[i] = myDuplicates[i][0];
          }

//          std::cout << "interfaceMapGIDs on proc " << myRank << ":      ";
//          for (auto gid : interfaceMapGIDs) {
//            std::cout << gid << ", ";
//          }
//          std::cout << std::endl;

          interfaceMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
              interfaceMapGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
        }

        RCP<Matrix> duplicateMapping = MatrixFactory::Build(interfaceMap, 2, Xpetra::DynamicProfile);

        // Fill the matrix
        RCP<Matrix> transDuplicateMapping = Teuchos::null;
        {
          size_t numRows = myDuplicates.size();
          std::vector<GlobalOrdinal> rowPtr(numRows);
          std::vector<Array<Scalar> > vals(numRows);
          std::vector<Array<LocalOrdinal> > colInds(numRows);

          for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
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
          for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
            duplicateMapping->insertGlobalValues(rowPtr[rowIdx], colInds[rowIdx], vals[rowIdx]);
          }

          duplicateMapping->fillComplete();
          TEUCHOS_ASSERT(duplicateMapping->isFillComplete());

//          sleep(1);
//          Comm.Barrier();
//          duplicateMapping->Print(std::cout);

          transDuplicateMapping = Utilities::Transpose(*duplicateMapping);

          Comm->barrier();

          TEUCHOS_ASSERT(!transDuplicateMapping.is_null());

          transDuplicateMapping->fillComplete();
          TEUCHOS_ASSERT(transDuplicateMapping->isFillComplete());

//          sleep(1);
//          std::cout << myRank << " | Printing the transpose ..." << std::endl;
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
        std::vector<std::vector<GlobalOrdinal> > myCoarseInterfaceDuplicates(transDuplicateMapping->getNodeNumRows());
        {
          for (size_t lRow = 0; lRow < transDuplicateMapping->getNodeNumRows(); ++lRow) {
            ArrayView<const Scalar> vals;
            ArrayView<const LocalOrdinal> inds;
            transDuplicateMapping->getLocalRowView(lRow, inds, vals);

            const size_t numEntries = inds.size();
            myCoarseInterfaceDuplicates[lRow].resize(numEntries);
            myCoarseInterfaceDuplicates[lRow].assign(vals.begin(), vals.end());

//            std::cout << myRank << " | myCoarseInterfaceDuplicates[" << i << "] = ";
//            for (auto id : myCoarseInterfaceDuplicates[i])
//              std::cout << id << ", ";
//            std::cout << std::endl;
          }

        }

        ////////////////////////////////////////////////////////////////////////
        // CREATE COARSE LEVEL MAPS
        ////////////////////////////////////////////////////////////////////////
        {
//          sleep(1);
//          std::cout << myRank << " | Printing regRowMaps[" << l+1 << "][0] ..." << std::endl;
//          Comm->barrier();
//          regRowMaps[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);
//          sleep(2);

          // create quasiRegional row map
          {
            std::vector<GlobalOrdinal> myQuasiRegGIDs;

            for (size_t i = 0; i < regRowMaps[l+1][0]->getNodeNumElements(); ++i) {
              // grab current regional GID to be processed
              GlobalOrdinal currGID = regRowMaps[l+1][0]->getGlobalElement(i);
              GlobalOrdinal quasiGID = currGID; // assign dummy value

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
                        quasiGID = std::min(quasiGID, Teuchos::as<const GlobalOrdinal>(myCoarseInterfaceDuplicates[k][kk]));

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

            quasiRegRowMaps[l+1][0] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                myQuasiRegGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
            TEUCHOS_ASSERT(!quasiRegRowMaps[l+1][0].is_null());

//            sleep(1);
//            std::cout << myRank << " | Printing quasiRegRowMaps[" << l+1 << "][0] ..." << std::endl;
//            Comm->barrier();
//            quasiRegRowMaps[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);
          }

          // create composite row map
          {
            std::vector<GlobalOrdinal> myCompGIDs;
            for (int j = 0; j < maxRegPerProc; j++) {
              for (size_t i = 0; i < quasiRegRowMaps[l+1][j]->getNodeNumElements(); ++i) {
                const GlobalOrdinal trialGID = quasiRegRowMaps[l+1][j]->getGlobalElement(i);

                if (regRowMaps[l+1][j]->isNodeGlobalElement(trialGID))
                  myCompGIDs.push_back(trialGID);
              }
            }

            compRowMaps[l+1] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                myCompGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
            TEUCHOS_ASSERT(!compRowMaps[l+1].is_null());

//            sleep(1);
//            std::cout << myRank << " | Printing compRowMaps["<< l+1 << "] ..." << std::endl;
//            Comm->barrier();
//            compRowMaps[l+1]->describe(*fos, Teuchos::VERB_HIGH);
          }

          // create regRowImporter
          for (int j = 0; j < maxRegPerProc; ++j) {
            regRowImporters[l+1][j] = ImportFactory::Build(compRowMaps[l+1], quasiRegRowMaps[l+1][j]);
            TEUCHOS_ASSERT(!regRowImporters[l+1][j].is_null());
          }

          // Create quasiRegional column map
          {
            std::vector<GlobalOrdinal> myQuasiRegGIDs;

            for (size_t i = 0; i < regColMaps[l+1][0]->getNodeNumElements(); ++i) {
              // grab current regional GID to be processed
              GlobalOrdinal currGID = regColMaps[l+1][0]->getGlobalElement(i);
              GlobalOrdinal quasiGID = currGID; // assign dummy value

              /* Find quasiRegional counterpart
               *
               * We should be able to use the same procedure as for the
               * quasiRegRowMap since duplicated GIDs only live on the interface.
               * Hence, even if the column map reaches across an interface,
               * GIDs from 'the other side of the interface' do not have to be
               * modified as they have not been duplicated (only those on the
               * interface).
               */
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
                        quasiGID = std::min(quasiGID, Teuchos::as<GlobalOrdinal>(myCoarseInterfaceDuplicates[k][kk]));

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

            quasiRegColMaps[l+1][0] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                myQuasiRegGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
            TEUCHOS_ASSERT(!quasiRegColMaps[l+1][0].is_null());

//            sleep(1);
//            std::cout << myRank << " | Printing quasiRegColMaps[" << l+1 << "][0] ..." << std::endl;
//            Comm->barrier();
//            quasiRegColMaps[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);
          }
        }
      }
    }

    Comm->barrier();

    // Form the composite coarse level operator
    {
      const int maxLevel = numLevels - 1;
      coarseCompOp = MatrixFactory::Build(compRowMaps[maxLevel], 3, Xpetra::DynamicProfile);
//      coarseCompOp->setAllToScalar(scalarZero);
//      coarseCompOp->describe(*fos, Teuchos::VERB_EXTREME);

      regionalToComposite(regMatrices[maxLevel], coarseCompOp, maxRegPerProc,
          quasiRegRowMaps[maxLevel], quasiRegColMaps[maxLevel],
          regRowImporters[maxLevel], Xpetra::ADD);

//      coarseCompOp->fillComplete(compRowMaps[maxLevel], compRowMaps[maxLevel]);
//      TEUCHOS_ASSERT(coarseCompOp->isFillComplete());
//
//      sleep(1);
//      std::cout << myRank << " | Printing coarseCompOp ..." << std::endl;
//      Comm->barrier();
//      coarseCompOp->describe(*fos, Teuchos::VERB_HIGH);
    }
  }

  Comm->barrier();

  // Make interface scaling factors recursively
  {
    std::cout << myRank << " | Computing interface scaling factors ..." << std::endl;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

    for (int l = 0; l < numLevels; l++)
    {
      // initialize region vector with all ones.
      for (int j = 0; j < maxRegPerProc; j++) {
        regInterfaceScalings[l][j] = VectorFactory::Build(regRowMaps[l][j]);
        regInterfaceScalings[l][j]->putScalar(scalarOne);
      }

      // transform to composite layout while adding interface values via the Export() combine mode
      RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(compRowMaps[l], true);
      regionalToComposite(regInterfaceScalings[l], compInterfaceScalingSum, maxRegPerProc, quasiRegRowMaps[l], regRowImporters[l], Xpetra::ADD);

      /* transform composite layout back to regional layout. Now, GIDs associated
       * with region interface should carry a scaling factor (!= 1).
       */
      std::vector<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc);
      compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
          regInterfaceScalings[l], maxRegPerProc, quasiRegRowMaps[l],
          regRowMaps[l], regRowImporters[l]);
    }
  }

  Comm->barrier();

  // Run V-cycle
  {
    std::cout << myRank << " | Running V-cycle ..." << std::endl;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

    /* We first use the non-level container variables to setup the fine grid problem.
     * This is ok since the initial setup just mimics the application and the outer
     * Krylov method.
     *
     * We switch to using the level container variables as soon as we enter the
     * recursive part of the algorithm.
     */

    // initial guess for solution
    compX = VectorFactory::Build(mapComp, true);

// -- NOT MIGRATED TO XPETRA YET -- START
// debugging using shadow.m
//double *z;
//compX->ExtractView(&z); for (int kk = 0; kk < compX->MyLength(); kk++) z[kk] = (double) kk*kk;
//compX->Print(std::cout);
// -- NOT MIGRATED TO XPETRA YET -- END
    // forcing vector
    RCP<Vector> compB = VectorFactory::Build(mapComp, true);

    if (doing1D)
    {
      // 1D
      {
        {
          compB->replaceGlobalValue(compB->getGlobalLength() - 1, 0, Teuchos::as<Scalar>(1.0e-3));
        }
//        {
//        compB->putScalar(scalarOne);
//        compB->replaceGlobalValue(Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), scalarOne);
//        }
//        {
//          compB->replaceGlobalValue(Teuchos::as<GlobalOrdinal>(15), scalarOne);
//        }
//        {
//          compB->replaceGlobalValue(Teuchos::as<GlobalOrdinal>(16), scalarOne);
//        }
      }
    }
    else //2D
    {
      std::vector<GlobalOrdinal> dirichletGIDs;

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
        compB->putScalar(Teuchos::as<Scalar>(1.0e-3));
        for (size_t i = 0; i < dirichletGIDs.size(); ++i)
          compB->replaceGlobalValue(dirichletGIDs[i], scalarZero);
      }
    }

    // residual vector
    RCP<Vector> compRes = VectorFactory::Build(mapComp, true);
    {
      AComp->apply(*compX, *compRes);
      compRes->update(1.0, *compB, -1.0);
    }

    // transform composite vectors to regional layout
    compositeToRegional(compX, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
        revisedRowMapPerGrp, rowImportPerGrp);

    std::vector<RCP<Vector> > quasiRegB(maxRegPerProc);
    std::vector<RCP<Vector> > regB(maxRegPerProc);
    compositeToRegional(compB, quasiRegB, regB, maxRegPerProc, rowMapPerGrp,
        revisedRowMapPerGrp, rowImportPerGrp);

//    printRegionalObject<Vector>("regB 0", regB, myRank, *fos);

    std::vector<RCP<Vector> > regRes(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) { // step 1
      regRes[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    }

    /////////////////////////////////////////////////////////////////////////
    // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
    /////////////////////////////////////////////////////////////////////////

    // define max iteration counts
    const int maxVCycle = 200;
    const int maxFineIter = 20;
    const int maxCoarseIter = 100;
    const double omega = 0.67;

    // Prepare output of residual norm to file
    RCP<std::ofstream> log;
    if (myRank == 0)
    {
      std::string s = "residual_norm.txt";
      log = rcp(new std::ofstream(s.c_str()));
      (*log) << "# num procs = " << Comm->getSize() << "\n"
             << "# iteration | res-norm\n"
             << "#\n";
    }

    // Richardson iterations
    for (int cycle = 0; cycle < maxVCycle; ++cycle) {
      // check for convergence
      {
        ////////////////////////////////////////////////////////////////////////
        // SWITCH BACK TO NON-LEVEL VARIABLES
        ////////////////////////////////////////////////////////////////////////

        computeResidual(regRes, regX, regB, regionGrpMats, mapComp,
            rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

//        printRegionalObject<Vector>("regB 1", regB, myRank, *fos);

        compRes = VectorFactory::Build(mapComp, true);
        regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
            rowImportPerGrp, Xpetra::ADD);
        Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();

        // Output current residual norm to screen (on proc 0 only)
        if (myRank == 0)
        {
          std::cout << cycle << "\t" << normRes << std::endl;
          (*log) << cycle << "\t" << normRes << "\n";
        }

        if (normRes < 1.0e-12)
          break;
      }

      /////////////////////////////////////////////////////////////////////////
      // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
      /////////////////////////////////////////////////////////////////////////

//      printRegionalObject<Vector>("regB 2", regB, myRank, *fos);
      vCycle(0, numLevels, maxFineIter, maxCoarseIter, omega, maxRegPerProc, regX, regB, regMatrices,
          regProlong, compRowMaps, quasiRegRowMaps, regRowMaps, regRowImporters,
          regInterfaceScalings, coarseCompOp);
    }

    ////////////////////////////////////////////////////////////////////////
    // SWITCH BACK TO NON-LEVEL VARIABLES
    ////////////////////////////////////////////////////////////////////////

    // -----------------------------------------------------------------------
    // Print fine-level solution
    // -----------------------------------------------------------------------
/*
    Comm->barrier();
    sleep(1);

    // ToDo (mayr.mt) Is this the right CombineMode?
    regionalToComposite(regX, compX, maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Xpetra::INSERT);

    std::cout << myRank << " | compX after V-cycle" << std::endl;
    sleep(1);
    compX->describe(*fos, Teuchos::VERB_HIGH);
    sleep(2);

    // Write solution to file for printing
    std::string outFileName = "compX.mm";
    Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write(outFileName, *compX);
*/
  }

  Comm->barrier();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

}

/* Returns local ID (within region curRegion) for the LIDcomp^th composite grid
 * point that myRank owns. If this grid point is not part of curRegion, then
 * -1 is returned.
 *
 * For fully structured 1D problems only.
 */
LocalOrdinal LIDregion(void *ptr, LocalOrdinal LIDcomp, int whichGrp)
{
   struct widget * myWidget = (struct widget *) ptr;

   int        *minGIDComp  = myWidget->minGIDComp;
   int        *maxGIDComp  = myWidget->maxGIDComp;
   int        *myRegions   = myWidget->myRegions;
   Map        *colMap      = myWidget->colMap;
   int        maxRegPerGID = myWidget->maxRegPerGID;
   Teuchos::RCP<MultiVector> regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

   int curRegion = myRegions[whichGrp];
   if (LIDcomp >=  Teuchos::as<int>(colMap->getNodeNumElements())) return(-1);

   const Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = colMap->getNodeElementList();

   if (colGIDsComp[LIDcomp] < minGIDComp[whichGrp]) return(-1);
   if (colGIDsComp[LIDcomp] > maxGIDComp[whichGrp]) return(-1);

   bool found = false;
   for (int j = 0; j <  maxRegPerGID; j++) {
     Teuchos::ArrayRCP<const Scalar> jthRegions = regionsPerGIDWithGhosts->getData(j);
      if  ( ((int) jthRegions[LIDcomp]) == curRegion ) {
         found = true;
         break;
      }
   }
   if (found == false) return(-1);

   return( colGIDsComp[LIDcomp] - minGIDComp[whichGrp] );
}

/* Returns local ID (within region curRegion) for the LIDcomp^th composite grid
 * point that myRank owns. If this grid point is not part of curRegion, then
 * -1 is returned.
 *
 * For fully structured 2D problems only.
 */
LocalOrdinal LID2Dregion(void *ptr, LocalOrdinal LIDcomp, int whichGrp)
{
   struct widget * myWidget = (struct widget *) ptr;

//   int        *minGIDComp  = myWidget->minGIDComp;
//   int        *maxGIDComp  = myWidget->maxGIDComp;
   int        *myRegions   = myWidget->myRegions;
   Map        *colMap      = myWidget->colMap;
   int        maxRegPerGID = myWidget->maxRegPerGID;
   int       *trueCornerx  = myWidget->trueCornerx; // global coords of region
   int       *trueCornery  = myWidget->trueCornery; // corner within entire 2D mesh
                                                    // across all regions/procs
   int       *relcornerx   = myWidget->relcornerx;  // coords of corner relative
   int       *relcornery   = myWidget->relcornery;  // to region corner
   int       *lDimx        = myWidget->lDimx;
   int       *lDimy        = myWidget->lDimy;
   Teuchos::RCP<MultiVector> regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

   int curRegion = myRegions[whichGrp];
   if (LIDcomp >=  Teuchos::as<int>(colMap->getNodeNumElements())) return(-1);

   const Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = colMap->getNodeElementList();

   int xGIDComp =  colGIDsComp[LIDcomp]%(myWidget->nx);
   int yGIDComp = (colGIDsComp[LIDcomp] - xGIDComp)/myWidget->nx;

   if (xGIDComp < trueCornerx[whichGrp]+relcornerx[whichGrp]) return(-1);
   if (yGIDComp < trueCornery[whichGrp]+relcornery[whichGrp]) return(-1);
   if (xGIDComp > trueCornerx[whichGrp]+relcornerx[whichGrp]+lDimx[whichGrp]-1) return(-1);
   if (yGIDComp > trueCornery[whichGrp]+relcornery[whichGrp]+lDimy[whichGrp]-1) return(-1);


   bool found = false;
   for (int j = 0; j <  maxRegPerGID; j++) {
     Teuchos::ArrayRCP<const Scalar> jthRegions = regionsPerGIDWithGhosts->getData(j);
      if  ( ((int) jthRegions[LIDcomp]) == curRegion ) {
         found = true;
         break;
      }
   }
   if (found == false) return(-1);

   return(
    (yGIDComp - relcornery[whichGrp]-trueCornery[whichGrp])*lDimx[whichGrp]+ (xGIDComp - relcornerx[whichGrp]-trueCornerx[whichGrp]));
}



void fillCircleSquareData(int ArowPtr[], int Acols[], int ownedX, int ownedY, int Cx, int Cy, int Rx, int Ry,
std::vector<int> &appData)
{
/* Fills the vector appData so that the function LIDregionCircleSquareWithUnstr2D() can
 * properly map CompositeLIDs to RegionLIDs. Specifically, appData's
 * contents will be given by
 *
 *   appData[inpData_ownedX]  x/y dimensions of rectangle owned by proc, which is
 *   appData[inpData_ownedY]  given as input parameters 'ownedX' and 'ownedY'
 *
 *   appData[inpData_regionX] x/y dimensions of proc's region, which is given
 *   appData[inpData_regionY] as input parameters 'Rx' and 'Ry'
 *
 *   appData[inpData_cornerX] Offset of the lower left corner defined by the
 *   appData[inpData_cornerY] rectangular region piece that is actually owned by
 *                            this processor (in the composite layout). Should be
 *                            either 0 or 1.  So, Cx = Cy=0 means that the processor
 *                            actually owns the lower left corner of the region.
 *                            This is given as input parameters 'Cx' and 'Cy'
 *
 *   appData[k]           Gives the region LID associated with the
 *                        (k-inpData_firstLIDsOfGhosts+ownedX*ownedY)^th
 *                        composite LID  .. for k >= inpData_firstLIDsOfGhosts,
 *                        which is the (k-inpData_firstLIDsOfGhosts)^th ghost
 *                        composite LID.
 *
 *
 * Before filling appData, fillCircleSquareData() first fills
 * ghostCompLIDs with ids of any ghosts associated with a region boundary.
 * The edges are done in the following order: bottom, top, left, right.
 * If a region boundary does not correspond to ghost unknowns in the composite layout,
 * then this boundary is not included in ghostCompLIDs. Once filled, ghostCompLIDs[]
 * should have the following contents:
 *
 *    ghostCompLIDs[startBot:startTop-1]   Ghost composite LIDs for bottom
 *                                         edge of region
 *
 *    ghostCompLIDs[startTop:startLft-1]   Ghost composite LIDs for top
 *                                         edge of region
 *
 *    ghostCompLIDs[startLft:startRgt-1]   Ghost composite LIDs for left edge
 *                                         of region. Does not include corners
 *                                         if they are already included
 *                                         with the bottom or top ghosts
 *
 *    ghostCompLIDs[startRgt:              Ghost composite LIDs for right edge
 *                     nRegionalGhosts ]   of region. Does not include corners
 *                                         if they are already included
 *                                         with the bottom or top ghosts
 *
 * Read comments for edgeGhosts(). The main assumption is that all stencils are full (9
 * point in the interior).
 *
 */

  /* Compute the number of ghosts that lie on the region boundary and where*/
  /* the different region boundaries will start in the vector.             */

  int nRegionalGhosts = 0;
  int startBot = nRegionalGhosts;
  if (Cy==1) {
     nRegionalGhosts= nRegionalGhosts+ownedX;
     if (Cx==1)           nRegionalGhosts=nRegionalGhosts+1;
     if (Rx==ownedX+Cx+1) nRegionalGhosts=nRegionalGhosts+1;
  }
  int startTop = nRegionalGhosts;

  if (Ry == ownedY+Cy+1) {
     nRegionalGhosts= nRegionalGhosts+ownedX;
     if (Cx==1)          nRegionalGhosts=nRegionalGhosts+1;
     if (Rx==ownedX+Cx+1)nRegionalGhosts=nRegionalGhosts+1;
  }
  int startLft = nRegionalGhosts;

  if (Cx==1)             nRegionalGhosts= nRegionalGhosts+ownedY;
  int startRgt = nRegionalGhosts;

  if (Rx == ownedX+Cx+1) nRegionalGhosts= nRegionalGhosts+ownedY;
  std::vector<int> ghostCompLIDs(nRegionalGhosts);

  /* insert ghosts for bottom, top, left, and right edges into ghostCompLIDs */

  int nGhostFound = 0;
  if (Cy==1)            edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(), Rx-Cx-1,1,     1,       2,2-Cx,ownedX, ownedY);
  if (Ry == ownedY+Cy+1)edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(), Rx-Cx-1,1,ownedY,ownedY-1,2-Cx,ownedX, ownedY);
  if (Cx==1)            edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(),ownedY-1,0,     1,       2,   2,ownedX, ownedY);
  if (Rx == ownedX+Cx+1)edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(),ownedY-1,0,ownedX,ownedX-1,   2,ownedX, ownedY);

  /* determine the largest ghost LID so that we can allocate enough space */
  int biggest = ownedX*ownedY-1;
  for (int k = 0; k < nRegionalGhosts; k++)
     if (ghostCompLIDs[k] > biggest) biggest = ghostCompLIDs[k];

  // fill appData

  appData.resize(inpData_firstLIDsOfGhosts+biggest-ownedX*ownedY);

  appData[inpData_isStructured] = 1;
  appData[inpData_ownedX] = ownedX;
  appData[inpData_ownedY] = ownedY;
  appData[inpData_regionX] = Rx;
  appData[inpData_regionY] = Ry;
  appData[inpData_cornerX] = Cx;
  appData[inpData_cornerY] = Cy;
  appData[inpData_nGhosts] = biggest-ownedX*ownedY;

  int offset = inpData_firstLIDsOfGhosts-ownedX*ownedY;

  for (int k = 0; k < biggest-ownedX*ownedY; k++) appData[inpData_firstLIDsOfGhosts+k] = -1;
  for (int k=startBot; k < startTop; k++)
     appData[ghostCompLIDs[k]+offset]=k;
  for (int k=startTop; k < startLft; k++)
     appData[ghostCompLIDs[k]+offset] = Rx*(Ry-1)+k-startTop;
  for (int k=startLft; k < startRgt; k++)
     appData[ghostCompLIDs[k]+offset] = Rx*(k+Cy-startLft);
  for (int k=startRgt; k < nRegionalGhosts; k++)
     appData[ghostCompLIDs[k]+offset] = Rx*(k+Cy-startRgt)+Rx-1;

//  for (int i = 0; i < ghostCompLIDs.size(); i++)
//    printf("ghostComp(%d)=%d ", i, ghostCompLIDs[i]);
//  printf("\n");
//  fflush (stdout);

  return;
}

void edgeGhosts(int ArowPtr[], int Acols[], int &nGhostFound, int ghostCompLIDs[], int edgeLength, int alongX, int ownedEdge, int interiorEdge, int start, int ownedX, int ownedY)
{
/*
%
%  Find the local region-oriented ids of a region's shared interface, which is
%  owned by another processor in the composite layout. The situation is basically
%  this
%
%            ?   ?   ?  ?  ?  ?  ?  ?
%             ======================
%            ?|| .   .  .  .  .  .||?
%            ?||                  ||?
%
%  The = and || denote the inter-processor boundary. We know that we have a
%  bunch of externally owned vertices, denoted by ?. The problem is that we
%  are not completely sure which remote local composite id is associated with
%  which vertex on the regular grid layout (as this was done automatically
%  by a fill complete on the composite matrix. To figure this out, this function
%  assumes a 9-pt stencil and that the chunk that each processor owns is at
%  least 3 wide in each dimension.
%
%  The way local region ids are computed is as follows:
%   1) We basically do a series of find(A(row,:)) for matrix rows corresponding
%      to the ownedEdge that is adjacent to the shared interface that we wish
%      to find region ids.
%   2) We assign any non-local columns (those > nOwned) an index, ii, indicating
%      that this column is adjacent to the ii^th point along ownedEdge. Some
%      columns might be adjacent to several points along ownedEdge. These
%      columns end up with the largest ii value, as we over-write the assigned
%      indices. Thus, a typical situation might look like this
%
%            1   2   3  4  5  6  6  6
%             ======================
%            1|| .   .  .  .  .  .||6
%            1||                  ||6
%
%      In this picture, locations corresponding to entries owned by other
%      processors have been assigned a number. The = and || symbols denote
%      the edge of the inter-processor boundary. The .'s lie along ownedEdge
%      So non-local dofs of the kth dot (counting from left to right) are
%      assigned the number k.
%   3) The two lower 1's and the two lower 6's are problematic as we do not
%      wish to consider these vertices as each function invocation is only
%      concerned with the one region edge. So we will stick large numbers
%      (equal to nOwned) in these locations. This is done by probing points
%      that are 1 away (provided externally as interiorEdge) from the corner
%      along orthogonal edges and assigning these a largeNumber. Our example
%      might now look like
%
%            1   2   3  4  5  6  6  6
%             ======================
%            L|| .   .  .  .  .  .||L
%            L|| *               *||L
%            L||                  ||L
%
%      where L denotes newly assigned large numbers and * are the one away
%      point that were just probed to assign large numbers.
%
%   4) We can now determine which remote vertices correspond to the 2D layout.
%      Specifically, we take the 1st vertex along ownedEdge and examine its
%      5 neighbors (e.g., col1,col2,col3,col4,col5). We look at assignedIndices
%      computed in steps 2 and 3. The column with the lowest assignedIndices
%      value (= 1) is the corner point. The column with the next lowest
%      (= 2) is adjacent to this corner point along the desired edge. The
%      column with the next lowest value (= 3) is the next point. We record
%      these 3 column indices in ghostCompLIDs and set assignedIndices for
%      them to a large number (so now 1 2 3 in the above picture would now be
%      replaced by L L L). We now examine the 2nd vertex's remote neighbors
%      which should have the values L L 3 and assign the column associated
%      with the smallest value (= 3) to ghostCompLIDs ... continuing along
%      until we have assigned the entire edge.
%
%  Note: In the above picture we are assuming that all region edges are owned
%  by remote processors. However, this function works as well when either
%  of the orthogonal edges (vertical edges in our example) are owned locally.
%  For example, we might have the case below
%
%                ?   ?  ?  ?  ?  ?  ?
%             ======================
%             || .   .  .  .  .  .||?
%             ||                  ||?
%  Here, the left edge might be a real physical boundary or it might be that
%  the processor owns the shared interface of the region (so it has remotes
%  to the left of the leftmost ||, but there is no need to assign them a
%  local region-oriented id.  In these cases start is not necessarily 1
%  and fullEdgeLength is not necessarily equal to edgeLength
*/

  if (ownedX < 3) { fprintf(stderr,"edges must be longer\n"); exit(1);}
  if (ownedY < 3) { fprintf(stderr,"edges must be longer\n"); exit(1);}

  int  nOwned      = ownedX*ownedY;
  int  largeNumber = nOwned+10;

  std::vector<int> assignedIndices( (ownedX+2)*(ownedY+2) + 1);

  /* Perform steps 1) and 2) described above. */

  int fullEdgeLength, row, *cols, nCols;

  if (alongX)  fullEdgeLength = ownedX;
  else         fullEdgeLength = ownedY;

  for (int ii=1; ii <= fullEdgeLength; ii++) {

    /* non-Ghosts are lexicographically ordered  */

    if   (alongX) row = ownedX*(ownedEdge-1) +     ii   ;
    else          row = ownedX*(   ii    -1) + ownedEdge;
    cols = &(Acols[ArowPtr[row-1]]);
    nCols = ArowPtr[row] - ArowPtr[row-1];
    for (int k = 0; k < nCols; k++) {
       if (cols[k] >= nOwned ) assignedIndices[cols[k]] = ii;
    }
  }

  /* Now assign large numbers (step 3) to 2 closest vertices */
  /* along each orthogonal edge                              */

  /* non-Ghosts are lexicographically ordered  */
  if  (alongX) row = ownedX*(interiorEdge-1)+     1     ;
  else         row =                        interiorEdge;

  cols = &(Acols[ArowPtr[row-1]]);
  nCols = ArowPtr[row] - ArowPtr[row-1];
  bool firstCornerHasGhosts = false;
  for (int k = 0; k < nCols; k++)
    if (cols[k] >= nOwned) {firstCornerHasGhosts = true; assignedIndices[cols[k]]=largeNumber;}

// This is a case not originally considered. I hope it fixes this bug.
// When coding this up, I failed to recognize the case when Cx is 0
// because
if ( (start==2) && (firstCornerHasGhosts == false)) start--;

  /* non-Ghosts are lexicographically ordered  */
  if  (alongX) row = ownedX*(interiorEdge-1) +  edgeLength ;
  else         row = ownedX*( edgeLength -1) + interiorEdge;

  cols = &(Acols[ArowPtr[row-1]]);
  nCols = ArowPtr[row] - ArowPtr[row-1];
  for (int k = 0; k < nCols; k++)
    if (cols[k] >= nOwned) assignedIndices[cols[k]]=largeNumber;

  /* Now perform step 4 by looking at smallest */
  /* assignedIndices along ownedEdge.          */

  int min, kk;

  for (int ii=1; ii <= edgeLength; ii++) {

    /* non-Ghosts are lexicographically ordered  */
    if  (alongX) row = ownedX*(ownedEdge-1) +    ii    ;
    else         row = ownedX*(   ii    -1) + ownedEdge;

    cols = &(Acols[ArowPtr[row-1]]);
    nCols= ArowPtr[row] - ArowPtr[row-1];

    min  = largeNumber-1;  kk = -1;

    for (int k = 0; k < nCols; k++) {
      if (cols[k] >= nOwned ) {
        if (assignedIndices[cols[k]] == min) { fprintf(stderr,"a tie?\n"); exit(1);}
        if (assignedIndices[cols[k]] < min) {
          kk = cols[k];
          min = assignedIndices[cols[k]];
        }
      }
    }
    if ((ii>=start) && (kk != -1)) ghostCompLIDs[nGhostFound++]= kk;
    if (kk != -1) assignedIndices[kk] = largeNumber;

    if (ii==1) {
      for (int kkkk = 1; kkkk <= 2; kkkk++) {
        min  = largeNumber-1;  kk = -1;
        for (int k = 0; k < nCols; k++) {
          if (cols[k] >= nOwned ) {
            if (assignedIndices[cols[k]] == min) { fprintf(stderr,"a Tie?\n"); exit(1);}
            if (assignedIndices[cols[k]] < min) {
              kk = cols[k];
              min = assignedIndices[cols[k]];
            }
          }
        }
        if (kk != -1) {
          ghostCompLIDs[nGhostFound++]= kk;
          assignedIndices[kk] = largeNumber;
        }
      } // for (int kkkk = 1; kkkk <= 2; kkkk++)
    } // if (ii==1)
  } // for (int ii=1; ii <= edgeLength; ii++)
}

LocalOrdinal LIDregionCircleSquare(void *ptr, LocalOrdinal compLID, int whichGrp)
{
   // Maps composite LIDs to region LIDs for the example
   // corresponding to a circle embedded within a box created
   // by the mkUnstrQuads Matlab functions. Assumes that
   // fillCircleSquareData() has been invoked.

   int* appData      = (int *) ptr;

   if (appData[inpData_isStructured] == 0) return(compLID);

   int* LIDsOfGhosts = (int *) &(appData[inpData_firstLIDsOfGhosts]);
   int  ownedX = appData[inpData_ownedX];
   int  ownedY = appData[inpData_ownedY];
   int  Rx     = appData[inpData_regionX];
   // int  Ry     = appData[inpData_regionY];
   int  Cx     = appData[inpData_cornerX];
   int  Cy     = appData[inpData_cornerY];

   int i,j,ii;
   // local composite ids are assumed to be lexicographical
   // on the owned rectangle. These need to be mapped to
   // the region rectangle.

   if (compLID < ownedX*ownedY) {
      i = (compLID+1)%ownedX;
      if (i==0) i=ownedX;

      j = (compLID+1 - i)/ownedX + 1;

      return(Rx*(j-1+Cy)+i+Cx -1);  // C-style ==> -1
   }
   else {
      ii = compLID - ownedX*ownedY;
      if (ii > appData[inpData_nGhosts] ) return(-1);
      if (LIDsOfGhosts[ii] == -1) return(-1);
      return(LIDsOfGhosts[ii]);
   }
}
