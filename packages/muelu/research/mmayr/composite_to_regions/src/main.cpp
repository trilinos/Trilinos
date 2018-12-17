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
#include <MueLu_CreateEpetraPreconditioner.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>

#include "Teuchos_Assert.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

int LIDregion(void *ptr, int LIDcomp, int whichGrp);
int LID2Dregion(void *ptr, int LIDcomp, int whichGrp);

extern void edgeGhosts(int ArowPtr[], int Acols[], int &nGhostFound, int ghostCompLIDs[], int edgeLength, int alongX, int ownedEdge, int interiorEdge, int start, int ownedX, int ownedY);

extern void fillCircleSquareData(int ArowPtr[], int Acols[], int ownedX, int ownedY, int Cx, int Cy, int Rx, int Ry, std::vector<int> &appData);

extern int LIDregionCircleSquare(void *ptr, int compLID,int whichGrp);

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

//! Print an Epetra_Vector in regional layout to screen
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

//! Print an Epetra_Map in regional layout to screen
void printRegionalMap(const std::string mapName, ///< string to be used for screen output
    const std::vector<Teuchos::RCP<Epetra_Map> > regMaps, ///< regional map to be printed to screen
    const int myRank ///< rank of calling proc
    )
{
//  sleep(myRank);
  for (int j = 0; j < (int) regMaps.size(); j++) {
    printf("%d: %s %d\n", myRank, mapName.c_str(), j);
    regMaps[j]->Print(std::cout);
  }
}

//! Print an Epetra_CrsMatrix in regional layout to screen
void printRegionalMatrix(const std::string matrixName, ///< string to be used for screen output
    const std::vector<Teuchos::RCP<Epetra_CrsMatrix> > regMats, ///< regional matrix to be printed to screen
    const int myRank ///< rank of calling proc
    )
{
  for (int j = 0; j < (int) regMats.size(); j++) {
    printf("%d: %s %d\n", myRank, matrixName.c_str(), j);
    regMats[j]->Print(std::cout);
  }
}

/*! \brief Transform composite vector to regional layout
 *
 *  Starting from a \c Epetra_Vector in composite layout, we
 *  1. import it into an auxiliary vector in the quasiRegional layout
 *  2. replace the quasiRegional map of the auxiliary vector with the regional map
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
 *
 *  Starting from a \c Epetra_Vector in regional layout, we
 *  1. replace the regional map with the quasiRegional map
 *  2. export it into a vector with composite layout using the \c CombineMode \c Add.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note The \c Add operation in the \c Export() can be done twofold, depending on the linear algebra package:
 *  - Epetra provides a \c CombineMode \c Eptra_AddLocalAlso that adds on-process values.
 *  - Tpetra does not provide such a capability (as it is not unique), so we perform the local summation
 *    manually.
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

/*! \brief Transform regional matrix to composite layout
 *
 *  Starting from a \c Epetra_CrsMatrix in regional layout, we
 *  1. copy data from regional matrix into quasiRegional matrix and set all maps
 *     to be quasiRegional maps
 *  2. export it into a Epetra_CrsMatrix with composite layout using the \c CombineMode \c Add.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note The \c Add operation in the \c Export() can be done twofold, depending on the linear algebra package:
 *  - Epetra provides a \c CombineMode \c Eptra_AddLocalAlso that adds on-process values.
 *  - Tpetra does not provide such a capability (as it is not unique), so we perform the local summation
 *    manually.
 */
void regionalToComposite(const std::vector<Teuchos::RCP<Epetra_CrsMatrix> >& regMat, ///< Matrix in region layout [in]
    Teuchos::RCP<Epetra_CrsMatrix> compMat, ///< Matrix in composite layout [in/out]
    const int maxRegPerProc, ///< max number of regions per proc
    const std::vector<Teuchos::RCP<Epetra_Map> > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
    const std::vector<Teuchos::RCP<Epetra_Map> > colMapPerGrp, ///< col maps in quasiRegion layout [in]
    const std::vector<Teuchos::RCP<Epetra_Import> > rowImportPerGrp, ///< row importer in region layout [in]
    const Epetra_CombineMode combineMode, ///< Combine mode for import/export [in]
    const bool addLocalManually ///< perform ADD of local values manually (not via Epetra CombineMode)
    )
{
  if (not addLocalManually) {
    /* Use the Eptra_AddLocalAlso combine mode to add processor-local values.
     * Note that such a combine mode is not available in Tpetra.
     */

    // Copy data from quasiRegionGrpMats, but into new map layout
    std::vector<Teuchos::RCP<Epetra_CrsMatrix> > quasiRegMat(maxRegPerProc);
    {
      for (int j = 0; j < maxRegPerProc; j++) {
        // create empty matrix with correct row and column map
        quasiRegMat[j] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(rowMapPerGrp[j]), *(colMapPerGrp[j]), 3));

        // extract pointers to crs arrays from region matrix
        Epetra_IntSerialDenseVector & rowPtr = regMat[j]->ExpertExtractIndexOffset();
        Epetra_IntSerialDenseVector & colInd = regMat[j]->ExpertExtractIndices();
        double *& vals = regMat[j]->ExpertExtractValues();

        // extract pointers to crs arrays from quasiRegion matrix
        Epetra_IntSerialDenseVector & qRowPtr = quasiRegMat[j]->ExpertExtractIndexOffset();
        Epetra_IntSerialDenseVector & qColInd = quasiRegMat[j]->ExpertExtractIndices();
        double *& qVals = quasiRegMat[j]->ExpertExtractValues();

        // assign array values from regional to quasiRegional matrices
        qRowPtr.Resize(rowPtr.Length());
        qColInd.Resize(colInd.Length());
        delete [] qVals;
        qVals = new double[qColInd.Length()];
        for (int i = 0; i < qRowPtr.Length(); ++i) qRowPtr[i] = rowPtr[i];
        for (int i = 0; i < qColInd.Length(); ++i) {
          qColInd[i] = colInd[i];
          qVals[i] = vals[i];
        }
      }

      /* add domain and range map to quasiRegion matrices (pass in quasiRow map, since
       * we assume that row map = range map = domain map)
       */
      for (int j = 0; j < maxRegPerProc; j++) {
        quasiRegMat[j]->ExpertStaticFillComplete(*(rowMapPerGrp[j]),*(rowMapPerGrp[j]));
      }

      for (int j = 0; j < maxRegPerProc; j++) {
        int err = compMat->Export(*quasiRegMat[j], *(rowImportPerGrp[j]), combineMode);
        TEUCHOS_ASSERT(err == 0);
      }

//      int err = compMat->FillComplete(compMat->RowMap(), compMat->RowMap());
      int err = compMat->FillComplete();
      TEUCHOS_ASSERT(err == 0);
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(false, "Not implemented, yet.");
//    /* Let's fake an ADD combine mode that also adds local values by
//     * 1. exporting quasiRegional vectors to auxiliary composite vectors (1 per group)
//     * 2. add all auxiliary vectors together
//     */
//
//    Teuchos::RCP<Epetra_MultiVector> partialCompVec = Teuchos::rcp(new Epetra_MultiVector(compVec->Map(), maxRegPerProc, true));
//    TEUCHOS_ASSERT(!partialCompVec.is_null());
//
//    std::vector<Teuchos::RCP<Epetra_Vector> > quasiRegVec(maxRegPerProc);
//    for (int j = 0; j < maxRegPerProc; j++) {
//      // copy vector and replace map
//      quasiRegVec[j] = Teuchos::rcp(new Epetra_Vector(*(regVec[j])));
//      TEUCHOS_ASSERT(!quasiRegVec[j].is_null());
//
//      int err = quasiRegVec[j]->ReplaceMap(*(rowMapPerGrp[j]));
//      TEUCHOS_ASSERT(err == 0);
//
//      err = (*partialCompVec)(j)->Export(*quasiRegVec[j], *(rowImportPerGrp[j]), Add);
//      TEUCHOS_ASSERT(err == 0);
//    }
//
//    compVec->PutScalar(0.0);
//    for (int j = 0; j < maxRegPerProc; j++) {
//      int err = compVec->Update(1.0, *(*partialCompVec)(j), 1.0);
//      TEUCHOS_ASSERT(err == 0);
//    }
  }

  return;
}

/*! \brief Sum region interface values
 *
 *  Sum values of interface GIDs using the underlying Export() routines. Technically, we perform the
 *  exchange/summation of interface data by exporting a regional vector to the composite layout and
 *  then immediately importing it back to the regional layout. The Export() involved when going to the
 *  composite layout takes care of the summation of interface values.
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

//! Create an empty Epetra_Vector in the regional layout
void createRegionalVector(std::vector<Teuchos::RCP<Epetra_Vector> >& regVecs, ///< regional vector to be filled
    const int maxRegPerProc, ///< max number of regions per process
    const std::vector<Teuchos::RCP<Epetra_Map> > revisedRowMapPerGrp ///< regional map
    )
{
  regVecs.resize(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++)
    regVecs[j] = Teuchos::rcp(new Epetra_Vector(*(revisedRowMapPerGrp[j]), true));
  return;
}

/*! \brief Compute the residual \f$r = b - Ax\f$
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute y = A*x in regional layout.
 *  2. Sum interface values of y to account for duplication of interface DOFs.
 *  3. Compute r = b - y
 */
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

/*! \brief Do Jacobi smoothing
 *
 *  Perform Jacobi smoothing in the region layout using the true diagonal value
 *  recovered from the splitted matrix.
 */
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

/*! \brief Create coarse level maps with continuous GIDs
 *
 *  The direct solver requires maps with continuous GIDs. Starting from the
 *  coarse level composite maps with discontinuous GIDs, we create a new row map
 *  and a matching column map.
 *
 *  Range and Domain map happen to correspond to the Row map, so we don't have
 *  to deal with them in particular.
 */
void createContinuousCoarseLevelMaps(const Epetra_Map& rowMap, ///< row map
    const Epetra_Map& colMap, ///< column map
    Teuchos::RCP<Epetra_Map>& contRowMap, ///< row map with continuous GIDs
    Teuchos::RCP<Epetra_Map>& contColMap ///< column map with continuous GIDs
    )
{
  // Create row map with continuous GIDs
  contRowMap = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),
      rowMap.NumMyElements(), 0, rowMap.Comm()));

  /* Create column map based on row map with continuous GIDs
   *
   * We use an Importer to create an auxiliary vector in colMap format containing
   * the GIDs of the contRowMap as its entries. By looping over its LIDs, we can
   * then form the contColMap.
   */
  Teuchos::RCP<Epetra_Import> rowColImport = Teuchos::rcp(new Epetra_Import(colMap, rowMap));
  Teuchos::RCP<Epetra_Vector> colGIDVec = Teuchos::rcp(new Epetra_Vector(rowMap, true));
  for (int i = 0; i < colGIDVec->MyLength(); ++i)
    (*colGIDVec)[i] = contRowMap->GID(i);
  Teuchos::RCP<Epetra_Vector> contColGIDVec = Teuchos::rcp(new Epetra_Vector(colMap, true));
  int err = contColGIDVec->Import(*colGIDVec, *rowColImport, Insert);
  TEUCHOS_ASSERT(err == 0);

  std::vector<int> contColGIDs;
  for (int i = 0; i < contColGIDVec->MyLength(); ++i) {
    contColGIDs.push_back((*contColGIDVec)[i]);
  }
  contColMap = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),
      contColGIDs.size(), contColGIDs.data(), 0, rowMap.Comm()));

  return;
}

//! Recursive V-cycle in region fashion
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
    Teuchos::Array<Teuchos::RCP<Epetra_Map> > compRowMaps, ///< composite maps
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Map> > > quasiRegRowMaps, ///< quasiRegional row maps
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Map> > > regRowMaps, ///< regional row maps
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Import> > > regRowImporters, ///< regional row importers
    Teuchos::Array<std::vector<Teuchos::RCP<Epetra_Vector> > > regInterfaceScalings, ///< regional interface scaling factors
    Teuchos::RCP<Epetra_CrsMatrix> coarseCompMat ///< Coarsest level composite operator
    )
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  if (l < numLevels - 1) { // fine or intermediate levels

    // pre-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
        regInterfaceScalings[l], maxRegPerProc, compRowMaps[l],
        quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    std::vector<Teuchos::RCP<Epetra_Vector> > regRes(maxRegPerProc);
    createRegionalVector(regRes, maxRegPerProc, regRowMaps[l]);
    computeResidual(regRes, fineRegX, fineRegB, regMatrices[l], compRowMaps[l],
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
    sumInterfaceValues(coarseRegB, compRowMaps[l+1], maxRegPerProc,
        quasiRegRowMaps[l+1], regRowMaps[l+1], regRowImporters[l+1]);

    // Call V-cycle recursively
    vCycle(l+1, numLevels, maxFineIter, maxCoarseIter, omega, maxRegPerProc,
        coarseRegX, coarseRegB, regMatrices, regProlong, compRowMaps,
        quasiRegRowMaps, regRowMaps, regRowImporters, regInterfaceScalings, coarseCompMat);

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
        regInterfaceScalings[l], maxRegPerProc, compRowMaps[l], quasiRegRowMaps[l],
        regRowMaps[l], regRowImporters[l]);
  }

  else {

    // create row and column maps with continuous GIDs
    Teuchos::RCP<Epetra_Map> contigRowMap = Teuchos::rcp(new Epetra_Map(coarseCompMat->RowMap()));
    Teuchos::RCP<Epetra_Map> contigColMap = Teuchos::rcp(new Epetra_Map(coarseCompMat->ColMap()));
    createContinuousCoarseLevelMaps(coarseCompMat->RowMap(),
        coarseCompMat->ColMap(), contigRowMap, contigColMap);

    // Store non contiguous maps for later
    RCP<Epetra_Map> noncontigRowMap = rcp(new Epetra_Map(coarseCompMat->RowMap()));
    RCP<Epetra_Map> noncontigColMap = rcp(new Epetra_Map(coarseCompMat->ColMap()));

    // Create composite error vector (zero initial guess)
    Teuchos::RCP<Epetra_Vector> compX = Teuchos::rcp(new Epetra_Vector(coarseCompMat->RowMap(), true));

    // Create composite right-hand side vector
    Teuchos::RCP<Epetra_Vector> compRhs = Teuchos::rcp(new Epetra_Vector(coarseCompMat->RowMap(), true));
    {
      for (int j = 0; j < maxRegPerProc; j++) {
        for (int i = 0; i < fineRegB[j]->MyLength(); ++i)
        (*fineRegB[j])[i] /= (*((regInterfaceScalings[l])[j]))[i];
      }

      regionalToComposite(fineRegB, compRhs, maxRegPerProc, quasiRegRowMaps[l],
        regRowImporters[l], Add, true);
    }

    // Replace non-continuos maps by continuous maps
    compX->ReplaceMap(*contigRowMap);
    compRhs->ReplaceMap(*contigRowMap);
    int err = coarseCompMat->ReplaceRowMap(*contigRowMap);
    TEUCHOS_ASSERT(err==0);
    err = coarseCompMat->ReplaceColMap(*contigColMap);
    TEUCHOS_ASSERT(err==0);
    err = coarseCompMat->ExpertStaticFillComplete(*contigRowMap, *contigRowMap);
    TEUCHOS_ASSERT(err==0);

    // create a linear problem object
    Epetra_LinearProblem problem(coarseCompMat.get(), &(*compX), &(*compRhs));

    // Direct solver
    {
      Teuchos::ParameterList pList;
      pList.set("PrintTiming",true);
      pList.set("PrintStatus",true);
      pList.set("MaxProcs", coarseCompMat->Comm().NumProc());

      Amesos Factory;
      Amesos_BaseSolver* solver = Factory.Create("Amesos_Umfpack", problem);
      TEUCHOS_ASSERT(solver!=NULL);

      solver->SetParameters(pList);
      solver->SetUseTranspose(false);

      solver->SymbolicFactorization();
      solver->NumericFactorization();
      solver->Solve();
    }

    // Replace maps with the original non-continuous maps
    compX->ReplaceMap(*noncontigColMap);
    compRhs->ReplaceMap(*noncontigRowMap);
  }

  return;
}

Teuchos::Array<int> setLocalNodesPerDim(const std::string& problemType,
    const bool doing1D, const Epetra_Map& rowMap, const int dimX, const int dimY, const int dimZ = 1)
{
  // Number of nodes per x/y/z-direction per processor
  Teuchos::Array<int> lNodesPerDim(3);

  if (problemType == "structured") {
    if (doing1D) { // One-dimensional problems
      lNodesPerDim[0] = rowMap.NumMyElements();
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
      if (rowMap.Comm().MyPID() == 4)
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
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  MPI_Group world_group;
  MPI_Comm_group(Comm.Comm(),&world_group);
#else
  Epetra_SerialComm Comm;
#endif
  int  myRank;
  FILE *fp;
  char command[40];
  bool doing1D = false;
  int  globalNx, globalNy;
  std::string xmlFileName;
  std::string problemType;
  std::string regionDataDirectory;

  myRank = Comm.MyPID();

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;

  Comm.Barrier();

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

  Comm.Barrier();

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
  Comm.Barrier();

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

  Comm.Barrier();

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
  Array<RCP<Epetra_Map> > compRowMaps; ///< composite row maps on each level
  Array<RCP<Epetra_Map> > compColMaps; ///< composite columns maps on each level
  Array<std::vector<RCP<Epetra_Map> > > regRowMaps; ///< regional row maps on each level
  Array<std::vector<RCP<Epetra_Map> > > quasiRegRowMaps; ///< quasiRegional row maps on each level
  Array<std::vector<RCP<Epetra_Map> > > regColMaps; ///< regional column maps on each level
  Array<std::vector<RCP<Epetra_Map> > > quasiRegColMaps; ///< quasiRegional column maps on each level
  Array<std::vector<RCP<Epetra_CrsMatrix> > > regMatrices; ///< regional matrices on each level
  Array<std::vector<RCP<Epetra_CrsMatrix> > > regProlong; ///< regional prolongators on each level
  Array<std::vector<RCP<Epetra_Import> > > regRowImporters; ///< regional row importers on each level
  Array<std::vector<RCP<Epetra_Vector> > > regInterfaceScalings; ///< regional interface scaling factors on each level

  Array<std::vector<std::vector<int> > > interfaceLIDs; // local IDs of interface nodes on each level in each group
  Array<std::vector<std::vector<std::vector<int> > > > interfaceGIDPairs; // pairs of GIDs of interface nodes on each level in each group

  Teuchos::RCP<Epetra_CrsMatrix> coarseCompOp = Teuchos::null;

  /* The actual computations start here. It's a sequence of operations to
   * - read region data from files
   * - create composite and region matrices
   * - create a region multigrid hierarchy
   * - run a region-type V-cycle
   */

  Comm.Barrier();

  // Load composite map
  {
    std::cout << myRank << " | Loading composite map ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myCompositeMap_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    if (fp == NULL)
      std::cout << std::endl << ">>> Check number of MPI ranks!" << std::endl << std::endl;
    TEUCHOS_ASSERT(fp!=NULL);

    std::vector<int> fileData; // composite GIDs
    int i;
    while (fscanf(fp, "%d", &i) != EOF)
      fileData.push_back(i);

    mapComp = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), (int) fileData.size(), fileData.data(), 0, Comm));
//    mapComp->Print(std::cout);
  }

  Comm.Barrier();

  // Read matrix from file
  {
    std::cout << myRank << " | Reading matrix from file ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/Amat.mm";

    Epetra_CrsMatrix* ACompPtr;
    EpetraExt::MatrixMarketFileToCrsMatrix(fileNameSS.str().c_str(), *mapComp, ACompPtr);
    AComp = Teuchos::rcp(new Epetra_CrsMatrix(*ACompPtr));
//    AComp->Print(std::cout);
  }

  Comm.Barrier();

  // Load and communicate region assignments
  {
    std::cout << myRank << " | Loading and communicating region assignments ..." << std::endl;

    regionsPerGID = Teuchos::rcp(new Epetra_MultiVector(AComp->RowMap(),maxRegPerGID));

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegionAssignment_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    int k;
    double *jthRegions; // pointer to jth column in regionsPerGID
    for (int i = 0; i < mapComp->NumMyElements(); i++) {
      for (int j = 0; j < maxRegPerGID; j++) {
        jthRegions = (*regionsPerGID)[j];

        if (fscanf(fp, "%d", &k) == EOF) {
          fprintf(stderr, "not enough region assignments\n");
          exit(1);
        }
        jthRegions[i] = (double) k;

        // identify interface DOFs. An interface DOF is assigned to at least two regions.
        if (j > 0 and k != -1)
          intIDs.push_back(i);
      }
    }

    // make extended Region Assignments
    Epetra_Import Importer(AComp->ColMap(), *mapComp);
    regionsPerGIDWithGhosts = Teuchos::rcp(new Epetra_MultiVector(AComp->ColMap(),maxRegPerGID));
    regionsPerGIDWithGhosts->Import(*regionsPerGID, Importer, Insert);
//    regionsPerGIDWithGhosts->Print(std::cout);
  }

  Comm.Barrier();

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

  Comm.Barrier();

  std::vector<int> genericVector;
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
      int isStructured;
      retval = fscanf(fp,"%d",&isStructured);
      if (isStructured == 0) {
         genericVector.resize(3);
         genericVector[inpData_isStructured] = 0;
         retval = fscanf(fp,"%d",&(genericVector[inpData_ownedX]));
         genericVector[inpData_ownedY] = 1;
      }
      else {
         int Cx, Cy, Rx, Ry;
         int ownedX,ownedY;
         retval = fscanf(fp,"%d%d",&Rx,&Ry);
         retval = fscanf(fp,"%d%d",&ownedX,&ownedY);
         retval = fscanf(fp,"%d%d",&Cx,&Cy);
         int *rowptr, *cols;  double *values;
         AComp->ExtractCrsDataPointers(rowptr, cols, values);
         fillCircleSquareData(rowptr, cols, ownedX,ownedY, Cx, Cy, Rx, Ry,genericVector);
      }
      fclose(fp);
    }

    appData.maxRegPerGID = maxRegPerGID;
    appData.myRegions = myRegions.data();
    appData.colMap = (Epetra_Map *) &(AComp->ColMap());
    appData.regionsPerGIDWithGhosts = regionsPerGIDWithGhosts;
    appData.myRank = myRank;
  }

  Comm.Barrier();

  // Make group region row maps
  {
    std::cout << myRank << " | Creating region group row maps ..." << std::endl;

    std::vector<int> rowGIDsReg;
    int *colGIDsComp = AComp->ColMap().MyGlobalElements();
//sleep(myRank*3);
    for (int k = 0; k < (int) myRegions.size(); k++) {
      rowGIDsReg.resize(0);
      std::vector<int> tempRegIDs(AComp->ColMap().NumMyElements());
      for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {

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

      for (int i = 0; i < AComp->ColMap().NumMyElements(); i++) {
        if (tempRegIDs[idx[i]] != -1)
          rowGIDsReg.push_back(colGIDsComp[idx[i]]);
      }
      rowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), (int) rowGIDsReg.size(),
          rowGIDsReg.data(), 0, Comm));
    }

    for (int k=(int) myRegions.size(); k < maxRegPerProc; k++) {
      rowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),0,NULL,0,Comm));
    }
  }

  Comm.Barrier();

  // Make group region column maps
  {
    std::cout << myRank << " | Creating region group column maps ..." << std::endl;

    if (whichCase == MultipleRegionsPerProc) {
      // clone rowMap
      for (int j=0; j < maxRegPerProc; j++) {
        if (j < (int) myRegions.size()) {
          colMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),
              rowMapPerGrp[j]->NumMyElements(),
              rowMapPerGrp[j]->MyGlobalElements(), 0, Comm));
        }
        else colMapPerGrp[j] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),0,NULL,0,Comm));
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
            jthRegions = (*regionsPerGIDWithGhosts)[j];
            if  ( ((int) jthRegions[i]) == myRegions[0]) {
              colIDsReg.push_back(colGIDsComp[i]);
              break;
            }
          }
        }
      }
      if ((int) myRegions.size() > 0) {
       colMapPerGrp[0] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), (int) colIDsReg.size(),
           colIDsReg.data(), 0, Comm));
      }
      else colMapPerGrp[0] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),0,NULL,0,Comm));
    }
    else { fprintf(stderr,"whichCase not set properly\n"); exit(1); }
  }

  Comm.Barrier();

  // Make extended group region maps
  {
    std::cout << myRank << " | Creating extended group region maps ..." << std::endl;

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

      revisedRowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),(int) revisedGIDs.size(),
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
      revisedColMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), (int) revisedGIDs.size(),
          revisedGIDs.data(), 0, Comm));
    }
    for (int k = (int) myRegions.size(); k < maxRegPerProc; k++) {
      revisedRowMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),0,NULL,0,Comm));
      revisedColMapPerGrp[k] = Teuchos::rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(),0,NULL,0,Comm));
    }

    // Setup importers
    for (int j = 0; j < maxRegPerProc; j++) {
      rowImportPerGrp[j] = Teuchos::rcp(new Epetra_Import(*(rowMapPerGrp[j]), *mapComp));
      colImportPerGrp[j] = Teuchos::rcp(new Epetra_Import(*(colMapPerGrp[j]), *mapComp));
    }
  }

  Comm.Barrier();

  // Make quasiRegion matrices
  {
    std::cout << myRank << " | Forming quasiRegion matrices ..." << std::endl;

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

  Comm.Barrier();

  // Make region matrices
  {
    std::cout << myRank << " | Forming region matrices ..." << std::endl;

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
        vals = new double[colInd.Length()];
        for (int i = 0; i < rowPtr.Length(); ++i) rowPtr[i] = qRowPtr[i];
        for (int i = 0; i < colInd.Length(); ++i) {
          colInd[i] = qColInd[i];
          vals[i] = qVals[i];
        }
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

  Comm.Barrier();

//  /* Form composite operator on fine level for debug purposes */
//  {
//    Teuchos::RCP<Epetra_CrsMatrix> compOp = Teuchos::rcp(new Epetra_CrsMatrix(Copy, AComp->RowMap(), AComp->ColMap(), 3));
//    regionalToComposite(regionGrpMats, compOp, maxRegPerProc, rowMapPerGrp, colMapPerGrp, rowImportPerGrp, Add, false);
//
//    sleep(1);
//    std::cout << myRank << " | Printing compOp ..." << std::endl;
//    Comm.Barrier();
//    compOp->Print(std::cout);
//
//    Teuchos::RCP<Epetra_CrsMatrix> diffOp = Teuchos::rcp(new Epetra_CrsMatrix(Copy, AComp->Graph()));
//    Epetra_CrsMatrix* diffOp_ptr = diffOp.get();
//    EpetraExt::MatrixMatrix::Add(*AComp, false, 1.0, *compOp, false, -1.0, diffOp_ptr);
//
//    sleep(1);
//    std::cout << myRank << " | Printing diffOp ..." << std::endl;
//    Comm.Barrier();
//    diffOp->Print(std::cout);
//
//    Comm.Barrier();
//    std::cout << myRank << " | Calling exit(0) ..." << std::endl;
//    Comm.Barrier();
//    exit(0);
//  }

  Comm.Barrier();

  // Create multigrid hierarchy
  {
    std::cout << myRank << " | Setting up MueLu hierarchies ..." << std::endl;

    typedef MueLu::Hierarchy<double,int,int,Xpetra::EpetraNode> Hierarchy;
    typedef MueLu::Utilities<double,int,int,Xpetra::EpetraNode> Utilities;
    typedef Xpetra::Map<int,int,Xpetra::EpetraNode> Map;
    typedef Xpetra::MultiVector<double,int,int,Xpetra::EpetraNode> MultiVector;
    typedef Xpetra::MultiVectorFactory<double,int,int,Xpetra::EpetraNode> MultiVectorFactory;
    typedef Xpetra::Matrix<double,int,int,Xpetra::EpetraNode> Matrix;

    using Teuchos::RCP; using Teuchos::rcp; using Teuchos::ParameterList;

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
      RCP<Epetra_Vector> nullspace = rcp(new Epetra_Vector(*revisedRowMapPerGrp[j]));
      nullspace->PutScalar(1.0);

      // create dummy coordinates vector
      RCP<Epetra_MultiVector> coordinates = rcp(new Epetra_MultiVector(*revisedRowMapPerGrp[j], 3));
      coordinates->PutScalar(1.0);

      // Read MueLu parameter list form xml file
      RCP<ParameterList> mueluParams = Teuchos::getParametersFromXmlFile(xmlFileName);

      // Insert region-specific data into parameter list
      const std::string userName = "user data";
      Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);
      userParamList.set<Array<int> >("Array<LO> lNodesPerDim", lNodesPerDim);
      userParamList.set<std::string>("string aggregationRegionType", aggregationRegionType);

      // Setup hierarchy
      RCP<MueLu::EpetraOperator> eH = MueLu::CreateEpetraPreconditioner(regionGrpMats[j], *mueluParams,
          coordinates, nullspace);
      regGrpHierarchy[j] = eH->GetHierarchy();
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
      std::vector<RCP<Epetra_CrsMatrix> > fineLevelProlong(maxRegPerProc);
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

        RCP<Matrix> regPXpetra = level->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
        regProlong[l][j] = Utilities::Op2NonConstEpetraCrs(regPXpetra);

        RCP<Matrix> regRAPXpetra = level->Get<RCP<Matrix> >("A", MueLu::NoFactory::get());
        regMatrices[l][j] = Utilities::Op2NonConstEpetraCrs(regRAPXpetra);

        regRowMaps[l][j] = Teuchos::rcp(new Epetra_Map(regMatrices[l][j]->RowMap())); // ToDo (mayr.mt) Do not copy!
        regColMaps[l][j] = Teuchos::rcp(new Epetra_Map(regMatrices[l][j]->ColMap())); // ToDo (mayr.mt) Do not copy!
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

      RCP<Epetra_CrsMatrix> summedDuplicateMatrix = Teuchos::null;

      // Find lists of duplicated GIDs locally
      interfaceLIDs.resize(numLevels);
      interfaceGIDPairs.resize(numLevels);
      for (int l = 0; l < numLevels; ++l) {
        interfaceLIDs[l].resize(maxRegPerProc);
        interfaceGIDPairs[l].resize(maxRegPerProc);
      }

      for (int l = 0; l < numLevels - 1; ++l) {
//        Comm.Barrier();
//        if (Comm.MyPID() == 0) {
//          std::cout << std::endl << std::endl
//              << "Processing GID pairs on level " << l
//              << std::endl << std::endl;
//        }
//        Comm.Barrier();

//        sleep(1);
//        std::cout << "Prolongator" << std::endl;
//        regProlong[l+1][0]->Print(std::cout);

        // create list of LIDs per group
        for (int j = 0; j < maxRegPerProc; j++) {
          for (int i = 0; i < regRowMaps[l][j]->NumMyElements(); ++i) {
            if (regRowMaps[l][j]->GID(i) != quasiRegRowMaps[l][j]->GID(i)) {
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
            interfaceGIDPairs[l][j][i].push_back(quasiRegRowMaps[l][j]->GID(interfaceLIDs[l][j][i]));
        for (int j = 0; j < maxRegPerProc; j++)
          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
            if (regRowMaps[l][j]->GID(interfaceLIDs[l][j][i]) != quasiRegRowMaps[l][j]->GID(interfaceLIDs[l][j][i]))
              interfaceGIDPairs[l][j][i].push_back(regRowMaps[l][j]->GID(interfaceLIDs[l][j][i]));

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

        std::vector<RCP<Epetra_Vector> > regDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
        std::vector<RCP<Epetra_Vector> > quasiRegDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
        for (int j = 0; j < maxRegPerProc; ++j) {
          regDupGIDVec[j] = rcp(new Epetra_Vector(*regRowMaps[l][j], true));
          quasiRegDupGIDVec[j] = rcp(new Epetra_Vector(*quasiRegRowMaps[l][j], true));

          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
            (*regDupGIDVec[j])[interfaceLIDs[l][j][i]] = 1.0;
        }

        RCP<Epetra_Vector> compDupGIDVec = rcp(new Epetra_Vector(*compRowMaps[l]), true);
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

//          sleep(1);
//          Comm.Barrier();
//          std::cout << myRank << " | duplicateMap:" << std::endl;
//          duplicateMap->Print(std::cout);
        }

        // create row/range/domain map for the transpose of the fine level duplicates mapping matrix
        RCP<Epetra_Map> fullDuplicateMap = Teuchos::null;
        {
          std::vector<int> myIntGIDs;
          for (int i = 0; i < regRowMaps[l][0]->NumMyElements(); ++i) {
            if ((*regDupGIDVec[0])[i] != 0)
              myIntGIDs.push_back(regRowMaps[l][0]->GID(i));
          }

//          std::cout << myRank << " | myIntGIDs = ";
//          for (auto gid : myIntGIDs)
//            std::cout << gid << ", ";
//          std::cout << std::endl;

          fullDuplicateMap = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myIntGIDs.size(), myIntGIDs.data(), 0, Comm));

//          sleep(1);
//          Comm.Barrier();
//          std::cout << myRank << " | fullDuplicateMap:" << std::endl;
//          fullDuplicateMap->Print(std::cout);
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

//        sleep(1);
//        Comm.Barrier();
//        std::cout << myRank << " | duplicateMatrix:" << std::endl;
//        duplicateMatrix->Print(std::cout);
//
//        sleep(1);
//        Comm.Barrier();
//        std::cout << myRank << " | duplicateMatrix->RangeMap():" << std::endl;
//        duplicateMatrix->OperatorRangeMap().Print(std::cout);
//
//        sleep(1);
//        Comm.Barrier();
//        std::cout << myRank << " | duplicateMatrix->DomainMap():" << std::endl;
//        duplicateMatrix->OperatorDomainMap().Print(std::cout);
//
//        sleep(1);
//        Comm.Barrier();
//        std::cout << myRank << " | duplicateMatrix->ColMap():" << std::endl;
//        duplicateMatrix->ColMap().Print(std::cout);

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
//        err = summedDuplicateMatrix->FillComplete();
        TEUCHOS_ASSERT(err == 0);

//        sleep(1);
//        Comm.Barrier();
//        std::cout << myRank << " | summedDuplicateMatrix:" << std::endl;
//        summedDuplicateMatrix->Print(std::cout);

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
//          Comm.Barrier();
//          regRowMaps[l+1][0]->Print(std::cout);
//          sleep(2);

          // create quasiRegional row map
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

//            sleep(1);
//            std::cout << myRank << " | Printing quasiRegRowMaps[" << l+1 << "][0] ..." << std::endl;
//            Comm.Barrier();
//            quasiRegRowMaps[l+1][0]->Print(std::cout);
          }

          // create composite row map
          {
            std::vector<int> myCompGIDs;
            for (int j = 0; j < maxRegPerProc; j++) {
              for (int i = 0; i < quasiRegRowMaps[l+1][j]->NumMyElements(); ++i) {
                const int trialGID = quasiRegRowMaps[l+1][j]->GID(i);

                if (regRowMaps[l+1][j]->MyGID(trialGID))
                  myCompGIDs.push_back(trialGID);
              }
            }

            compRowMaps[l+1] = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myCompGIDs.size(), myCompGIDs.data(), 0, Comm));
            TEUCHOS_ASSERT(!compRowMaps[l+1].is_null());

//            sleep(1);
//            std::cout << myRank << " | Printing compRowMaps["<< l+1 << "] ..." << std::endl;
//            Comm.Barrier();
//            compRowMaps[l+1]->Print(std::cout);
          }

          // create regRowImporter
          for (int j = 0; j < maxRegPerProc; ++j) {
            regRowImporters[l+1][j] =
                rcp(new Epetra_Import(*quasiRegRowMaps[l+1][j], *compRowMaps[l+1]));
            TEUCHOS_ASSERT(!regRowImporters[l+1][j].is_null());
          }

          // Create quasiRegional column map
          {
            std::vector<int> myQuasiRegGIDs;

            for (int i = 0; i < regColMaps[l+1][0]->NumMyElements(); ++i) {
              // grab current regional GID to be processed
              int currGID = regColMaps[l+1][0]->GID(i);
              int quasiGID = currGID; // assign dummy value

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

            quasiRegColMaps[l+1][0] = rcp(new Epetra_Map(Teuchos::OrdinalTraits<int>::invalid(), myQuasiRegGIDs.size(), myQuasiRegGIDs.data(), 0, Comm));
            TEUCHOS_ASSERT(!quasiRegColMaps[l+1][0].is_null());

//            sleep(1);
//            std::cout << myRank << " | Printing quasiRegColMaps[" << l+1 << "][0] ..." << std::endl;
//            Comm.Barrier();
//            quasiRegColMaps[l+1][0]->Print(std::cout);
          }
        }
      }
    }

    // Form the composite coarse level operator
    {
      const int maxLevel = numLevels - 1;
      coarseCompOp = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *compRowMaps[maxLevel], 3));
      regionalToComposite(regMatrices[maxLevel], coarseCompOp, maxRegPerProc,
          quasiRegRowMaps[maxLevel], quasiRegColMaps[maxLevel],
          regRowImporters[maxLevel], Add, false);

      int err = coarseCompOp->FillComplete(*compRowMaps[maxLevel], *compRowMaps[maxLevel]);
      TEUCHOS_ASSERT(err==0);


//      sleep(1);
//      std::cout << myRank << " | Printing coarseCompOp ..." << std::endl;
//      Comm.Barrier();
//      coarseCompOp->Print(std::cout);
    }
  }

  Comm.Barrier();

  // Make interface scaling factors recursively
  {
    std::cout << myRank << " | Computing interface scaling factors ..." << std::endl;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probaly, numLevel has not been set, yet.");

    for (int l = 0; l < numLevels; l++)
    {
      // initialize region vector with all ones.
      for (int j = 0; j < maxRegPerProc; j++) {
        regInterfaceScalings[l][j] = rcp(new Epetra_Vector(*regRowMaps[l][j], true));
        regInterfaceScalings[l][j]->PutScalar(1.0);
      }

      // transform to composite layout while adding interface values via the Export() combine mode
      RCP<Epetra_Vector> compInterfaceScalingSum = rcp(new Epetra_Vector(*compRowMaps[l], true));
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

  Comm.Barrier();

  // Run V-cycle
  {
    std::cout << myRank << " | Running V-cycle ..." << std::endl;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probaly, numLevel has not been set, yet.");

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

    if (doing1D)
    {
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
    }
    else //2D
    {
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

    // Prepare output of residual norm to file
    RCP<std::ofstream> log;
    if (myRank == 0)
    {
      std::string s = "residual_norm.txt";
      log = rcp(new std::ofstream(s.c_str()));
      (*log) << "# num procs = " << AComp->Comm().NumProc() << "\n"
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

        compRes = Teuchos::rcp(new Epetra_Vector(*mapComp, true));
        regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
            rowImportPerGrp, Add, true);
        double normRes = 0.0;
        compRes->Norm2(&normRes);

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
*/
  }

  Comm.Barrier();

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

/* Returns local ID (within region curRegion) for the LIDcomp^th composite grid
 * point that myRank owns. If this grid point is not part of curRegion, then
 * -1 is returned.
 *
 * For fully structured 2D problems only.
 */
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

int LIDregionCircleSquare(void *ptr, int compLID, int whichGrp)
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
   int  Ry     = appData[inpData_regionY];
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
