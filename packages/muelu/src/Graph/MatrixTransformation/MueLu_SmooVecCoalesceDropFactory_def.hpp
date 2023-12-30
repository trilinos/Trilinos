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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_SmootherBase.hpp"
#ifndef MUELU_SMOOVECCOALESCEDROPFACTORY_DEF_HPP
#define MUELU_SMOOVECCOALESCEDROPFACTORY_DEF_HPP

#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "MueLu_SmooVecCoalesceDropFactory_decl.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"

#include <Xpetra_IO.hpp>

#include <algorithm>
#include <cstdlib>
#include <string>

// If defined, read environment variables.
// Should be removed once we are confident that this works.
// #define DJS_READ_ENV_VARIABLES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define poly0thOrderCoef 0
#define poly1stOrderCoef 1
#define poly2ndOrderCoef 2
#define poly3rdOrderCoef 3
#define poly4thOrderCoef 4

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: drop scheme");
  {
    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
    validParamList->getEntry("aggregation: drop scheme").setValidator(rcp(new validatorType(Teuchos::tuple<std::string>("unsupported vector smoothing"), "aggregation: drop scheme")));
  }
  SET_VALID_ENTRY("aggregation: number of random vectors");
  SET_VALID_ENTRY("aggregation: number of times to pre or post smooth");
  SET_VALID_ENTRY("aggregation: penalty parameters");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("PreSmoother", Teuchos::null, "Generating factory of the PreSmoother");
  validParamList->set<RCP<const FactoryBase> >("PostSmoother", Teuchos::null, "Generating factory of the PostSmoother");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SmooVecCoalesceDropFactory()
  : predrop_(Teuchos::null) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  if (currentLevel.IsAvailable("PreSmoother")) {        // rst: totally unsure that this is legal
    Input(currentLevel, "PreSmoother");                 //      my guess is that this is not yet available
  }                                                     //      so this always comes out false.
  else if (currentLevel.IsAvailable("PostSmoother")) {  // perhaps we can look on the param list?
    Input(currentLevel, "PostSmoother");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  typedef Teuchos::ScalarTraits<SC> STS;

  if (predrop_ != Teuchos::null)
    GetOStream(Parameters0) << predrop_->description();

  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");

  const ParameterList& pL = GetParameterList();

  LO nPDEs = A->GetFixedBlockSize();

  RCP<MultiVector> testVecs;
  RCP<MultiVector> nearNull;

#ifdef takeOut
  testVecs = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector("TpetraTVecs.mm", A->getRowMap());
#endif
  size_t numRandom = as<size_t>(pL.get<int>("aggregation: number of random vectors"));
  testVecs         = MultiVectorFactory::Build(A->getRowMap(), numRandom, true);
  // use random test vectors but should be positive in order to not get
  // crummy results ... so take abs() of randomize().
  testVecs->randomize();
  for (size_t kk = 0; kk < testVecs->getNumVectors(); kk++) {
    Teuchos::ArrayRCP<Scalar> curVec = testVecs->getDataNonConst(kk);
    for (size_t ii = kk; ii < as<size_t>(A->getRowMap()->getLocalNumElements()); ii++) curVec[ii] = Teuchos::ScalarTraits<SC>::magnitude(curVec[ii]);
  }
  nearNull = MultiVectorFactory::Build(A->getRowMap(), nPDEs, true);

  // initialize null space to constants
  for (size_t kk = 0; kk < nearNull->getNumVectors(); kk++) {
    Teuchos::ArrayRCP<Scalar> curVec = nearNull->getDataNonConst(kk);
    for (size_t ii = kk; ii < as<size_t>(A->getRowMap()->getLocalNumElements()); ii += nearNull->getNumVectors()) curVec[ii] = Teuchos::ScalarTraits<Scalar>::one();
  }

  RCP<MultiVector> zeroVec_TVecs;
  RCP<MultiVector> zeroVec_Null;

  zeroVec_TVecs = MultiVectorFactory::Build(A->getRowMap(), testVecs->getNumVectors(), true);
  zeroVec_Null  = MultiVectorFactory::Build(A->getRowMap(), nPDEs, true);
  zeroVec_TVecs->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
  zeroVec_Null->putScalar(Teuchos::ScalarTraits<Scalar>::zero());

  size_t nInvokeSmoother = as<size_t>(pL.get<int>("aggregation: number of times to pre or post smooth"));
  if (currentLevel.IsAvailable("PreSmoother")) {
    RCP<SmootherBase> preSmoo = currentLevel.Get<RCP<SmootherBase> >("PreSmoother");
    for (size_t ii = 0; ii < nInvokeSmoother; ii++) preSmoo->Apply(*testVecs, *zeroVec_TVecs, false);
    for (size_t ii = 0; ii < nInvokeSmoother; ii++) preSmoo->Apply(*nearNull, *zeroVec_Null, false);
  } else if (currentLevel.IsAvailable("PostSmoother")) {
    RCP<SmootherBase> postSmoo = currentLevel.Get<RCP<SmootherBase> >("PostSmoother");
    for (size_t ii = 0; ii < nInvokeSmoother; ii++) postSmoo->Apply(*testVecs, *zeroVec_TVecs, false);
    for (size_t ii = 0; ii < nInvokeSmoother; ii++) postSmoo->Apply(*nearNull, *zeroVec_Null, false);
  } else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Must set a smoother");

  Teuchos::ArrayRCP<Scalar> penaltyPolyCoef(5);
  Teuchos::ArrayView<const double> inputPolyCoef;

  penaltyPolyCoef[poly0thOrderCoef] = 12.;
  penaltyPolyCoef[poly1stOrderCoef] = -.2;
  penaltyPolyCoef[poly2ndOrderCoef] = 0.0;
  penaltyPolyCoef[poly3rdOrderCoef] = 0.0;
  penaltyPolyCoef[poly4thOrderCoef] = 0.0;

  if (pL.isParameter("aggregation: penalty parameters") && pL.get<Teuchos::Array<double> >("aggregation: penalty parameters").size() > 0) {
    if (pL.get<Teuchos::Array<double> >("aggregation: penalty parameters").size() > penaltyPolyCoef.size())
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Number of penalty parameters must be " << penaltyPolyCoef.size() << " or less");
    inputPolyCoef = pL.get<Teuchos::Array<double> >("aggregation: penalty parameters")();

    for (size_t i = 0; i < as<size_t>(inputPolyCoef.size()); i++) penaltyPolyCoef[i] = as<Scalar>(inputPolyCoef[i]);
    for (size_t i = as<size_t>(inputPolyCoef.size()); i < as<size_t>(penaltyPolyCoef.size()); i++) penaltyPolyCoef[i] = Teuchos::ScalarTraits<Scalar>::zero();
  }

  RCP<GraphBase> filteredGraph;
  badGuysCoalesceDrop(*A, penaltyPolyCoef, nPDEs, *testVecs, *nearNull, filteredGraph);

#ifdef takeOut
  /* write out graph for serial debugging purposes only. */

  FILE* fp = fopen("codeOutput", "w");
  fprintf(fp, "%d %d %d\n", (int)filteredGraph->GetNodeNumVertices(), (int)filteredGraph->GetNodeNumVertices(),
          (int)filteredGraph->GetNodeNumEdges());
  for (size_t i = 0; i < filteredGraph->GetNodeNumVertices(); i++) {
    ArrayView<const LO> inds = filteredGraph->getNeighborVertices(as<LO>(i));
    for (size_t j = 0; j < as<size_t>(inds.size()); j++) {
      fprintf(fp, "%d %d 1.00e+00\n", (int)i + 1, (int)inds[j] + 1);
    }
  }
  fclose(fp);
#endif

  SC threshold = .01;
  Set<bool>(currentLevel, "Filtering", (threshold != STS::zero()));
  Set(currentLevel, "Graph", filteredGraph);
  Set(currentLevel, "DofsPerNode", 1);

}  // Build

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::badGuysCoalesceDrop(const Matrix& Amat, Teuchos::ArrayRCP<Scalar>& penaltyPolyCoef, LO nPDEs, const MultiVector& testVecs, const MultiVector& nearNull, RCP<GraphBase>& filteredGraph) const {
  /*
   * Compute coalesce/drop graph (in filteredGraph) for A. The basic idea is to
   * balance trade-offs associated with
   *
   *           (I - P inv(R P) R ) testVecs
   *
   * being worse for larger aggregates (less dropping) while MG cycle costs are
   * cheaper with larger aggregates. MG costs are "penalties" in the
   * optimization while (I - P inv(R P) R ) is the "fit" (how well a
   * a fine grid function can be approximated on the coarse grid).
   *
   * For MG costs, we don't actually approximate the cost. Instead, we
   * have just hardwired penalties below. Specifically,
   *
   *     penalties[j] is the cost if aggregates are of size j+1, where right
   *                  now a linear function of the form const*(60-j) is used.
   *
   * (I - P inv(P^T P) P^T ) testVecs is estimated by just looking locally at
   * the vector portion corresponding to a possible aggregate defined by
   * all non-dropped connections in the ith row.  A tentative prolognator is
   * used for P. This prolongator corresponds to a null space vector given
   * by 'nearNull', which is provided to dropper(). In initial testing, nearNull is
   * first set as a vector of all 1's and then smoothed with a relaxation
   * method applied to a nice matrix (with the same sparsity pattern as A).
   * Originally, nearNull was used to handle Dir bcs where relaxation of the
   * vector of 1's has a more pronounced effect.
   *
   * For PDE systems, fit only considers the same dof at each node. That is,
   * it effectively assumes that we have a tentative prolongator with no
   * coupling between different dof types. When checking the fit for the kth
   * dof at a paritcular node, it only considers the kth dof of this node
   * and neighboring nodes.
   *
   * Note: testVecs is supplied by the user, but normally is the result of
   * applying a relaxation scheme to Au = 0 where u is initial random.
   */

  GO numMyNnz = Teuchos::as<GO>(Amat.getLocalNumEntries());
  size_t nLoc = Amat.getRowMap()->getLocalNumElements();

  size_t nBlks = nLoc / nPDEs;
  if (nBlks * nPDEs != nLoc)
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Number of local dofs not divisible by BlkSize");

  Teuchos::ArrayRCP<LO> newRowPtr(nBlks + 1); /* coalesce & drop matrix     */
  Teuchos::ArrayRCP<LO> newCols(numMyNnz);    /* arrays                     */

  Teuchos::ArrayRCP<LO> bcols(nBlks);       /* returned by dropfun(j,...) */
  Teuchos::ArrayRCP<bool> keepOrNot(nBlks); /* gives cols for jth row and */
                                            /* whether or not entry is    */
                                            /* kept or dropped.           */

  LO maxNzPerRow = 200;
  Teuchos::ArrayRCP<Scalar> penalties(maxNzPerRow); /* Penalty function */
                                                    /* described above. */

  Teuchos::ArrayRCP<bool> keepStatus(nBlks, true); /* accumulated keepOrNot info */
  Teuchos::ArrayRCP<LO> bColList(nBlks);           /* accumulated bcols info     */
                                                   /* for an entire block as     */
                                                   /* opposed to a single row    */
                                                   /* Additionally, keepOrNot[j] */
                                                   /* refers to status of jth    */
                                                   /* entry in a row while       */
                                                   /* keepStatus[j] refers to    */
                                                   /* whether the jth block is   */
                                                   /* kept within the block row. */

  Teuchos::ArrayRCP<bool> alreadyOnBColList(nBlks, false); /* used to avoid recording the*/
                                                           /* same block column when     */
                                                           /* processing different pt    */
                                                           /* rows within a block.       */

  Teuchos::ArrayRCP<bool> boundaryNodes(nBlks, false);

  for (LO i = 0; i < maxNzPerRow; i++)
    penalties[i] = penaltyPolyCoef[poly0thOrderCoef] +
                   penaltyPolyCoef[poly1stOrderCoef] * (as<Scalar>(i)) +
                   penaltyPolyCoef[poly2ndOrderCoef] * (as<Scalar>(i * i)) +
                   (penaltyPolyCoef[poly3rdOrderCoef] * (as<Scalar>(i * i)) * (as<Scalar>(i))) +  // perhaps avoids overflow?
                   (penaltyPolyCoef[poly4thOrderCoef] * (as<Scalar>(i * i)) * (as<Scalar>(i * i)));

  LO nzTotal = 0, numBCols = 0, row = -1, Nbcols, bcol;
  newRowPtr[0] = 0;

  /* proceed block by block */
  for (LO i = 0; i < as<LO>(nBlks); i++) {
    newRowPtr[i + 1] = newRowPtr[i];
    for (LO j = 0; j < nPDEs; j++) {
      row = row + 1;

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;

      Amat.getLocalRowView(row, indices, vals);

      if (indices.size() > maxNzPerRow) {
        LO oldSize  = maxNzPerRow;
        maxNzPerRow = indices.size() + 100;
        penalties.resize(as<size_t>(maxNzPerRow), 0.0);
        for (LO k = oldSize; k < maxNzPerRow; k++)
          penalties[k] = penaltyPolyCoef[poly0thOrderCoef] +
                         penaltyPolyCoef[poly1stOrderCoef] * (as<Scalar>(i)) +
                         penaltyPolyCoef[poly2ndOrderCoef] * (as<Scalar>(i * i)) +
                         (penaltyPolyCoef[poly3rdOrderCoef] * (as<Scalar>(i * i)) * (as<Scalar>(i))) +
                         (penaltyPolyCoef[poly4thOrderCoef] * (as<Scalar>(i * i)) * (as<Scalar>(i * i)));
      }
      badGuysDropfunc(row, indices, vals, testVecs, nPDEs, penalties, nearNull, bcols, keepOrNot, Nbcols, nLoc);
      for (LO k = 0; k < Nbcols; k++) {
        bcol = bcols[k];

        /* add to bColList if not already on it */

        if (alreadyOnBColList[bcol] == false) { /* for PDE systems only record */
          bColList[numBCols++]    = bcol;       /* neighboring block one time  */
          alreadyOnBColList[bcol] = true;
        }
        /* drop if any pt row within block indicates entry should be dropped */

        if (keepOrNot[k] == false) keepStatus[bcol] = false;

      } /* for (k=0; k < Nbcols; k++)   */
    }   /* for (j = 0; i < nPDEs; j++)  */

    /* finished with block row. Now record block entries that we keep */
    /* and reset keepStatus, bColList, and alreadyOnBColList.      */

    if (numBCols < 2) boundaryNodes[i] = true;
    for (LO j = 0; j < numBCols; j++) {
      bcol = bColList[j];
      if (keepStatus[bcol] == true) {
        newCols[nzTotal] = bColList[j];
        newRowPtr[i + 1]++;
        nzTotal = nzTotal + 1;
      }
      keepStatus[bcol]        = true;
      alreadyOnBColList[bcol] = false;
      bColList[j]             = 0;
    }
    numBCols = 0;
  } /* for (i = 0; i < nBlks; i++)  */

  /* create array of the correct size and copy over newCols to it */

  Teuchos::ArrayRCP<LO> finalCols(nzTotal);
  for (LO i = 0; i < nzTotal; i++) finalCols[i] = newCols[i];

  // Not using column map because we do not allow for any off-proc stuff.
  // Not sure if this is okay.  FIXME

  RCP<const Map> rowMap = Amat.getRowMap();  // , colMap = Amat.getColMap();

  LO nAmalgNodesOnProc = rowMap->getLocalNumElements() / nPDEs;
  Teuchos::Array<GO> nodalGIDs(nAmalgNodesOnProc);
  typename Teuchos::ScalarTraits<Scalar>::coordinateType temp;
  for (size_t i = 0; i < as<size_t>(nAmalgNodesOnProc); i++) {
    GO gid       = rowMap->getGlobalElement(i * nPDEs);
    temp         = ((typename Teuchos::ScalarTraits<Scalar>::coordinateType)(gid)) / ((typename Teuchos::ScalarTraits<Scalar>::coordinateType)(nPDEs));
    nodalGIDs[i] = as<GO>(floor(temp));
  }
  GO nAmalgNodesGlobal = rowMap->getGlobalNumElements();
  GO nBlkGlobal        = nAmalgNodesGlobal / nPDEs;
  if (nBlkGlobal * nPDEs != nAmalgNodesGlobal)
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Number of global dofs not divisible by BlkSize");

  Teuchos::RCP<Map> AmalgRowMap = MapFactory::Build(rowMap->lib(), nBlkGlobal,
                                                    nodalGIDs(), 0, rowMap->getComm());

  filteredGraph = rcp(new LWGraph(newRowPtr, finalCols, AmalgRowMap, AmalgRowMap, "thresholded graph of A"));
  filteredGraph->SetBoundaryNodeMap(boundaryNodes);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::badGuysDropfunc(LO row, const Teuchos::ArrayView<const LocalOrdinal>& cols, const Teuchos::ArrayView<const Scalar>& vals, const MultiVector& testVecs, LO nPDEs, Teuchos::ArrayRCP<Scalar>& penalties, const MultiVector& nearNull, Teuchos::ArrayRCP<LO>& Bcols, Teuchos::ArrayRCP<bool>& keepOrNot, LO& Nbcols, LO nLoc) const {
  using TST = Teuchos::ScalarTraits<Scalar>;

  LO nLeng = cols.size();
  typename TST::coordinateType temp;
  temp      = ((typename TST::coordinateType)(row)) / ((typename TST::coordinateType)(nPDEs));
  LO blkRow = as<LO>(floor(temp));
  Teuchos::ArrayRCP<Scalar> badGuy(nLeng, 0.0);
  Teuchos::ArrayRCP<Scalar> subNull(nLeng, 0.0); /* subset of nearNull       */
                                                 /* associated with current  */
                                                 /* dof within node.         */

  /* Only consider testVecs associated with same dof & on processor. Further   */
  /* collapse testVecs to a single badGuy vector by basically taking the worst */
  /* (least smooth) values for each of the off diags. In particular, we look at*/
  /* the ratio of each off-diag test value / diag test value and compare this  */
  /* with the nearNull vector ratio. The further the testVec ratio is from the */
  /* nearNull ratio, the harder is will be to accurately interpolate is these  */
  /* two guys are aggregated. So, the biggest ratio mismatch is used to choose */
  /* the testVec entry associated with each off-diagonal entry.                */

  for (LO i = 0; i < nLeng; i++) keepOrNot[i] = false;

  LO diagInd                              = -1;
  Nbcols                                  = 0;
  LO rowDof                               = row - blkRow * nPDEs;
  Teuchos::ArrayRCP<const Scalar> oneNull = nearNull.getData(as<size_t>(rowDof));

  for (LO i = 0; i < nLeng; i++) {
    if ((cols[i] < nLoc) && (TST::magnitude(vals[i]) != 0.0)) { /* on processor */
      temp      = ((typename TST::coordinateType)(cols[i])) / ((typename TST::coordinateType)(nPDEs));
      LO colDof = cols[i] - (as<LO>(floor(temp))) * nPDEs;
      if (colDof == rowDof) { /* same dof within node as row */
        Bcols[Nbcols]   = (cols[i] - colDof) / nPDEs;
        subNull[Nbcols] = oneNull[cols[i]];

        if (cols[i] != row) { /* not diagonal  */
          Scalar worstRatio  = -TST::one();
          Scalar targetRatio = subNull[Nbcols] / oneNull[row];
          Scalar actualRatio;
          for (size_t kk = 0; kk < testVecs.getNumVectors(); kk++) {
            Teuchos::ArrayRCP<const Scalar> curVec = testVecs.getData(kk);
            actualRatio                            = curVec[cols[i]] / curVec[row];
            if (TST::magnitude(actualRatio - targetRatio) > TST::magnitude(worstRatio)) {
              badGuy[Nbcols] = actualRatio;
              worstRatio     = Teuchos::ScalarTraits<SC>::magnitude(actualRatio - targetRatio);
            }
          }
        } else {
          badGuy[Nbcols]    = 1.;
          keepOrNot[Nbcols] = true;
          diagInd           = Nbcols;
        }
        (Nbcols)++;
      }
    }
  }

  /* Make sure that diagonal entry is in block col list */

  if (diagInd == -1) {
    Bcols[Nbcols]     = (row - rowDof) / nPDEs;
    subNull[Nbcols]   = 1.;
    badGuy[Nbcols]    = 1.;
    keepOrNot[Nbcols] = true;
    diagInd           = Nbcols;
    (Nbcols)++;
  }

  Scalar currentRP           = oneNull[row] * oneNull[row];
  Scalar currentRTimesBadGuy = oneNull[row] * badGuy[diagInd];
  Scalar currentScore        = penalties[0]; /* (I - P inv(R*P)*R )=0 for size */
                                             /* size 1 agg, so fit is perfect  */

  /* starting from a set that only includes the diagonal entry consider adding */
  /* one off-diagonal at a time until the fitValue exceeds the penalty term.   */
  /* Here, the fit value is  (I - P inv(R P) R ) and we always consider the    */
  /* lowest fitValue that is not currently in the set. R and P correspond to   */
  /* a simple tentaive grid transfer associated with an aggregate that         */
  /* includes the diagonal, all already determined neighbors, and the potential*/
  /* new neighbor                                                              */

  LO nKeep = 1, flag      = 1, minId;
  Scalar minFit, minFitRP = 0., minFitRTimesBadGuy = 0.;
  Scalar newRP, newRTimesBadGuy;

  while (flag == 1) {
    /* compute a fit for each possible off-diagonal neighbor */
    /* that has not already been added as a neighbor         */

    minFit = 1000000.;
    minId  = -1;

    for (LO i = 0; i < Nbcols; i++) {
      if (keepOrNot[i] == false) {
        keepOrNot[i]    = true; /* temporarily view i as non-dropped neighbor */
        newRP           = currentRP + subNull[i] * subNull[i];
        newRTimesBadGuy = currentRTimesBadGuy + subNull[i] * badGuy[i];
        Scalar ratio    = newRTimesBadGuy / newRP;

        Scalar newFit = 0.0;
        for (LO k = 0; k < Nbcols; k++) {
          if (keepOrNot[k] == true) {
            Scalar diff = badGuy[k] - ratio * subNull[k];
            newFit      = newFit + diff * diff;
          }
        }
        if (Teuchos::ScalarTraits<SC>::magnitude(newFit) < Teuchos::ScalarTraits<SC>::magnitude(minFit)) {
          minId              = i;
          minFit             = newFit;
          minFitRP           = newRP;
          minFitRTimesBadGuy = newRTimesBadGuy;
        }
        keepOrNot[i] = false;
      }
    }
    if (minId == -1)
      flag = 0;
    else {
      minFit          = sqrt(minFit);
      Scalar newScore = penalties[nKeep] + minFit;
      if (Teuchos::ScalarTraits<SC>::magnitude(newScore) < Teuchos::ScalarTraits<SC>::magnitude(currentScore)) {
        nKeep               = nKeep + 1;
        keepOrNot[minId]    = true;
        currentScore        = newScore;
        currentRP           = minFitRP;
        currentRTimesBadGuy = minFitRTimesBadGuy;
      } else
        flag = 0;
    }
  }
}

}  // namespace MueLu

#endif  // MUELU_SMOOVECCOALESCEDROPFACTORY_DEF_HPP
