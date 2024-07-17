// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BLACKBOXPFACTORY_DEF_HPP
#define MUELU_BLACKBOXPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>

// #include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Xpetra_CrsMatrixUtils.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_BlackBoxPFactory_decl.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  // Coarsen can come in two forms, either a single char that will be interpreted as an integer
  // which is used as the coarsening rate in every spatial dimentions,
  // or it can be a longer string that will then be interpreted as an array of integers.
  // By default coarsen is set as "{2}", hence a coarsening rate of 2 in every spatial dimension
  // is the default setting!
  validParamList->set<std::string>("Coarsen", "{3}", "Coarsening rate in all spatial dimensions");
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Generating factory for coorindates");
  validParamList->set<RCP<const FactoryBase> >("gNodesPerDim", Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
  validParamList->set<RCP<const FactoryBase> >("lNodesPerDim", Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
  validParamList->set<std::string>("stencil type", "full", "You can use two type of stencils: full and reduced, that correspond to 27 and 7 points stencils respectively in 3D.");
  validParamList->set<std::string>("block strategy", "coupled", "The strategy used to handle systems of PDEs can be: coupled or uncoupled.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel,
                                                                               Level& /* coarseLevel */)
    const {
  Input(fineLevel, "A");
  Input(fineLevel, "Nullspace");
  Input(fineLevel, "Coordinates");
  // Request the global number of nodes per dimensions
  if (fineLevel.GetLevelID() == 0) {
    if (fineLevel.IsAvailable("gNodesPerDim", NoFactory::get())) {
      fineLevel.DeclareInput("gNodesPerDim", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "gNodesPerDim was not provided by the user on level0!");
    }
  } else {
    Input(fineLevel, "gNodesPerDim");
  }

  // Request the local number of nodes per dimensions
  if (fineLevel.GetLevelID() == 0) {
    if (fineLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
      fineLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "lNodesPerDim was not provided by the user on level0!");
    }
  } else {
    Input(fineLevel, "lNodesPerDim");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel,
                                                                        Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel,
                                                                         Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  // Get parameter list
  const ParameterList& pL = GetParameterList();

  // obtain general variables
  RCP<Matrix> A                  = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<MultiVector> fineNullspace = Get<RCP<MultiVector> >(fineLevel, "Nullspace");
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates =
      Get<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > >(fineLevel, "Coordinates");
  LO numDimensions = coordinates->getNumVectors();
  LO BlkSize       = A->GetFixedBlockSize();

  // Get fine level geometric data: g(l)FineNodesPerDir and g(l)NumFineNodes
  Array<GO> gFineNodesPerDir(3);
  Array<LO> lFineNodesPerDir(3);
  // Get the number of points in each direction
  if (fineLevel.GetLevelID() == 0) {
    gFineNodesPerDir = fineLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
    lFineNodesPerDir = fineLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
  } else {
    // Loading global number of nodes per diretions
    gFineNodesPerDir = Get<Array<GO> >(fineLevel, "gNodesPerDim");

    // Loading local number of nodes per diretions
    lFineNodesPerDir = Get<Array<LO> >(fineLevel, "lNodesPerDim");
  }
  for (LO i = 0; i < 3; ++i) {
    if (gFineNodesPerDir[i] == 0) {
      GetOStream(Runtime0) << "gNodesPerDim in direction " << i << " is set to 1 from 0"
                           << std::endl;
      gFineNodesPerDir[i] = 1;
    }
    if (lFineNodesPerDir[i] == 0) {
      GetOStream(Runtime0) << "lNodesPerDim in direction " << i << " is set to 1 from 0"
                           << std::endl;
      lFineNodesPerDir[i] = 1;
    }
  }
  GO gNumFineNodes = gFineNodesPerDir[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
  LO lNumFineNodes = lFineNodesPerDir[2] * lFineNodesPerDir[1] * lFineNodesPerDir[0];

  // Get the coarsening rate
  std::string coarsenRate = pL.get<std::string>("Coarsen");
  Array<LO> coarseRate(3);
  {
    Teuchos::Array<LO> crates;
    try {
      crates = Teuchos::fromStringToArray<LO>(coarsenRate);
    } catch (const Teuchos::InvalidArrayStringRepresentation& e) {
      GetOStream(Errors, -1) << " *** Coarsen must be a string convertible into an array! *** "
                             << std::endl;
      throw e;
    }
    TEUCHOS_TEST_FOR_EXCEPTION((crates.size() > 1) && (crates.size() < numDimensions),
                               Exceptions::RuntimeError,
                               "Coarsen must have at least as many components as the number of"
                               " spatial dimensions in the problem.");
    for (LO i = 0; i < 3; ++i) {
      if (i < numDimensions) {
        if (crates.size() == 1) {
          coarseRate[i] = crates[0];
        } else if (i < crates.size()) {
          coarseRate[i] = crates[i];
        } else {
          GetOStream(Errors, -1) << " *** Coarsen must be at least as long as the number of"
                                    " spatial dimensions! *** "
                                 << std::endl;
          throw Exceptions::RuntimeError(
              " *** Coarsen must be at least as long as the number of"
              " spatial dimensions! *** \n");
        }
      } else {
        coarseRate[i] = 1;
      }
    }
  }  // End of scope for crates

  // Get the stencil type used for discretization
  const std::string stencilType = pL.get<std::string>("stencil type");
  if (stencilType != "full" && stencilType != "reduced") {
    GetOStream(Errors, -1) << " *** stencil type must be set to: full or reduced, any other value "
                              "is trated as an error! *** "
                           << std::endl;
    throw Exceptions::RuntimeError(" *** stencil type is neither full, nor reduced! *** \n");
  }

  // Get the strategy for PDE systems
  const std::string blockStrategy = pL.get<std::string>("block strategy");
  if (blockStrategy != "coupled" && blockStrategy != "uncoupled") {
    GetOStream(Errors, -1) << " *** block strategy must be set to: coupled or uncoupled, any other "
                              "value is trated as an error! *** "
                           << std::endl;
    throw Exceptions::RuntimeError(" *** block strategy is neither coupled, nor uncoupled! *** \n");
  }

  GO gNumCoarseNodes = 0;
  LO lNumCoarseNodes = 0;
  Array<GO> gIndices(3), gCoarseNodesPerDir(3), ghostGIDs, coarseNodesGIDs, colGIDs;
  Array<LO> myOffset(3), lCoarseNodesPerDir(3), glCoarseNodesPerDir(3), endRate(3);
  Array<bool> ghostInterface(6);
  Array<int> boundaryFlags(3);
  ArrayRCP<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coarseNodes(numDimensions);
  Array<ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > fineNodes(numDimensions);
  for (LO dim = 0; dim < numDimensions; ++dim) {
    fineNodes[dim] = coordinates->getData(dim)();
  }

  // This struct stores PIDs, LIDs and GIDs on the fine mesh and GIDs on the coarse mesh.
  RCP<NodesIDs> ghostedCoarseNodes = rcp(new NodesIDs{});

  GetGeometricData(coordinates, coarseRate, gFineNodesPerDir, lFineNodesPerDir, BlkSize,  // inputs
                   gIndices, myOffset, ghostInterface, endRate, gCoarseNodesPerDir,       // outputs
                   lCoarseNodesPerDir, glCoarseNodesPerDir, ghostGIDs, coarseNodesGIDs, colGIDs,
                   gNumCoarseNodes, lNumCoarseNodes, coarseNodes, boundaryFlags,
                   ghostedCoarseNodes);

  // Create the MultiVector of coarse coordinates
  Xpetra::UnderlyingLib lib      = coordinates->getMap()->lib();
  RCP<const Map> coarseCoordsMap = MapFactory::Build(lib,
                                                     gNumCoarseNodes,
                                                     coarseNodesGIDs.view(0, lNumCoarseNodes),
                                                     coordinates->getMap()->getIndexBase(),
                                                     coordinates->getMap()->getComm());
  Array<ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coarseCoords(numDimensions);
  for (LO dim = 0; dim < numDimensions; ++dim) {
    coarseCoords[dim] = coarseNodes[dim]();
  }
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coarseCoordinates =
      Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(coarseCoordsMap, coarseCoords(),
                                                                                                           numDimensions);

  // Now create a new matrix: Aghost that contains all the data
  // locally needed to compute the local part of the prolongator.
  // Here assuming that all the coarse nodes o and fine nodes +
  // are local then all the data associated with the coarse
  // nodes O and the fine nodes * needs to be imported.
  //
  //                  *--*--*--*--*--*--*--*
  //                  |  |  |  |  |  |  |  |
  //                  o--+--+--o--+--+--O--*
  //                  |  |  |  |  |  |  |  |
  //                  +--+--+--+--+--+--*--*
  //                  |  |  |  |  |  |  |  |
  //                  +--+--+--+--+--+--*--*
  //                  |  |  |  |  |  |  |  |
  //                  o--+--+--o--+--+--O--*
  //                  |  |  |  |  |  |  |  |
  //                  +--+--+--+--+--+--*--*
  //                  |  |  |  |  |  |  |  |
  //                  *--*--*--*--*--*--*--*
  //                  |  |  |  |  |  |  |  |
  //                  O--*--*--O--*--*--O--*
  //
  // Creating that local matrix is easy enough using proper range
  // and domain maps to import data from A. Note that with this
  // approach we reorder the local entries using the domain map and
  // can subsequently compute the prolongator without reordering.
  // As usual we need to be careful about any coarsening rate
  // change at the boundary!

  // The ingredients needed are an importer, a range map and a domain map
  Array<GO> ghostRowGIDs, ghostColGIDs, nodeSteps(3);
  nodeSteps[0] = 1;
  nodeSteps[1] = gFineNodesPerDir[0];
  nodeSteps[2] = gFineNodesPerDir[0] * gFineNodesPerDir[1];
  Array<LO> glFineNodesPerDir(3);
  GO startingGID = A->getRowMap()->getMinGlobalIndex();
  for (LO dim = 0; dim < 3; ++dim) {
    LO numCoarseNodes = 0;
    if (dim < numDimensions) {
      startingGID -= myOffset[dim] * nodeSteps[dim];
      numCoarseNodes = lCoarseNodesPerDir[dim];
      if (ghostInterface[2 * dim]) {
        ++numCoarseNodes;
      }
      if (ghostInterface[2 * dim + 1]) {
        ++numCoarseNodes;
      }
      if (gIndices[dim] + lFineNodesPerDir[dim] == gFineNodesPerDir[dim]) {
        glFineNodesPerDir[dim] = (numCoarseNodes - 2) * coarseRate[dim] + endRate[dim] + 1;
      } else {
        glFineNodesPerDir[dim] = (numCoarseNodes - 1) * coarseRate[dim] + 1;
      }
    } else {
      glFineNodesPerDir[dim] = 1;
    }
  }
  ghostRowGIDs.resize(glFineNodesPerDir[0] * glFineNodesPerDir[1] * glFineNodesPerDir[2] * BlkSize);
  for (LO k = 0; k < glFineNodesPerDir[2]; ++k) {
    for (LO j = 0; j < glFineNodesPerDir[1]; ++j) {
      for (LO i = 0; i < glFineNodesPerDir[0]; ++i) {
        for (LO l = 0; l < BlkSize; ++l) {
          ghostRowGIDs[(k * glFineNodesPerDir[1] * glFineNodesPerDir[0] + j * glFineNodesPerDir[0] + i) * BlkSize + l] = startingGID + (k * gFineNodesPerDir[1] * gFineNodesPerDir[0] + j * gFineNodesPerDir[0] + i) * BlkSize + l;
        }
      }
    }
  }

  // Looking at the above loops it is easy to find startingGID for the ghostColGIDs
  Array<GO> startingGlobalIndices(numDimensions), dimStride(numDimensions);
  Array<GO> startingColIndices(numDimensions), finishingColIndices(numDimensions);
  GO colMinGID = 0;
  Array<LO> colRange(numDimensions);
  dimStride[0] = 1;
  for (int dim = 1; dim < numDimensions; ++dim) {
    dimStride[dim] = dimStride[dim - 1] * gFineNodesPerDir[dim - 1];
  }
  {
    GO tmp = startingGID;
    for (int dim = numDimensions; dim > 0; --dim) {
      startingGlobalIndices[dim - 1] = tmp / dimStride[dim - 1];
      tmp                            = tmp % dimStride[dim - 1];

      if (startingGlobalIndices[dim - 1] > 0) {
        startingColIndices[dim - 1] = startingGlobalIndices[dim - 1] - 1;
      }
      if (startingGlobalIndices[dim - 1] + glFineNodesPerDir[dim - 1] < gFineNodesPerDir[dim - 1]) {
        finishingColIndices[dim - 1] = startingGlobalIndices[dim - 1] + glFineNodesPerDir[dim - 1];
      } else {
        finishingColIndices[dim - 1] = startingGlobalIndices[dim - 1] + glFineNodesPerDir[dim - 1] - 1;
      }
      colRange[dim - 1] = finishingColIndices[dim - 1] - startingColIndices[dim - 1] + 1;
      colMinGID += startingColIndices[dim - 1] * dimStride[dim - 1];
    }
  }
  ghostColGIDs.resize(colRange[0] * colRange[1] * colRange[2] * BlkSize);
  for (LO k = 0; k < colRange[2]; ++k) {
    for (LO j = 0; j < colRange[1]; ++j) {
      for (LO i = 0; i < colRange[0]; ++i) {
        for (LO l = 0; l < BlkSize; ++l) {
          ghostColGIDs[(k * colRange[1] * colRange[0] + j * colRange[0] + i) * BlkSize + l] = colMinGID + (k * gFineNodesPerDir[1] * gFineNodesPerDir[0] + j * gFineNodesPerDir[0] + i) * BlkSize + l;
        }
      }
    }
  }

  RCP<const Map> ghostedRowMap    = Xpetra::MapFactory<LO, GO, NO>::Build(A->getRowMap()->lib(),
                                                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                          ghostRowGIDs(),
                                                                          A->getRowMap()->getIndexBase(),
                                                                          A->getRowMap()->getComm());
  RCP<const Map> ghostedColMap    = Xpetra::MapFactory<LO, GO, NO>::Build(A->getRowMap()->lib(),
                                                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                          ghostColGIDs(),
                                                                          A->getRowMap()->getIndexBase(),
                                                                          A->getRowMap()->getComm());
  RCP<const Import> ghostImporter = Xpetra::ImportFactory<LO, GO, NO>::Build(A->getRowMap(),
                                                                             ghostedRowMap);
  RCP<const Matrix> Aghost        = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(A, *ghostImporter,
                                                                                 ghostedRowMap,
                                                                                 ghostedRowMap);

  // Create the maps and data structures for the projection matrix
  RCP<const Map> rowMapP = A->getDomainMap();

  RCP<const Map> domainMapP;

  RCP<const Map> colMapP;

  // At this point we need to create the column map which is a delicate operation.
  // The entries in that map need to be ordered as follows:
  //         1) first owned entries ordered by LID
  //         2) second order the remaining entries by PID
  //         3) entries with the same remote PID are ordered by GID.
  // One piece of good news: lNumCoarseNodes is the number of ownedNodes and lNumGhostNodes
  // is the number of remote nodes. The sorting can be limited to remote nodes
  // as the owned ones are alreaded ordered by LID!

  LO lNumGhostedNodes = ghostedCoarseNodes->GIDs.size();
  {
    std::vector<NodeID> colMapOrdering(lNumGhostedNodes);
    for (LO ind = 0; ind < lNumGhostedNodes; ++ind) {
      colMapOrdering[ind].GID = ghostedCoarseNodes->GIDs[ind];
      if (ghostedCoarseNodes->PIDs[ind] == rowMapP->getComm()->getRank()) {
        colMapOrdering[ind].PID = -1;
      } else {
        colMapOrdering[ind].PID = ghostedCoarseNodes->PIDs[ind];
      }
      colMapOrdering[ind].LID     = ghostedCoarseNodes->LIDs[ind];
      colMapOrdering[ind].lexiInd = ind;
    }
    std::sort(colMapOrdering.begin(), colMapOrdering.end(),
              [](NodeID a, NodeID b) -> bool {
                return ((a.PID < b.PID) || ((a.PID == b.PID) && (a.GID < b.GID)));
              });

    colGIDs.resize(BlkSize * lNumGhostedNodes);
    for (LO ind = 0; ind < lNumGhostedNodes; ++ind) {
      // Save the permutation calculated to go from Lexicographic indexing to column map indexing
      ghostedCoarseNodes->colInds[colMapOrdering[ind].lexiInd] = ind;
      for (LO dof = 0; dof < BlkSize; ++dof) {
        colGIDs[BlkSize * ind + dof] = BlkSize * colMapOrdering[ind].GID + dof;
      }
    }
    domainMapP = Xpetra::MapFactory<LO, GO, NO>::Build(rowMapP->lib(),
                                                       BlkSize * gNumCoarseNodes,
                                                       colGIDs.view(0, BlkSize * lNumCoarseNodes),
                                                       rowMapP->getIndexBase(),
                                                       rowMapP->getComm());
    colMapP    = Xpetra::MapFactory<LO, GO, NO>::Build(rowMapP->lib(),
                                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                       colGIDs.view(0, colGIDs.size()),
                                                       rowMapP->getIndexBase(),
                                                       rowMapP->getComm());
  }  // End of scope for colMapOrdering and colGIDs

  std::vector<size_t> strideInfo(1);
  strideInfo[0]                    = BlkSize;
  RCP<const Map> stridedDomainMapP = Xpetra::StridedMapFactory<LO, GO, NO>::Build(domainMapP,
                                                                                  strideInfo);

  GO gnnzP = 0;
  LO lnnzP = 0;
  // coarse points have one nnz per row
  gnnzP += gNumCoarseNodes;
  lnnzP += lNumCoarseNodes;
  // add all other points multiplying by 2^numDimensions
  gnnzP += (gNumFineNodes - gNumCoarseNodes) * std::pow(2, numDimensions);
  lnnzP += (lNumFineNodes - lNumCoarseNodes) * std::pow(2, numDimensions);
  // finally multiply by the number of dofs per node
  gnnzP = gnnzP * BlkSize;
  lnnzP = lnnzP * BlkSize;

  // Create the matrix itself using the above maps
  RCP<Matrix> P;
  P                   = rcp(new CrsMatrixWrap(rowMapP, colMapP, 0));
  RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  ArrayRCP<size_t> iaP;
  ArrayRCP<LO> jaP;
  ArrayRCP<SC> valP;

  PCrs->allocateAllValues(lnnzP, iaP, jaP, valP);

  ArrayView<size_t> ia = iaP();
  ArrayView<LO> ja     = jaP();
  ArrayView<SC> val    = valP();
  ia[0]                = 0;

  LO numCoarseElements = 1;
  Array<LO> lCoarseElementsPerDir(3);
  for (LO dim = 0; dim < numDimensions; ++dim) {
    lCoarseElementsPerDir[dim] = lCoarseNodesPerDir[dim];
    if (ghostInterface[2 * dim]) {
      ++lCoarseElementsPerDir[dim];
    }
    if (!ghostInterface[2 * dim + 1]) {
      --lCoarseElementsPerDir[dim];
    }
    numCoarseElements = numCoarseElements * lCoarseElementsPerDir[dim];
  }

  for (LO dim = numDimensions; dim < 3; ++dim) {
    lCoarseElementsPerDir[dim] = 1;
  }

  // Loop over the coarse elements
  Array<int> elementFlags(3);
  Array<LO> elemInds(3), elementNodesPerDir(3), glElementRefTuple(3);
  Array<LO> glElementRefTupleCG(3), glElementCoarseNodeCG(8);
  const int numCoarseNodesInElement = std::pow(2, numDimensions);
  const int nnzPerCoarseNode        = (blockStrategy == "coupled") ? BlkSize : 1;
  const int numRowsPerPoint         = BlkSize;
  for (elemInds[2] = 0; elemInds[2] < lCoarseElementsPerDir[2]; ++elemInds[2]) {
    for (elemInds[1] = 0; elemInds[1] < lCoarseElementsPerDir[1]; ++elemInds[1]) {
      for (elemInds[0] = 0; elemInds[0] < lCoarseElementsPerDir[0]; ++elemInds[0]) {
        elementFlags[0] = 0;
        elementFlags[1] = 0;
        elementFlags[2] = 0;
        for (int dim = 0; dim < 3; ++dim) {
          // Detect boundary conditions on the element and set corresponding flags.
          if (elemInds[dim] == 0 && elemInds[dim] == lCoarseElementsPerDir[dim] - 1) {
            elementFlags[dim] = boundaryFlags[dim];
          } else if (elemInds[dim] == 0 && (boundaryFlags[dim] == 1 || boundaryFlags[dim] == 3)) {
            elementFlags[dim] += 1;
          } else if ((elemInds[dim] == lCoarseElementsPerDir[dim] - 1) && (boundaryFlags[dim] == 2 || boundaryFlags[dim] == 3)) {
            elementFlags[dim] += 2;
          } else {
            elementFlags[dim] = 0;
          }

          // Compute the number of nodes in the current element.
          if (dim < numDimensions) {
            if ((elemInds[dim] == lCoarseElementsPerDir[dim]) && (gIndices[dim] + lFineNodesPerDir[dim] == gFineNodesPerDir[dim])) {
              elementNodesPerDir[dim] = endRate[dim] + 1;
            } else {
              elementNodesPerDir[dim] = coarseRate[dim] + 1;
            }
          } else {
            elementNodesPerDir[dim] = 1;
          }

          // Get the lowest tuple of the element using the ghosted local coordinate system
          glElementRefTuple[dim]   = elemInds[dim] * coarseRate[dim];
          glElementRefTupleCG[dim] = elemInds[dim];
        }

        // Now get the column map indices corresponding to the dofs associated with the current
        // element's coarse nodes.
        for (typename Array<LO>::size_type elem = 0; elem < glElementCoarseNodeCG.size(); ++elem) {
          glElementCoarseNodeCG[elem] = glElementRefTupleCG[2] * glCoarseNodesPerDir[1] * glCoarseNodesPerDir[0] + glElementRefTupleCG[1] * glCoarseNodesPerDir[0] + glElementRefTupleCG[0];
        }
        glElementCoarseNodeCG[4] += glCoarseNodesPerDir[1] * glCoarseNodesPerDir[0];
        glElementCoarseNodeCG[5] += glCoarseNodesPerDir[1] * glCoarseNodesPerDir[0];
        glElementCoarseNodeCG[6] += glCoarseNodesPerDir[1] * glCoarseNodesPerDir[0];
        glElementCoarseNodeCG[7] += glCoarseNodesPerDir[1] * glCoarseNodesPerDir[0];

        glElementCoarseNodeCG[2] += glCoarseNodesPerDir[0];
        glElementCoarseNodeCG[3] += glCoarseNodesPerDir[0];
        glElementCoarseNodeCG[6] += glCoarseNodesPerDir[0];
        glElementCoarseNodeCG[7] += glCoarseNodesPerDir[0];

        glElementCoarseNodeCG[1] += 1;
        glElementCoarseNodeCG[3] += 1;
        glElementCoarseNodeCG[5] += 1;
        glElementCoarseNodeCG[7] += 1;

        LO numNodesInElement = elementNodesPerDir[0] * elementNodesPerDir[1] * elementNodesPerDir[2];
        // LO elementOffset = elemInds[2]*coarseRate[2]*glFineNodesPerDir[1]*glFineNodesPerDir[0]
        //   + elemInds[1]*coarseRate[1]*glFineNodesPerDir[0] + elemInds[0]*coarseRate[0];

        // Compute the element prolongator
        Teuchos::SerialDenseMatrix<LO, SC> Pi, Pf, Pe;
        Array<LO> dofType(numNodesInElement * BlkSize), lDofInd(numNodesInElement * BlkSize);
        ComputeLocalEntries(Aghost, coarseRate, endRate, BlkSize, elemInds, lCoarseElementsPerDir,
                            numDimensions, glFineNodesPerDir, gFineNodesPerDir, gIndices,
                            lCoarseNodesPerDir, ghostInterface, elementFlags, stencilType,
                            blockStrategy, elementNodesPerDir, numNodesInElement, colGIDs,
                            Pi, Pf, Pe, dofType, lDofInd);

        // Find ghosted LID associated with nodes in the element and eventually which of these
        // nodes are ghosts, this information is used to fill the local prolongator.
        Array<LO> lNodeLIDs(numNodesInElement);
        {
          Array<LO> lNodeTuple(3), nodeInd(3);
          for (nodeInd[2] = 0; nodeInd[2] < elementNodesPerDir[2]; ++nodeInd[2]) {
            for (nodeInd[1] = 0; nodeInd[1] < elementNodesPerDir[1]; ++nodeInd[1]) {
              for (nodeInd[0] = 0; nodeInd[0] < elementNodesPerDir[0]; ++nodeInd[0]) {
                int stencilLength = 0;
                if ((nodeInd[0] == 0 || nodeInd[0] == elementNodesPerDir[0] - 1) &&
                    (nodeInd[1] == 0 || nodeInd[1] == elementNodesPerDir[1] - 1) &&
                    (nodeInd[2] == 0 || nodeInd[2] == elementNodesPerDir[2] - 1)) {
                  stencilLength = 1;
                } else {
                  stencilLength = std::pow(2, numDimensions);
                }
                LO nodeElementInd = nodeInd[2] * elementNodesPerDir[1] * elementNodesPerDir[1] + nodeInd[1] * elementNodesPerDir[0] + nodeInd[0];
                for (int dim = 0; dim < 3; ++dim) {
                  lNodeTuple[dim] = glElementRefTuple[dim] - myOffset[dim] + nodeInd[dim];
                }
                if (lNodeTuple[0] < 0 || lNodeTuple[1] < 0 || lNodeTuple[2] < 0 ||
                    lNodeTuple[0] > lFineNodesPerDir[0] - 1 ||
                    lNodeTuple[1] > lFineNodesPerDir[1] - 1 ||
                    lNodeTuple[2] > lFineNodesPerDir[2] - 1) {
                  // This flags the ghosts nodes used for prolongator calculation but for which
                  // we do not own the row, hence we won't fill these values on this rank.
                  lNodeLIDs[nodeElementInd] = -1;
                } else if ((nodeInd[0] == 0 && elemInds[0] != 0) ||
                           (nodeInd[1] == 0 && elemInds[1] != 0) ||
                           (nodeInd[2] == 0 && elemInds[2] != 0)) {
                  // This flags nodes that are owned but common to two coarse elements and that
                  // were already filled by another element, we don't want to fill them twice so
                  // we skip them
                  lNodeLIDs[nodeElementInd] = -2;
                } else {
                  // The remaining nodes are locally owned and we need to fill the coresponding
                  // rows of the prolongator

                  // First we need to find which row needs to be filled
                  lNodeLIDs[nodeElementInd] = lNodeTuple[2] * lFineNodesPerDir[1] * lFineNodesPerDir[0] + lNodeTuple[1] * lFineNodesPerDir[0] + lNodeTuple[0];

                  // We also compute the row offset using arithmetic to ensure that we can loop
                  // easily over the nodes in a macro-element as well as facilitate on-node
                  // parallelization. The node serial code could be rewritten with two loops over
                  // the local part of the mesh to avoid the costly integer divisions...
                  Array<LO> refCoarsePointTuple(3);
                  for (int dim = 2; dim > -1; --dim) {
                    if (dim == 0) {
                      refCoarsePointTuple[dim] = (lNodeTuple[dim] + myOffset[dim]) / coarseRate[dim];
                      if (myOffset[dim] == 0) {
                        ++refCoarsePointTuple[dim];
                      }
                    } else {
                      // Note:  no need for magnitudeType here, just use double because these things are LO's
                      refCoarsePointTuple[dim] =
                          std::ceil(static_cast<double>(lNodeTuple[dim] + myOffset[dim]) / coarseRate[dim]);
                    }
                    if ((lNodeTuple[dim] + myOffset[dim]) % coarseRate[dim] > 0) {
                      break;
                    }
                  }
                  const LO numCoarsePoints = refCoarsePointTuple[0] + refCoarsePointTuple[1] * lCoarseNodesPerDir[0] + refCoarsePointTuple[2] * lCoarseNodesPerDir[1] * lCoarseNodesPerDir[0];
                  const LO numFinePoints   = lNodeLIDs[nodeElementInd] + 1;

                  // The below formula computes the rowPtr for row 'n+BlkSize' and we are about to
                  // fill row 'n' to 'n+BlkSize'.
                  size_t rowPtr = (numFinePoints - numCoarsePoints) * numRowsPerPoint * numCoarseNodesInElement * nnzPerCoarseNode + numCoarsePoints * numRowsPerPoint;
                  if (dofType[nodeElementInd * BlkSize] == 0) {
                    // Check if current node is a coarse point
                    rowPtr = rowPtr - numRowsPerPoint;
                  } else {
                    rowPtr = rowPtr - numRowsPerPoint * numCoarseNodesInElement * nnzPerCoarseNode;
                  }

                  for (int dof = 0; dof < BlkSize; ++dof) {
                    // Now we loop over the stencil and fill the column indices and row values
                    int cornerInd = 0;
                    switch (dofType[nodeElementInd * BlkSize + dof]) {
                      case 0:  // Corner node
                        if (nodeInd[2] == elementNodesPerDir[2] - 1) {
                          cornerInd += 4;
                        }
                        if (nodeInd[1] == elementNodesPerDir[1] - 1) {
                          cornerInd += 2;
                        }
                        if (nodeInd[0] == elementNodesPerDir[0] - 1) {
                          cornerInd += 1;
                        }
                        ia[lNodeLIDs[nodeElementInd] * BlkSize + dof + 1] = rowPtr + dof + 1;
                        ja[rowPtr + dof]                                  = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[cornerInd]] * BlkSize + dof;
                        val[rowPtr + dof]                                 = 1.0;
                        break;

                      case 1:  // Edge node
                        ia[lNodeLIDs[nodeElementInd] * BlkSize + dof + 1] = rowPtr + (dof + 1) * numCoarseNodesInElement * nnzPerCoarseNode;
                        for (int ind1 = 0; ind1 < stencilLength; ++ind1) {
                          if (blockStrategy == "coupled") {
                            for (int ind2 = 0; ind2 < BlkSize; ++ind2) {
                              size_t lRowPtr = rowPtr + dof * numCoarseNodesInElement * nnzPerCoarseNode + ind1 * BlkSize + ind2;
                              ja[lRowPtr]    = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[ind1]] * BlkSize + ind2;
                              val[lRowPtr]   = Pe(lDofInd[nodeElementInd * BlkSize + dof],
                                                  ind1 * BlkSize + ind2);
                            }
                          } else if (blockStrategy == "uncoupled") {
                            size_t lRowPtr   = rowPtr + dof * numCoarseNodesInElement * nnzPerCoarseNode + ind1;
                            ja[rowPtr + dof] = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[ind1]] * BlkSize + dof;
                            val[lRowPtr]     = Pe(lDofInd[nodeElementInd * BlkSize + dof],
                                                  ind1 * BlkSize + dof);
                          }
                        }
                        break;

                      case 2:  // Face node
                        ia[lNodeLIDs[nodeElementInd] * BlkSize + dof + 1] = rowPtr + (dof + 1) * numCoarseNodesInElement * nnzPerCoarseNode;
                        for (int ind1 = 0; ind1 < stencilLength; ++ind1) {
                          if (blockStrategy == "coupled") {
                            for (int ind2 = 0; ind2 < BlkSize; ++ind2) {
                              size_t lRowPtr = rowPtr + dof * numCoarseNodesInElement * nnzPerCoarseNode + ind1 * BlkSize + ind2;
                              ja[lRowPtr]    = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[ind1]] * BlkSize + ind2;
                              val[lRowPtr]   = Pf(lDofInd[nodeElementInd * BlkSize + dof],
                                                  ind1 * BlkSize + ind2);
                            }
                          } else if (blockStrategy == "uncoupled") {
                            size_t lRowPtr = rowPtr + dof * numCoarseNodesInElement * nnzPerCoarseNode + ind1;
                            // ja [lRowPtr] = colGIDs[glElementCoarseNodeCG[ind1]*BlkSize + dof];
                            ja[lRowPtr]  = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[ind1]] * BlkSize + dof;
                            val[lRowPtr] = Pf(lDofInd[nodeElementInd * BlkSize + dof],
                                              ind1 * BlkSize + dof);
                          }
                        }
                        break;

                      case 3:  // Interior node
                        ia[lNodeLIDs[nodeElementInd] * BlkSize + dof + 1] = rowPtr + (dof + 1) * numCoarseNodesInElement * nnzPerCoarseNode;
                        for (int ind1 = 0; ind1 < stencilLength; ++ind1) {
                          if (blockStrategy == "coupled") {
                            for (int ind2 = 0; ind2 < BlkSize; ++ind2) {
                              size_t lRowPtr = rowPtr + dof * numCoarseNodesInElement * nnzPerCoarseNode + ind1 * BlkSize + ind2;
                              ja[lRowPtr]    = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[ind1]] * BlkSize + ind2;
                              val[lRowPtr]   = Pi(lDofInd[nodeElementInd * BlkSize + dof],
                                                  ind1 * BlkSize + ind2);
                            }
                          } else if (blockStrategy == "uncoupled") {
                            size_t lRowPtr = rowPtr + dof * numCoarseNodesInElement * nnzPerCoarseNode + ind1;
                            ja[lRowPtr]    = ghostedCoarseNodes->colInds[glElementCoarseNodeCG[ind1]] * BlkSize + dof;
                            val[lRowPtr]   = Pi(lDofInd[nodeElementInd * BlkSize + dof],
                                                ind1 * BlkSize + dof);
                          }
                        }
                        break;
                    }
                  }
                }
              }
            }
          }
        }  // End of scopt for lNodeTuple and nodeInd
      }
    }
  }

  // Sort all row's column indicies and entries by LID
  Xpetra::CrsMatrixUtils<SC, LO, GO, NO>::sortCrsEntries(ia, ja, val, rowMapP->lib());

  // Set the values of the prolongation operators into the CrsMatrix P and call FillComplete
  PCrs->setAllValues(iaP, jaP, valP);
  PCrs->expertStaticFillComplete(domainMapP, rowMapP);

  // set StridingInformation of P
  if (A->IsView("stridedMaps") == true) {
    P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), stridedDomainMapP);
  } else {
    P->CreateView("stridedMaps", P->getRangeMap(), stridedDomainMapP);
  }

  // store the transfer operator and node coordinates on coarse level
  Set(coarseLevel, "P", P);
  Set(coarseLevel, "coarseCoordinates", coarseCoordinates);
  Set<Array<GO> >(coarseLevel, "gCoarseNodesPerDim", gCoarseNodesPerDir);
  Set<Array<LO> >(coarseLevel, "lCoarseNodesPerDim", lCoarseNodesPerDir);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetGeometricData(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> >& coordinates,
                     const Array<LO> coarseRate, const Array<GO> gFineNodesPerDir,
                     const Array<LO> lFineNodesPerDir, const LO BlkSize, Array<GO>& gIndices,
                     Array<LO>& myOffset, Array<bool>& ghostInterface, Array<LO>& endRate,
                     Array<GO>& gCoarseNodesPerDir, Array<LO>& lCoarseNodesPerDir,
                     Array<LO>& glCoarseNodesPerDir, Array<GO>& ghostGIDs, Array<GO>& coarseNodesGIDs,
                     Array<GO>& colGIDs, GO& gNumCoarseNodes, LO& lNumCoarseNodes,
                     ArrayRCP<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coarseNodes, Array<int>& boundaryFlags,
                     RCP<NodesIDs> ghostedCoarseNodes) const {
  // This function is extracting the geometric information from the coordinates
  // and creates the necessary data/formatting to perform locally the calculation
  // of the pronlongator.
  //
  // inputs:

  RCP<const Map> coordinatesMap = coordinates->getMap();
  LO numDimensions              = coordinates->getNumVectors();

  // Using the coarsening rate and the fine level data,
  // compute coarse level data

  //                              Phase 1                               //
  // ------------------------------------------------------------------ //
  // We first start by finding small informations on the mesh such as   //
  // the number of coarse nodes (local and global) and the number of    //
  // ghost nodes / the end rate of coarsening.                          //
  // ------------------------------------------------------------------ //
  GO minGlobalIndex = coordinatesMap->getMinGlobalIndex();
  {
    GO tmp;
    gIndices[2] = minGlobalIndex / (gFineNodesPerDir[1] * gFineNodesPerDir[0]);
    tmp         = minGlobalIndex % (gFineNodesPerDir[1] * gFineNodesPerDir[0]);
    gIndices[1] = tmp / gFineNodesPerDir[0];
    gIndices[0] = tmp % gFineNodesPerDir[0];

    myOffset[2] = gIndices[2] % coarseRate[2];
    myOffset[1] = gIndices[1] % coarseRate[1];
    myOffset[0] = gIndices[0] % coarseRate[0];
  }

  for (int dim = 0; dim < 3; ++dim) {
    if (gIndices[dim] == 0) {
      boundaryFlags[dim] += 1;
    }
    if (gIndices[dim] + lFineNodesPerDir[dim] == gFineNodesPerDir[dim]) {
      boundaryFlags[dim] += 2;
    }
  }

  // Check whether ghost nodes are needed in each direction
  for (LO i = 0; i < numDimensions; ++i) {
    if ((gIndices[i] != 0) && (gIndices[i] % coarseRate[i] > 0)) {
      ghostInterface[2 * i] = true;
    }
    if (((gIndices[i] + lFineNodesPerDir[i]) != gFineNodesPerDir[i]) && ((gIndices[i] + lFineNodesPerDir[i] - 1) % coarseRate[i] > 0)) {
      ghostInterface[2 * i + 1] = true;
    }
  }

  for (LO i = 0; i < 3; ++i) {
    if (i < numDimensions) {
      lCoarseNodesPerDir[i] = (lFineNodesPerDir[i] + myOffset[i] - 1) / coarseRate[i];
      if (myOffset[i] == 0) {
        ++lCoarseNodesPerDir[i];
      }
      gCoarseNodesPerDir[i] = (gFineNodesPerDir[i] - 1) / coarseRate[i];
      endRate[i]            = (gFineNodesPerDir[i] - 1) % coarseRate[i];
      if (endRate[i] == 0) {
        ++gCoarseNodesPerDir[i];
        endRate[i] = coarseRate[i];
      }
    } else {
      // Most quantities need to be set to 1 for extra dimensions
      // this is rather logical, an x-y plane is like a single layer
      // of nodes in the z direction...
      gCoarseNodesPerDir[i] = 1;
      lCoarseNodesPerDir[i] = 1;
      endRate[i]            = 1;
    }
  }

  gNumCoarseNodes = gCoarseNodesPerDir[0] * gCoarseNodesPerDir[1] * gCoarseNodesPerDir[2];
  lNumCoarseNodes = lCoarseNodesPerDir[0] * lCoarseNodesPerDir[1] * lCoarseNodesPerDir[2];

  // For each direction, determine how many ghost points are required.
  LO lNumGhostNodes = 0;
  Array<GO> startGhostedCoarseNode(3);
  {
    const int complementaryIndices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    LO tmp                               = 0;
    for (LO i = 0; i < 3; ++i) {
      // The first branch of this if-statement will be used if the rank contains only one layer
      // of nodes in direction i, that layer must also coincide with the boundary of the mesh
      // and coarseRate[i] == endRate[i]...
      if ((gIndices[i] == gFineNodesPerDir[i] - 1) && (gIndices[i] % coarseRate[i] == 0)) {
        startGhostedCoarseNode[i] = gIndices[i] / coarseRate[i] - 1;
      } else {
        startGhostedCoarseNode[i] = gIndices[i] / coarseRate[i];
      }
      // First we load the number of locally owned coarse nodes
      glCoarseNodesPerDir[i] = lCoarseNodesPerDir[i];

      // Check whether a face in direction i needs ghost nodes
      if (ghostInterface[2 * i] || ghostInterface[2 * i + 1]) {
        ++glCoarseNodesPerDir[i];
        if (i == 0) {
          tmp = lCoarseNodesPerDir[1] * lCoarseNodesPerDir[2];
        }
        if (i == 1) {
          tmp = lCoarseNodesPerDir[0] * lCoarseNodesPerDir[2];
        }
        if (i == 2) {
          tmp = lCoarseNodesPerDir[0] * lCoarseNodesPerDir[1];
        }
      }
      // If both faces in direction i need nodes, double the number of ghost nodes
      if (ghostInterface[2 * i] && ghostInterface[2 * i + 1]) {
        ++glCoarseNodesPerDir[i];
        tmp = 2 * tmp;
      }
      lNumGhostNodes += tmp;

      // The corners and edges need to be checked in 2D / 3D to add more ghosts nodes
      for (LO j = 0; j < 2; ++j) {
        for (LO k = 0; k < 2; ++k) {
          // Check if two adjoining faces need ghost nodes and then add their common edge
          if (ghostInterface[2 * complementaryIndices[i][0] + j] && ghostInterface[2 * complementaryIndices[i][1] + k]) {
            lNumGhostNodes += lCoarseNodesPerDir[i];
            // Add corners if three adjoining faces need ghost nodes,
            // but add them only once! Hence when i == 0.
            if (ghostInterface[2 * i] && (i == 0)) {
              lNumGhostNodes += 1;
            }
            if (ghostInterface[2 * i + 1] && (i == 0)) {
              lNumGhostNodes += 1;
            }
          }
        }
      }
      tmp = 0;
    }
  }  // end of scope for tmp and complementaryIndices

  const LO lNumGhostedNodes = glCoarseNodesPerDir[0] * glCoarseNodesPerDir[1] * glCoarseNodesPerDir[2];
  ghostedCoarseNodes->PIDs.resize(lNumGhostedNodes);
  ghostedCoarseNodes->LIDs.resize(lNumGhostedNodes);
  ghostedCoarseNodes->GIDs.resize(lNumGhostedNodes);
  ghostedCoarseNodes->coarseGIDs.resize(lNumGhostedNodes);
  ghostedCoarseNodes->colInds.resize(lNumGhostedNodes);

  // We loop over all ghosted coarse nodes by increasing global lexicographic order
  Array<LO> coarseNodeCoarseIndices(3), coarseNodeFineIndices(3), ijk(3);
  LO currentIndex = -1;
  for (ijk[2] = 0; ijk[2] < glCoarseNodesPerDir[2]; ++ijk[2]) {
    for (ijk[1] = 0; ijk[1] < glCoarseNodesPerDir[1]; ++ijk[1]) {
      for (ijk[0] = 0; ijk[0] < glCoarseNodesPerDir[0]; ++ijk[0]) {
        currentIndex               = ijk[2] * glCoarseNodesPerDir[1] * glCoarseNodesPerDir[0] + ijk[1] * glCoarseNodesPerDir[0] + ijk[0];
        coarseNodeCoarseIndices[0] = startGhostedCoarseNode[0] + ijk[0];
        coarseNodeFineIndices[0]   = coarseNodeCoarseIndices[0] * coarseRate[0];
        if (coarseNodeFineIndices[0] > gFineNodesPerDir[0] - 1) {
          coarseNodeFineIndices[0] = gFineNodesPerDir[0] - 1;
        }
        coarseNodeCoarseIndices[1] = startGhostedCoarseNode[1] + ijk[1];
        coarseNodeFineIndices[1]   = coarseNodeCoarseIndices[1] * coarseRate[1];
        if (coarseNodeFineIndices[1] > gFineNodesPerDir[1] - 1) {
          coarseNodeFineIndices[1] = gFineNodesPerDir[1] - 1;
        }
        coarseNodeCoarseIndices[2] = startGhostedCoarseNode[2] + ijk[2];
        coarseNodeFineIndices[2]   = coarseNodeCoarseIndices[2] * coarseRate[2];
        if (coarseNodeFineIndices[2] > gFineNodesPerDir[2] - 1) {
          coarseNodeFineIndices[2] = gFineNodesPerDir[2] - 1;
        }
        GO myGID = 0, myCoarseGID = -1;
        GO factor[3] = {};
        factor[2]    = gFineNodesPerDir[1] * gFineNodesPerDir[0];
        factor[1]    = gFineNodesPerDir[0];
        factor[0]    = 1;
        for (int dim = 0; dim < 3; ++dim) {
          if (dim < numDimensions) {
            if (gIndices[dim] - myOffset[dim] + ijk[dim] * coarseRate[dim] < gFineNodesPerDir[dim] - 1) {
              myGID += (gIndices[dim] - myOffset[dim] + ijk[dim] * coarseRate[dim]) * factor[dim];
            } else {
              myGID += (gIndices[dim] - myOffset[dim] + (ijk[dim] - 1) * coarseRate[dim] + endRate[dim]) * factor[dim];
            }
          }
        }
        myCoarseGID                                  = coarseNodeCoarseIndices[0] + coarseNodeCoarseIndices[1] * gCoarseNodesPerDir[0] + coarseNodeCoarseIndices[2] * gCoarseNodesPerDir[1] * gCoarseNodesPerDir[0];
        ghostedCoarseNodes->GIDs[currentIndex]       = myGID;
        ghostedCoarseNodes->coarseGIDs[currentIndex] = myCoarseGID;
      }
    }
  }
  coordinatesMap->getRemoteIndexList(ghostedCoarseNodes->GIDs(),
                                     ghostedCoarseNodes->PIDs(),
                                     ghostedCoarseNodes->LIDs());

  //                              Phase 2                               //
  // ------------------------------------------------------------------ //
  // Build a list of GIDs to import the required ghost nodes.           //
  // The ordering of the ghosts nodes will be as natural as possible,   //
  // i.e. it should follow the GID ordering of the mesh.                //
  // ------------------------------------------------------------------ //

  // Saddly we have to more or less redo what was just done to figure out the value of
  // lNumGhostNodes, there might be some optimization possibility here...
  ghostGIDs.resize(lNumGhostNodes);
  LO countGhosts = 0;
  // Get the GID of the first point on the current partition.
  GO startingGID = minGlobalIndex;
  Array<GO> startingIndices(3);
  // We still want ghost nodes even if have with a 0 offset,
  // except when we are on a boundary
  if (ghostInterface[4] && (myOffset[2] == 0)) {
    if (gIndices[2] + coarseRate[2] > gFineNodesPerDir[2]) {
      startingGID -= endRate[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
    } else {
      startingGID -= coarseRate[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
    }
  } else {
    startingGID -= myOffset[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
  }
  if (ghostInterface[2] && (myOffset[1] == 0)) {
    if (gIndices[1] + coarseRate[1] > gFineNodesPerDir[1]) {
      startingGID -= endRate[1] * gFineNodesPerDir[0];
    } else {
      startingGID -= coarseRate[1] * gFineNodesPerDir[0];
    }
  } else {
    startingGID -= myOffset[1] * gFineNodesPerDir[0];
  }
  if (ghostInterface[0] && (myOffset[0] == 0)) {
    if (gIndices[0] + coarseRate[0] > gFineNodesPerDir[0]) {
      startingGID -= endRate[0];
    } else {
      startingGID -= coarseRate[0];
    }
  } else {
    startingGID -= myOffset[0];
  }

  {  // scope for tmp
    GO tmp;
    startingIndices[2] = startingGID / (gFineNodesPerDir[1] * gFineNodesPerDir[0]);
    tmp                = startingGID % (gFineNodesPerDir[1] * gFineNodesPerDir[0]);
    startingIndices[1] = tmp / gFineNodesPerDir[0];
    startingIndices[0] = tmp % gFineNodesPerDir[0];
  }

  GO ghostOffset[3] = {0, 0, 0};
  LO lengthZero     = lCoarseNodesPerDir[0];
  LO lengthOne      = lCoarseNodesPerDir[1];
  LO lengthTwo      = lCoarseNodesPerDir[2];
  if (ghostInterface[0]) {
    ++lengthZero;
  }
  if (ghostInterface[1]) {
    ++lengthZero;
  }
  if (ghostInterface[2]) {
    ++lengthOne;
  }
  if (ghostInterface[3]) {
    ++lengthOne;
  }
  if (ghostInterface[4]) {
    ++lengthTwo;
  }
  if (ghostInterface[5]) {
    ++lengthTwo;
  }

  // First check the bottom face as it will have the lowest GIDs
  if (ghostInterface[4]) {
    ghostOffset[2] = startingGID;
    for (LO j = 0; j < lengthOne; ++j) {
      if ((j == lengthOne - 1) && (startingIndices[1] + j * coarseRate[1] + 1 > gFineNodesPerDir[1])) {
        ghostOffset[1] = ((j - 1) * coarseRate[1] + endRate[1]) * gFineNodesPerDir[0];
      } else {
        ghostOffset[1] = j * coarseRate[1] * gFineNodesPerDir[0];
      }
      for (LO k = 0; k < lengthZero; ++k) {
        if ((k == lengthZero - 1) && (startingIndices[0] + k * coarseRate[0] + 1 > gFineNodesPerDir[0])) {
          ghostOffset[0] = (k - 1) * coarseRate[0] + endRate[0];
        } else {
          ghostOffset[0] = k * coarseRate[0];
        }
        // If the partition includes a changed rate at the edge, ghost nodes need to be picked
        // carefully.
        // This if statement is repeated each time a ghost node is selected.
        ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
        ++countGhosts;
      }
    }
  }

  // Sweep over the lCoarseNodesPerDir[2] coarse layers in direction 2 and gather necessary ghost
  // nodes located on these layers.
  for (LO i = 0; i < lengthTwo; ++i) {
    // Exclude the cases where ghost nodes exists on the faces in directions 2, these faces are
    // swept seperatly for efficiency.
    if (!((i == lengthTwo - 1) && ghostInterface[5]) && !((i == 0) && ghostInterface[4])) {
      // Set the ghostOffset in direction 2 taking into account a possible endRate different
      // from the regular coarseRate.
      if ((i == lengthTwo - 1) && !ghostInterface[5]) {
        ghostOffset[2] = startingGID + ((i - 1) * coarseRate[2] + endRate[2]) * gFineNodesPerDir[1] * gFineNodesPerDir[0];
      } else {
        ghostOffset[2] = startingGID + i * coarseRate[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
      }
      for (LO j = 0; j < lengthOne; ++j) {
        if ((j == 0) && ghostInterface[2]) {
          for (LO k = 0; k < lengthZero; ++k) {
            if ((k == lengthZero - 1) && (startingIndices[0] + k * coarseRate[0] + 1 > gFineNodesPerDir[0])) {
              if (k == 0) {
                ghostOffset[0] = endRate[0];
              } else {
                ghostOffset[0] = (k - 1) * coarseRate[0] + endRate[0];
              }
            } else {
              ghostOffset[0] = k * coarseRate[0];
            }
            // In this case j == 0 so ghostOffset[1] is zero and can be ignored
            ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[0];
            ++countGhosts;
          }
        } else if ((j == lengthOne - 1) && ghostInterface[3]) {
          // Set the ghostOffset in direction 1 taking into account a possible endRate different
          // from the regular coarseRate.
          if ((j == lengthOne - 1) && (startingIndices[1] + j * coarseRate[1] + 1 > gFineNodesPerDir[1])) {
            ghostOffset[1] = ((j - 1) * coarseRate[1] + endRate[1]) * gFineNodesPerDir[0];
          } else {
            ghostOffset[1] = j * coarseRate[1] * gFineNodesPerDir[0];
          }
          for (LO k = 0; k < lengthZero; ++k) {
            if ((k == lengthZero - 1) && (startingIndices[0] + k * coarseRate[0] + 1 > gFineNodesPerDir[0])) {
              ghostOffset[0] = (k - 1) * coarseRate[0] + endRate[0];
            } else {
              ghostOffset[0] = k * coarseRate[0];
            }
            ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
            ++countGhosts;
          }
        } else {
          // Set ghostOffset[1] for side faces sweep
          if ((j == lengthOne - 1) && (startingIndices[1] + j * coarseRate[1] + 1 > gFineNodesPerDir[1])) {
            ghostOffset[1] = ((j - 1) * coarseRate[1] + endRate[1]) * gFineNodesPerDir[0];
          } else {
            ghostOffset[1] = j * coarseRate[1] * gFineNodesPerDir[0];
          }

          // Set ghostOffset[0], ghostGIDs and countGhosts
          if (ghostInterface[0]) {  // In that case ghostOffset[0]==0, so we can ignore it
            ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1];
            ++countGhosts;
          }
          if (ghostInterface[1]) {  // Grab ghost point at the end of direction 0.
            if ((startingIndices[0] + (lengthZero - 1) * coarseRate[0]) > gFineNodesPerDir[0] - 1) {
              if (lengthZero > 2) {
                ghostOffset[0] = (lengthZero - 2) * coarseRate[0] + endRate[0];
              } else {
                ghostOffset[0] = endRate[0];
              }
            } else {
              ghostOffset[0] = (lengthZero - 1) * coarseRate[0];
            }
            ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
            ++countGhosts;
          }
        }
      }
    }
  }

  // Finally check the top face
  if (ghostInterface[5]) {
    if (startingIndices[2] + (lengthTwo - 1) * coarseRate[2] + 1 > gFineNodesPerDir[2]) {
      ghostOffset[2] = startingGID + ((lengthTwo - 2) * coarseRate[2] + endRate[2]) * gFineNodesPerDir[1] * gFineNodesPerDir[0];
    } else {
      ghostOffset[2] = startingGID + (lengthTwo - 1) * coarseRate[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
    }
    for (LO j = 0; j < lengthOne; ++j) {
      if ((j == lengthOne - 1) && (startingIndices[1] + j * coarseRate[1] + 1 > gFineNodesPerDir[1])) {
        ghostOffset[1] = ((j - 1) * coarseRate[1] + endRate[1]) * gFineNodesPerDir[0];
      } else {
        ghostOffset[1] = j * coarseRate[1] * gFineNodesPerDir[0];
      }
      for (LO k = 0; k < lengthZero; ++k) {
        if ((k == lengthZero - 1) && (startingIndices[0] + k * coarseRate[0] + 1 > gFineNodesPerDir[0])) {
          ghostOffset[0] = (k - 1) * coarseRate[0] + endRate[0];
        } else {
          ghostOffset[0] = k * coarseRate[0];
        }
        ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
        ++countGhosts;
      }
    }
  }

  //                              Phase 3                               //
  // ------------------------------------------------------------------ //
  // Final phase of this function, lists are being built to create the  //
  // column and domain maps of the projection as well as the map and    //
  // arrays for the coarse coordinates multivector.                     //
  // ------------------------------------------------------------------ //

  Array<GO> gCoarseNodesGIDs(lNumCoarseNodes);
  LO currentNode, offset2, offset1, offset0;
  // Find the GIDs of the coarse nodes on the partition.
  for (LO ind2 = 0; ind2 < lCoarseNodesPerDir[2]; ++ind2) {
    if (myOffset[2] == 0) {
      offset2 = startingIndices[2] + myOffset[2];
    } else {
      if (startingIndices[2] + endRate[2] == gFineNodesPerDir[2] - 1) {
        offset2 = startingIndices[2] + endRate[2];
      } else {
        offset2 = startingIndices[2] + coarseRate[2];
      }
    }
    if (offset2 + ind2 * coarseRate[2] > gFineNodesPerDir[2] - 1) {
      offset2 += (ind2 - 1) * coarseRate[2] + endRate[2];
    } else {
      offset2 += ind2 * coarseRate[2];
    }
    offset2 = offset2 * gFineNodesPerDir[1] * gFineNodesPerDir[0];

    for (LO ind1 = 0; ind1 < lCoarseNodesPerDir[1]; ++ind1) {
      if (myOffset[1] == 0) {
        offset1 = startingIndices[1] + myOffset[1];
      } else {
        if (startingIndices[1] + endRate[1] == gFineNodesPerDir[1] - 1) {
          offset1 = startingIndices[1] + endRate[1];
        } else {
          offset1 = startingIndices[1] + coarseRate[1];
        }
      }
      if (offset1 + ind1 * coarseRate[1] > gFineNodesPerDir[1] - 1) {
        offset1 += (ind1 - 1) * coarseRate[1] + endRate[1];
      } else {
        offset1 += ind1 * coarseRate[1];
      }
      offset1 = offset1 * gFineNodesPerDir[0];
      for (LO ind0 = 0; ind0 < lCoarseNodesPerDir[0]; ++ind0) {
        offset0 = startingIndices[0];
        if (myOffset[0] == 0) {
          offset0 += ind0 * coarseRate[0];
        } else {
          offset0 += (ind0 + 1) * coarseRate[0];
        }
        if (offset0 > gFineNodesPerDir[0] - 1) {
          offset0 += endRate[0] - coarseRate[0];
        }

        currentNode                   = ind2 * lCoarseNodesPerDir[1] * lCoarseNodesPerDir[0] + ind1 * lCoarseNodesPerDir[0] + ind0;
        gCoarseNodesGIDs[currentNode] = offset2 + offset1 + offset0;
      }
    }
  }

  // Actual loop over all the coarse/ghost nodes to find their index on the coarse mesh
  // and the corresponding dofs that will need to be added to colMapP.
  colGIDs.resize(BlkSize * (lNumCoarseNodes + lNumGhostNodes));
  coarseNodesGIDs.resize(lNumCoarseNodes);
  for (LO i = 0; i < numDimensions; ++i) {
    coarseNodes[i].resize(lNumCoarseNodes);
  }
  GO fineNodesPerCoarseSlab    = coarseRate[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
  GO fineNodesEndCoarseSlab    = endRate[2] * gFineNodesPerDir[1] * gFineNodesPerDir[0];
  GO fineNodesPerCoarsePlane   = coarseRate[1] * gFineNodesPerDir[0];
  GO fineNodesEndCoarsePlane   = endRate[1] * gFineNodesPerDir[0];
  GO coarseNodesPerCoarseLayer = gCoarseNodesPerDir[1] * gCoarseNodesPerDir[0];
  GO gCoarseNodeOnCoarseGridGID;
  LO gInd[3], lCol;
  Array<int> ghostPIDs(lNumGhostNodes);
  Array<LO> ghostLIDs(lNumGhostNodes);
  Array<LO> ghostPermut(lNumGhostNodes);
  for (LO k = 0; k < lNumGhostNodes; ++k) {
    ghostPermut[k] = k;
  }
  coordinatesMap->getRemoteIndexList(ghostGIDs, ghostPIDs, ghostLIDs);
  sh_sort_permute(ghostPIDs.begin(), ghostPIDs.end(), ghostPermut.begin(), ghostPermut.end());

  {  // scope for tmpInds, tmpVars and tmp.
    GO tmpInds[3], tmpVars[2];
    LO tmp;
    // Loop over the coarse nodes of the partition and add them to colGIDs
    // that will be used to construct the column and domain maps of P as well
    // as to construct the coarse coordinates map.
    // for(LO col = 0; col < lNumCoarseNodes; ++col) { // This should most likely be replaced
    // by loops of lCoarseNodesPerDir[] to simplify arithmetics
    LO col = 0;
    LO firstCoarseNodeInds[3], currentCoarseNode;
    for (LO dim = 0; dim < 3; ++dim) {
      if (myOffset[dim] == 0) {
        firstCoarseNodeInds[dim] = 0;
      } else {
        firstCoarseNodeInds[dim] = coarseRate[dim] - myOffset[dim];
      }
    }
    Array<ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > fineNodes(numDimensions);
    for (LO dim = 0; dim < numDimensions; ++dim) {
      fineNodes[dim] = coordinates->getData(dim);
    }
    for (LO k = 0; k < lCoarseNodesPerDir[2]; ++k) {
      for (LO j = 0; j < lCoarseNodesPerDir[1]; ++j) {
        for (LO i = 0; i < lCoarseNodesPerDir[0]; ++i) {
          col = k * lCoarseNodesPerDir[1] * lCoarseNodesPerDir[0] + j * lCoarseNodesPerDir[0] + i;

          // Check for endRate
          currentCoarseNode = 0;
          if (firstCoarseNodeInds[0] + i * coarseRate[0] > lFineNodesPerDir[0] - 1) {
            currentCoarseNode += firstCoarseNodeInds[0] + (i - 1) * coarseRate[0] + endRate[0];
          } else {
            currentCoarseNode += firstCoarseNodeInds[0] + i * coarseRate[0];
          }
          if (firstCoarseNodeInds[1] + j * coarseRate[1] > lFineNodesPerDir[1] - 1) {
            currentCoarseNode += (firstCoarseNodeInds[1] + (j - 1) * coarseRate[1] + endRate[1]) * lFineNodesPerDir[0];
          } else {
            currentCoarseNode += (firstCoarseNodeInds[1] + j * coarseRate[1]) * lFineNodesPerDir[0];
          }
          if (firstCoarseNodeInds[2] + k * coarseRate[2] > lFineNodesPerDir[2] - 1) {
            currentCoarseNode += (firstCoarseNodeInds[2] + (k - 1) * coarseRate[2] + endRate[2]) * lFineNodesPerDir[1] * lFineNodesPerDir[0];
          } else {
            currentCoarseNode += (firstCoarseNodeInds[2] + k * coarseRate[2]) * lFineNodesPerDir[1] * lFineNodesPerDir[0];
          }
          // Load coordinates
          for (LO dim = 0; dim < numDimensions; ++dim) {
            coarseNodes[dim][col] = fineNodes[dim][currentCoarseNode];
          }

          if ((endRate[2] != coarseRate[2]) && (gCoarseNodesGIDs[col] > (gCoarseNodesPerDir[2] - 2) * fineNodesPerCoarseSlab + fineNodesEndCoarseSlab - 1)) {
            tmpInds[2] = gCoarseNodesGIDs[col] / fineNodesPerCoarseSlab + 1;
            tmpVars[0] = gCoarseNodesGIDs[col] - (tmpInds[2] - 1) * fineNodesPerCoarseSlab - fineNodesEndCoarseSlab;
          } else {
            tmpInds[2] = gCoarseNodesGIDs[col] / fineNodesPerCoarseSlab;
            tmpVars[0] = gCoarseNodesGIDs[col] % fineNodesPerCoarseSlab;
          }
          if ((endRate[1] != coarseRate[1]) && (tmpVars[0] > (gCoarseNodesPerDir[1] - 2) * fineNodesPerCoarsePlane + fineNodesEndCoarsePlane - 1)) {
            tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane + 1;
            tmpVars[1] = tmpVars[0] - (tmpInds[1] - 1) * fineNodesPerCoarsePlane - fineNodesEndCoarsePlane;
          } else {
            tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane;
            tmpVars[1] = tmpVars[0] % fineNodesPerCoarsePlane;
          }
          if (tmpVars[1] == gFineNodesPerDir[0] - 1) {
            tmpInds[0] = gCoarseNodesPerDir[0] - 1;
          } else {
            tmpInds[0] = tmpVars[1] / coarseRate[0];
          }
          gInd[2]                    = col / (lCoarseNodesPerDir[1] * lCoarseNodesPerDir[0]);
          tmp                        = col % (lCoarseNodesPerDir[1] * lCoarseNodesPerDir[0]);
          gInd[1]                    = tmp / lCoarseNodesPerDir[0];
          gInd[0]                    = tmp % lCoarseNodesPerDir[0];
          lCol                       = gInd[2] * (lCoarseNodesPerDir[1] * lCoarseNodesPerDir[0]) + gInd[1] * lCoarseNodesPerDir[0] + gInd[0];
          gCoarseNodeOnCoarseGridGID = tmpInds[2] * coarseNodesPerCoarseLayer + tmpInds[1] * gCoarseNodesPerDir[0] + tmpInds[0];
          coarseNodesGIDs[lCol]      = gCoarseNodeOnCoarseGridGID;
          for (LO dof = 0; dof < BlkSize; ++dof) {
            colGIDs[BlkSize * lCol + dof] = BlkSize * gCoarseNodeOnCoarseGridGID + dof;
          }
        }
      }
    }
    // Now loop over the ghost nodes of the partition to add them to colGIDs
    // since they will need to be included in the column map of P
    for (col = lNumCoarseNodes; col < lNumCoarseNodes + lNumGhostNodes; ++col) {
      if ((endRate[2] != coarseRate[2]) && (ghostGIDs[ghostPermut[col - lNumCoarseNodes]] > (gCoarseNodesPerDir[2] - 2) * fineNodesPerCoarseSlab + fineNodesEndCoarseSlab - 1)) {
        tmpInds[2] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] / fineNodesPerCoarseSlab + 1;
        tmpVars[0] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] - (tmpInds[2] - 1) * fineNodesPerCoarseSlab - fineNodesEndCoarseSlab;
      } else {
        tmpInds[2] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] / fineNodesPerCoarseSlab;
        tmpVars[0] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] % fineNodesPerCoarseSlab;
      }
      if ((endRate[1] != coarseRate[1]) && (tmpVars[0] > (gCoarseNodesPerDir[1] - 2) * fineNodesPerCoarsePlane + fineNodesEndCoarsePlane - 1)) {
        tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane + 1;
        tmpVars[1] = tmpVars[0] - (tmpInds[1] - 1) * fineNodesPerCoarsePlane - fineNodesEndCoarsePlane;
      } else {
        tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane;
        tmpVars[1] = tmpVars[0] % fineNodesPerCoarsePlane;
      }
      if (tmpVars[1] == gFineNodesPerDir[0] - 1) {
        tmpInds[0] = gCoarseNodesPerDir[0] - 1;
      } else {
        tmpInds[0] = tmpVars[1] / coarseRate[0];
      }
      gCoarseNodeOnCoarseGridGID = tmpInds[2] * coarseNodesPerCoarseLayer + tmpInds[1] * gCoarseNodesPerDir[0] + tmpInds[0];
      for (LO dof = 0; dof < BlkSize; ++dof) {
        colGIDs[BlkSize * col + dof] = BlkSize * gCoarseNodeOnCoarseGridGID + dof;
      }
    }
  }  // End of scope for tmpInds, tmpVars and tmp

}  // GetGeometricData()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ComputeLocalEntries(const RCP<const Matrix>& Aghost, const Array<LO> coarseRate,
                        const Array<LO> /* endRate */, const LO BlkSize, const Array<LO> elemInds,
                        const Array<LO> /* lCoarseElementsPerDir */, const LO numDimensions,
                        const Array<LO> lFineNodesPerDir, const Array<GO> /* gFineNodesPerDir */,
                        const Array<GO> /* gIndices */, const Array<LO> /* lCoarseNodesPerDir */,
                        const Array<bool> ghostInterface, const Array<int> elementFlags,
                        const std::string stencilType, const std::string /* blockStrategy */,
                        const Array<LO> elementNodesPerDir, const LO numNodesInElement,
                        const Array<GO> /* colGIDs */,
                        Teuchos::SerialDenseMatrix<LO, SC>& Pi, Teuchos::SerialDenseMatrix<LO, SC>& Pf,
                        Teuchos::SerialDenseMatrix<LO, SC>& Pe, Array<LO>& dofType,
                        Array<LO>& lDofInd) const {
  // First extract data from Aghost and move it to the corresponding dense matrix
  // This step requires to compute the number of nodes (resp dofs) in the current
  // coarse element, then allocate storage for local dense matrices, finally loop
  // over the dofs and extract the corresponding row in Aghost before dispactching
  // its data to the dense matrices.
  // At the same time we want to compute a mapping from the element indexing to
  // group indexing. The groups are the following: corner, edge, face and interior
  // nodes. We uses these groups to operate on the dense matrices but need to
  // store the nodes in their original order after groupd operations are completed.
  LO countInterior = 0, countFace = 0, countEdge = 0, countCorner = 0;
  Array<LO> collapseDir(numNodesInElement * BlkSize);
  for (LO ke = 0; ke < elementNodesPerDir[2]; ++ke) {
    for (LO je = 0; je < elementNodesPerDir[1]; ++je) {
      for (LO ie = 0; ie < elementNodesPerDir[0]; ++ie) {
        // Check for corner node
        if ((ke == 0 || ke == elementNodesPerDir[2] - 1) && (je == 0 || je == elementNodesPerDir[1] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1)) {
          for (LO dof = 0; dof < BlkSize; ++dof) {
            dofType[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 0;
            lDofInd[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = BlkSize * countCorner + dof;
          }
          ++countCorner;

          // Check for edge node
        } else if (((ke == 0 || ke == elementNodesPerDir[2] - 1) && (je == 0 || je == elementNodesPerDir[1] - 1)) || ((ke == 0 || ke == elementNodesPerDir[2] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1)) || ((je == 0 || je == elementNodesPerDir[1] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1))) {
          for (LO dof = 0; dof < BlkSize; ++dof) {
            dofType[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 1;
            lDofInd[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = BlkSize * countEdge + dof;
            if ((ke == 0 || ke == elementNodesPerDir[2] - 1) && (je == 0 || je == elementNodesPerDir[1] - 1)) {
              collapseDir[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 0;
            } else if ((ke == 0 || ke == elementNodesPerDir[2] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1)) {
              collapseDir[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 1;
            } else if ((je == 0 || je == elementNodesPerDir[1] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1)) {
              collapseDir[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 2;
            }
          }
          ++countEdge;

          // Check for face node
        } else if ((ke == 0 || ke == elementNodesPerDir[2] - 1) || (je == 0 || je == elementNodesPerDir[1] - 1) || (ie == 0 || ie == elementNodesPerDir[0] - 1)) {
          for (LO dof = 0; dof < BlkSize; ++dof) {
            dofType[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 2;
            lDofInd[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = BlkSize * countFace + dof;
            if (ke == 0 || ke == elementNodesPerDir[2] - 1) {
              collapseDir[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 2;
            } else if (je == 0 || je == elementNodesPerDir[1] - 1) {
              collapseDir[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 1;
            } else if (ie == 0 || ie == elementNodesPerDir[0] - 1) {
              collapseDir[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 0;
            }
          }
          ++countFace;

          // Otherwise it is an interior node
        } else {
          for (LO dof = 0; dof < BlkSize; ++dof) {
            dofType[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = 3;
            lDofInd[BlkSize * (ke * elementNodesPerDir[1] * elementNodesPerDir[0] + je * elementNodesPerDir[0] + ie) + dof] = BlkSize * countInterior + dof;
          }
          ++countInterior;
        }
      }
    }
  }

  LO numInteriorNodes = 0, numFaceNodes = 0, numEdgeNodes = 0, numCornerNodes = 8;
  numInteriorNodes = (elementNodesPerDir[0] - 2) * (elementNodesPerDir[1] - 2) * (elementNodesPerDir[2] - 2);
  numFaceNodes     = 2 * (elementNodesPerDir[0] - 2) * (elementNodesPerDir[1] - 2) + 2 * (elementNodesPerDir[0] - 2) * (elementNodesPerDir[2] - 2) + 2 * (elementNodesPerDir[1] - 2) * (elementNodesPerDir[2] - 2);
  numEdgeNodes     = 4 * (elementNodesPerDir[0] - 2) + 4 * (elementNodesPerDir[1] - 2) + 4 * (elementNodesPerDir[2] - 2);
  // Diagonal blocks
  Teuchos::SerialDenseMatrix<LO, SC> Aii(BlkSize * numInteriorNodes, BlkSize * numInteriorNodes);
  Teuchos::SerialDenseMatrix<LO, SC> Aff(BlkSize * numFaceNodes, BlkSize * numFaceNodes);
  Teuchos::SerialDenseMatrix<LO, SC> Aee(BlkSize * numEdgeNodes, BlkSize * numEdgeNodes);
  // Upper triangular blocks
  Teuchos::SerialDenseMatrix<LO, SC> Aif(BlkSize * numInteriorNodes, BlkSize * numFaceNodes);
  Teuchos::SerialDenseMatrix<LO, SC> Aie(BlkSize * numInteriorNodes, BlkSize * numEdgeNodes);
  Teuchos::SerialDenseMatrix<LO, SC> Afe(BlkSize * numFaceNodes, BlkSize * numEdgeNodes);
  // Coarse nodes "right hand sides" and local prolongators
  Teuchos::SerialDenseMatrix<LO, SC> Aic(BlkSize * numInteriorNodes, BlkSize * numCornerNodes);
  Teuchos::SerialDenseMatrix<LO, SC> Afc(BlkSize * numFaceNodes, BlkSize * numCornerNodes);
  Teuchos::SerialDenseMatrix<LO, SC> Aec(BlkSize * numEdgeNodes, BlkSize * numCornerNodes);
  Pi.shape(BlkSize * numInteriorNodes, BlkSize * numCornerNodes);
  Pf.shape(BlkSize * numFaceNodes, BlkSize * numCornerNodes);
  Pe.shape(BlkSize * numEdgeNodes, BlkSize * numCornerNodes);

  ArrayView<const LO> rowIndices;
  ArrayView<const SC> rowValues;
  LO idof, iInd, jInd;
  int iType = 0, jType = 0;
  int orientation      = -1;
  int collapseFlags[3] = {};
  Array<SC> stencil((std::pow(3, numDimensions)) * BlkSize);
  // LBV Note: we could skip the extraction of rows corresponding to coarse nodes
  //           we might want to fuse the first three loops and do integer arithmetic
  //           to have more optimization during compilation...
  for (LO ke = 0; ke < elementNodesPerDir[2]; ++ke) {
    for (LO je = 0; je < elementNodesPerDir[1]; ++je) {
      for (LO ie = 0; ie < elementNodesPerDir[0]; ++ie) {
        collapseFlags[0] = 0;
        collapseFlags[1] = 0;
        collapseFlags[2] = 0;
        if ((elementFlags[0] == 1 || elementFlags[0] == 3) && ie == 0) {
          collapseFlags[0] += 1;
        }
        if ((elementFlags[0] == 2 || elementFlags[0] == 3) && ie == elementNodesPerDir[0] - 1) {
          collapseFlags[0] += 2;
        }
        if ((elementFlags[1] == 1 || elementFlags[1] == 3) && je == 0) {
          collapseFlags[1] += 1;
        }
        if ((elementFlags[1] == 2 || elementFlags[1] == 3) && je == elementNodesPerDir[1] - 1) {
          collapseFlags[1] += 2;
        }
        if ((elementFlags[2] == 1 || elementFlags[2] == 3) && ke == 0) {
          collapseFlags[2] += 1;
        }
        if ((elementFlags[2] == 2 || elementFlags[2] == 3) && ke == elementNodesPerDir[2] - 1) {
          collapseFlags[2] += 2;
        }

        // Based on (ie, je, ke) we detect the type of node we are dealing with
        GetNodeInfo(ie, je, ke, elementNodesPerDir, &iType, iInd, &orientation);
        for (LO dof0 = 0; dof0 < BlkSize; ++dof0) {
          // Compute the dof index for each dof at point (ie, je, ke)
          idof = ((elemInds[2] * coarseRate[2] + ke) * lFineNodesPerDir[1] * lFineNodesPerDir[0] + (elemInds[1] * coarseRate[1] + je) * lFineNodesPerDir[0] + elemInds[0] * coarseRate[0] + ie) * BlkSize + dof0;
          Aghost->getLocalRowView(idof, rowIndices, rowValues);
          FormatStencil(BlkSize, ghostInterface, ie, je, ke, rowValues,
                        elementNodesPerDir, collapseFlags, stencilType, stencil);
          LO io, jo, ko;
          if (iType == 3) {  // interior node, no stencil collapse needed
            for (LO interactingNode = 0; interactingNode < 27; ++interactingNode) {
              // Find (if, jf, kf) the indices associated with the interacting node
              ko = ke + (interactingNode / 9 - 1);
              {
                LO tmp = interactingNode % 9;
                jo     = je + (tmp / 3 - 1);
                io     = ie + (tmp % 3 - 1);
                int dummy;
                GetNodeInfo(io, jo, ko, elementNodesPerDir, &jType, jInd, &dummy);
              }
              if ((io > -1 && io < elementNodesPerDir[0]) &&
                  (jo > -1 && jo < elementNodesPerDir[1]) &&
                  (ko > -1 && ko < elementNodesPerDir[2])) {
                for (LO dof1 = 0; dof1 < BlkSize; ++dof1) {
                  if (jType == 3) {
                    Aii(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        stencil[interactingNode * BlkSize + dof1];
                  } else if (jType == 2) {
                    Aif(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        stencil[interactingNode * BlkSize + dof1];
                  } else if (jType == 1) {
                    Aie(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        stencil[interactingNode * BlkSize + dof1];
                  } else if (jType == 0) {
                    Aic(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        -stencil[interactingNode * BlkSize + dof1];
                  }
                }
              }
            }
          } else if (iType == 2) {  // Face node, collapse stencil along face normal (*orientation)
            CollapseStencil(2, orientation, collapseFlags, stencil);
            for (LO interactingNode = 0; interactingNode < 27; ++interactingNode) {
              // Find (if, jf, kf) the indices associated with the interacting node
              ko = ke + (interactingNode / 9 - 1);
              {
                LO tmp = interactingNode % 9;
                jo     = je + (tmp / 3 - 1);
                io     = ie + (tmp % 3 - 1);
                int dummy;
                GetNodeInfo(io, jo, ko, elementNodesPerDir, &jType, jInd, &dummy);
              }
              if ((io > -1 && io < elementNodesPerDir[0]) &&
                  (jo > -1 && jo < elementNodesPerDir[1]) &&
                  (ko > -1 && ko < elementNodesPerDir[2])) {
                for (LO dof1 = 0; dof1 < BlkSize; ++dof1) {
                  if (jType == 2) {
                    Aff(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        stencil[interactingNode * BlkSize + dof1];
                  } else if (jType == 1) {
                    Afe(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        stencil[interactingNode * BlkSize + dof1];
                  } else if (jType == 0) {
                    Afc(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        -stencil[interactingNode * BlkSize + dof1];
                  }
                }
              }
            }
          } else if (iType == 1) {  // Edge node, collapse stencil perpendicular to edge direction
            CollapseStencil(1, orientation, collapseFlags, stencil);
            for (LO interactingNode = 0; interactingNode < 27; ++interactingNode) {
              // Find (if, jf, kf) the indices associated with the interacting node
              ko = ke + (interactingNode / 9 - 1);
              {
                LO tmp = interactingNode % 9;
                jo     = je + (tmp / 3 - 1);
                io     = ie + (tmp % 3 - 1);
                int dummy;
                GetNodeInfo(io, jo, ko, elementNodesPerDir, &jType, jInd, &dummy);
              }
              if ((io > -1 && io < elementNodesPerDir[0]) &&
                  (jo > -1 && jo < elementNodesPerDir[1]) &&
                  (ko > -1 && ko < elementNodesPerDir[2])) {
                for (LO dof1 = 0; dof1 < BlkSize; ++dof1) {
                  if (jType == 1) {
                    Aee(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        stencil[interactingNode * BlkSize + dof1];
                  } else if (jType == 0) {
                    Aec(iInd * BlkSize + dof0, jInd * BlkSize + dof1) =
                        -stencil[interactingNode * BlkSize + dof1];
                  }
                }  // Check if interacting node is in element or not
              }    // Pc is the identity so no need to treat iType == 0
            }      // Loop over interacting nodes within a row
          }        // Check for current node type
        }          // Loop over degrees of freedom associated to a given node
      }            // Loop over i
    }              // Loop over j
  }                // Loop over k

  // Compute the projection operators for edge and interior nodes
  //
  //         [P_i] = [A_{ii} & A_{if} & A_{ie}]^{-1} [A_{ic}]
  //         [P_f] = [  0    & A_{ff} & A_{fe}]      [A_{fc}]
  //         [P_e] = [  0    &   0    & A_{ee}]      [A_{ec}]
  //         [P_c] =  I_c
  //
  {  // Compute P_e
    // We need to solve for P_e in: A_{ee}*P_e = A_{fc}
    // So we simple do P_e = A_{ee}^{-1)*A_{ec}
    Teuchos::SerialDenseSolver<LO, SC> problem;
    problem.setMatrix(Teuchos::rcp(&Aee, false));
    problem.setVectors(Teuchos::rcp(&Pe, false), Teuchos::rcp(&Aec, false));
    problem.factorWithEquilibration(true);
    problem.solve();
    problem.unequilibrateLHS();
  }

  {  // Compute P_f
    // We need to solve for P_f in: A_{ff}*P_f + A_{fe}*P_e = A_{fc}
    // Step one: A_{fc} = 1.0*A_{fc} + (-1.0)*A_{fe}*P_e
    Afc.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, Afe, Pe, 1.0);
    // Step two: P_f = A_{ff}^{-1}*A_{fc}
    Teuchos::SerialDenseSolver<LO, SC> problem;
    problem.setMatrix(Teuchos::rcp(&Aff, false));
    problem.setVectors(Teuchos::rcp(&Pf, false), Teuchos::rcp(&Afc, false));
    problem.factorWithEquilibration(true);
    problem.solve();
    problem.unequilibrateLHS();
  }

  {  // Compute P_i
    // We need to solve for P_i in: A_{ii}*P_i + A_{if}*P_f + A_{ie}*P_e = A_{ic}
    // Step one: A_{ic} = 1.0*A_{ic} + (-1.0)*A_{ie}*P_e
    Aic.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, Aie, Pe, 1.0);
    // Step one: A_{ic} = 1.0*A_{ic} + (-1.0)*A_{if}*P_f
    Aic.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, Aif, Pf, 1.0);
    // Step two: P_i = A_{ii}^{-1}*A_{ic}
    Teuchos::SerialDenseSolver<LO, SC> problem;
    problem.setMatrix(Teuchos::rcp(&Aii, false));
    problem.setVectors(Teuchos::rcp(&Pi, false), Teuchos::rcp(&Aic, false));
    problem.factorWithEquilibration(true);
    problem.solve();
    problem.unequilibrateLHS();
  }

}  // ComputeLocalEntries()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CollapseStencil(const int type, const int orientation, const int /* collapseFlags */[3],
                    Array<SC>& stencil) const {
  if (type == 2) {  // Face stencil collapse
    // Let's hope things vectorize well inside this if statement
    if (orientation == 0) {
      for (LO i = 0; i < 9; ++i) {
        stencil[3 * i + 1] = stencil[3 * i] + stencil[3 * i + 1] + stencil[3 * i + 2];
        stencil[3 * i]     = 0;
        stencil[3 * i + 2] = 0;
      }
    } else if (orientation == 1) {
      for (LO i = 0; i < 3; ++i) {
        stencil[9 * i + 3] = stencil[9 * i + 0] + stencil[9 * i + 3] + stencil[9 * i + 6];
        stencil[9 * i + 0] = 0;
        stencil[9 * i + 6] = 0;
        stencil[9 * i + 4] = stencil[9 * i + 1] + stencil[9 * i + 4] + stencil[9 * i + 7];
        stencil[9 * i + 1] = 0;
        stencil[9 * i + 7] = 0;
        stencil[9 * i + 5] = stencil[9 * i + 2] + stencil[9 * i + 5] + stencil[9 * i + 8];
        stencil[9 * i + 2] = 0;
        stencil[9 * i + 8] = 0;
      }
    } else if (orientation == 2) {
      for (LO i = 0; i < 9; ++i) {
        stencil[i + 9]  = stencil[i] + stencil[i + 9] + stencil[i + 18];
        stencil[i]      = 0;
        stencil[i + 18] = 0;
      }
    }
  } else if (type == 1) {  // Edge stencil collapse
    SC tmp1 = 0, tmp2 = 0, tmp3 = 0;
    if (orientation == 0) {
      for (LO i = 0; i < 9; ++i) {
        tmp1 += stencil[0 + 3 * i];
        stencil[0 + 3 * i] = 0;
        tmp2 += stencil[1 + 3 * i];
        stencil[1 + 3 * i] = 0;
        tmp3 += stencil[2 + 3 * i];
        stencil[2 + 3 * i] = 0;
      }
      stencil[12] = tmp1;
      stencil[13] = tmp2;
      stencil[14] = tmp3;
    } else if (orientation == 1) {
      for (LO i = 0; i < 3; ++i) {
        tmp1 += stencil[0 + i] + stencil[9 + i] + stencil[18 + i];
        stencil[0 + i]  = 0;
        stencil[9 + i]  = 0;
        stencil[18 + i] = 0;
        tmp2 += stencil[3 + i] + stencil[12 + i] + stencil[21 + i];
        stencil[3 + i]  = 0;
        stencil[12 + i] = 0;
        stencil[21 + i] = 0;
        tmp3 += stencil[6 + i] + stencil[15 + i] + stencil[24 + i];
        stencil[6 + i]  = 0;
        stencil[15 + i] = 0;
        stencil[24 + i] = 0;
      }
      stencil[10] = tmp1;
      stencil[13] = tmp2;
      stencil[16] = tmp3;
    } else if (orientation == 2) {
      for (LO i = 0; i < 9; ++i) {
        tmp1 += stencil[i];
        stencil[i] = 0;
        tmp2 += stencil[i + 9];
        stencil[i + 9] = 0;
        tmp3 += stencil[i + 18];
        stencil[i + 18] = 0;
      }
      stencil[4]  = tmp1;
      stencil[13] = tmp2;
      stencil[22] = tmp3;
    }
  }
}  // CollapseStencil

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FormatStencil(const LO BlkSize, const Array<bool> /* ghostInterface */, const LO /* ie */, const LO /* je */,
                  const LO /* ke */, const ArrayView<const SC> rowValues, const Array<LO> /* elementNodesPerDir */,
                  const int collapseFlags[3], const std::string stencilType, Array<SC>& stencil)
        const {
  if (stencilType == "reduced") {
    Array<int> fullStencilInds(7);
    fullStencilInds[0] = 4;
    fullStencilInds[1] = 10;
    fullStencilInds[2] = 12;
    fullStencilInds[3] = 13;
    fullStencilInds[4] = 14;
    fullStencilInds[5] = 16;
    fullStencilInds[6] = 22;

    // Create a mask array to automate mesh boundary processing
    Array<int> mask(7);
    int stencilSize = rowValues.size();
    if (collapseFlags[0] + collapseFlags[1] + collapseFlags[2] > 0) {
      if (stencilSize == 1) {
        // The assumption is made that if only one non-zero exists in the row
        // then this represent a Dirichlet BC being enforced.
        // One might want to zero out the corresponding row in the prolongator.
        mask[0] = 1;
        mask[1] = 1;
        mask[2] = 1;
        mask[4] = 1;
        mask[5] = 1;
        mask[6] = 1;
      } else {
        // Here we are looking at Neumann like BC where only a flux is prescribed
        // and less non-zeros will be present in the corresponding row.
        if (collapseFlags[2] == 1 || collapseFlags[2] == 3) {
          mask[0] = 1;
        }
        if (collapseFlags[2] == 2 || collapseFlags[2] == 3) {
          mask[6] = 1;
        }
        if (collapseFlags[1] == 1 || collapseFlags[1] == 3) {
          mask[1] = 1;
        }
        if (collapseFlags[1] == 2 || collapseFlags[1] == 3) {
          mask[5] = 1;
        }
        if (collapseFlags[0] == 1 || collapseFlags[0] == 3) {
          mask[2] = 1;
        }
        if (collapseFlags[0] == 2 || collapseFlags[0] == 3) {
          mask[4] = 1;
        }
      }
    }

    int offset = 0;
    for (int ind = 0; ind < 7; ++ind) {
      if (mask[ind] == 1) {
        for (int dof = 0; dof < BlkSize; ++dof) {
          stencil[BlkSize * fullStencilInds[ind] + dof] = 0.0;
        }
        ++offset;  // The offset is used to stay within rowValues bounds
      } else {
        for (int dof = 0; dof < BlkSize; ++dof) {
          stencil[BlkSize * fullStencilInds[ind] + dof] = rowValues[BlkSize * (ind - offset) + dof];
        }
      }
    }
  } else if (stencilType == "full") {
    // Create a mask array to automate mesh boundary processing
    Array<int> mask(27);
    if (collapseFlags[2] == 1 || collapseFlags[2] == 3) {
      for (int ind = 0; ind < 9; ++ind) {
        mask[ind] = 1;
      }
    }
    if (collapseFlags[2] == 2 || collapseFlags[2] == 3) {
      for (int ind = 0; ind < 9; ++ind) {
        mask[18 + ind] = 1;
      }
    }
    if (collapseFlags[1] == 1 || collapseFlags[1] == 3) {
      for (int ind1 = 0; ind1 < 3; ++ind1) {
        for (int ind2 = 0; ind2 < 3; ++ind2) {
          mask[ind1 * 9 + ind2] = 1;
        }
      }
    }
    if (collapseFlags[1] == 2 || collapseFlags[1] == 3) {
      for (int ind1 = 0; ind1 < 3; ++ind1) {
        for (int ind2 = 0; ind2 < 3; ++ind2) {
          mask[ind1 * 9 + ind2 + 6] = 1;
        }
      }
    }
    if (collapseFlags[0] == 1 || collapseFlags[0] == 3) {
      for (int ind = 0; ind < 9; ++ind) {
        mask[3 * ind] = 1;
      }
    }
    if (collapseFlags[0] == 2 || collapseFlags[0] == 3) {
      for (int ind = 0; ind < 9; ++ind) {
        mask[3 * ind + 2] = 1;
      }
    }

    // Now the stencil is extracted and formated
    int offset = 0;
    for (LO ind = 0; ind < 27; ++ind) {
      if (mask[ind] == 0) {
        for (int dof = 0; dof < BlkSize; ++dof) {
          stencil[BlkSize * ind + dof] = 0.0;
        }
        ++offset;  // The offset is used to stay within rowValues bounds
      } else {
        for (int dof = 0; dof < BlkSize; ++dof) {
          stencil[BlkSize * ind + dof] = rowValues[BlkSize * (ind - offset) + dof];
        }
      }
    }
  }  // stencilTpye
}  // FormatStencil()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetNodeInfo(
    const LO ie, const LO je, const LO ke,
    const Array<LO> elementNodesPerDir,
    int* type, LO& ind, int* orientation) const {
  // Initialize flags with -1 to be able to catch issues easily
  *type = -1, ind = 0, *orientation = -1;
  if ((ke == 0 || ke == elementNodesPerDir[2] - 1) && (je == 0 || je == elementNodesPerDir[1] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1)) {
    // Corner node
    *type = 0;
    ind   = 4 * ke / (elementNodesPerDir[2] - 1) + 2 * je / (elementNodesPerDir[1] - 1) + ie / (elementNodesPerDir[0] - 1);
  } else if (((ke == 0 || ke == elementNodesPerDir[2] - 1) && (je == 0 || je == elementNodesPerDir[1] - 1)) || ((ke == 0 || ke == elementNodesPerDir[2] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1)) || ((je == 0 || je == elementNodesPerDir[1] - 1) && (ie == 0 || ie == elementNodesPerDir[0] - 1))) {
    // Edge node
    *type = 1;
    if (ke > 0) {
      ind += 2 * (elementNodesPerDir[0] - 2 + elementNodesPerDir[1] - 2);
    }
    if (ke == elementNodesPerDir[2] - 1) {
      ind += 4 * (elementNodesPerDir[2] - 2);
    }
    if ((ke == 0) || (ke == elementNodesPerDir[2] - 1)) {  // je or ie edges
      if (je == 0) {                                       // jlo ie edge
        *orientation = 0;
        ind += ie - 1;
      } else if (je == elementNodesPerDir[1] - 1) {  // jhi ie edge
        *orientation = 0;
        ind += 2 * (elementNodesPerDir[1] - 2) + elementNodesPerDir[0] - 2 + ie - 1;
      } else {  // ilo or ihi je edge
        *orientation = 1;
        ind += elementNodesPerDir[0] - 2 + 2 * (je - 1) + ie / (elementNodesPerDir[0] - 1);
      }
    } else {  // ke edges
      *orientation = 2;
      ind += 4 * (ke - 1) + 2 * (je / (elementNodesPerDir[1] - 1)) + ie / (elementNodesPerDir[0] - 1);
    }
  } else if ((ke == 0 || ke == elementNodesPerDir[2] - 1) || (je == 0 || je == elementNodesPerDir[1] - 1) || (ie == 0 || ie == elementNodesPerDir[0] - 1)) {
    // Face node
    *type = 2;
    if (ke == 0) {  // current node is on "bottom" face
      *orientation = 2;
      ind          = (je - 1) * (elementNodesPerDir[0] - 2) + ie - 1;
    } else {
      // add nodes from "bottom" face
      ind += (elementNodesPerDir[1] - 2) * (elementNodesPerDir[0] - 2);
      // add nodes from side faces
      ind += 2 * (ke - 1) * (elementNodesPerDir[1] - 2 + elementNodesPerDir[0] - 2);
      if (ke == elementNodesPerDir[2] - 1) {  // current node is on "top" face
        *orientation = 2;
        ind += (je - 1) * (elementNodesPerDir[0] - 2) + ie - 1;
      } else {  // current node is on a side face
        if (je == 0) {
          *orientation = 1;
          ind += ie - 1;
        } else if (je == elementNodesPerDir[1] - 1) {
          *orientation = 1;
          ind += 2 * (elementNodesPerDir[1] - 2) + elementNodesPerDir[0] - 2 + ie - 1;
        } else {
          *orientation = 0;
          ind += elementNodesPerDir[0] - 2 + 2 * (je - 1) + ie / (elementNodesPerDir[0] - 1);
        }
      }
    }
  } else {
    // Interior node
    *type = 3;
    ind   = (ke - 1) * (elementNodesPerDir[1] - 2) * (elementNodesPerDir[0] - 2) + (je - 1) * (elementNodesPerDir[0] - 2) + ie - 1;
  }
}  // GetNodeInfo()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sh_sort_permute(
    const typename Teuchos::Array<LocalOrdinal>::iterator& first1,
    const typename Teuchos::Array<LocalOrdinal>::iterator& last1,
    const typename Teuchos::Array<LocalOrdinal>::iterator& first2,
    const typename Teuchos::Array<LocalOrdinal>::iterator& /* last2 */) const {
  typedef typename std::iterator_traits<typename Teuchos::Array<LocalOrdinal>::iterator>::difference_type DT;
  DT n = last1 - first1;
  DT m = n / 2;
  DT z = Teuchos::OrdinalTraits<DT>::zero();
  while (m > z) {
    DT max = n - m;
    for (DT j = 0; j < max; j++) {
      for (DT k = j; k >= 0; k -= m) {
        if (first1[first2[k + m]] >= first1[first2[k]])
          break;
        std::swap(first2[k + m], first2[k]);
      }
    }
    m = m / 2;
  }
}

}  // namespace MueLu

#define MUELU_BLACKBOXPFACTORY_SHORT
#endif  // MUELU_BLACKBOXPFACTORY_DEF_HPP
