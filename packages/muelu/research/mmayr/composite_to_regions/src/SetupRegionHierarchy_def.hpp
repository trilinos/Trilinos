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
#ifndef MUELU_SETUPREGIONHIERARCHY_DEF_HPP
#define MUELU_SETUPREGIONHIERARCHY_DEF_HPP

#include <vector>
#include <iostream>

#include <Kokkos_DefaultNode.hpp>

#include <Teuchos_RCP.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"
#include "SetupRegionSmoothers_def.hpp"


#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>
#endif

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ParameterList;

/*! \brief Create coarse level maps with continuous GIDs
 *
 *  The direct solver requires maps with continuous GIDs. Starting from the
 *  coarse level composite maps with discontinuous GIDs, we create a new row map
 *  and a matching column map.
 *
 *  Range and Domain map happen to correspond to the Row map, so we don't have
 *  to deal with them in particular.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void createContinuousCoarseLevelMaps(const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap, ///< row map
                                     const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap, ///< column map
                                     RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& contRowMap, ///< row map with continuous GIDs
                                     RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& contColMap ///< column map with continuous GIDs
                                     )
{
#include "Xpetra_UseShortNamesOrdinal.hpp"
  //   /!\ This function is pure ordinal, no scalar type is passed as input
  //       This means that use only three template paramters and that we
  //       do not use the Scalar dependent short names!

  // Create row map with continuous GIDs
  contRowMap = MapFactory::Build(rowMap->lib(),
                                 rowMap->getGlobalNumElements(),
                                 rowMap->getNodeNumElements(),
                                 rowMap->getIndexBase(),
                                 rowMap->getComm());

  return;
} // createContinuousCoarseLevelMaps



/* Reconstruct coarse-level maps (assuming fully structured grids)
 *
 * We know the regional map on the coarse levels since they are just the
 * row maps of the coarse level operators. Though, we need to re-construct
 * the quasiRegional and composite maps ourselves.
 *
 * We ultimately are only interested in the composite map on the coarsest level.
 * Intermediate levels are dealt with along the way, because we go through all
 * levels recursively.
 *
 * Assumptions:
 * - fully structured grid
 * - only on region per proc and group
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeCoarseLevelMaps(const int maxRegPerGID,
                         Teuchos::ArrayView<LocalOrdinal>  compositeToRegionLIDsFinest,
                         Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regProlong,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regColMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegColMaps,
                         Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                         Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                         Array<Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >& interfaceGIDs,
                         Array<Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >& regionsPerGIDWithGhosts,
                         Array<Teuchos::ArrayRCP<LocalOrdinal> >& regionMatVecLIDs,
                         Array<Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >& regionInterfaceImporter) {

#include "Xpetra_UseShortNames.hpp"

  using MT = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  const GO GO_INV = Teuchos::OrdinalTraits<GO>::invalid();
  const int numLevels = regProlong.size();

  // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // Teuchos::FancyOStream& out = *fancy;
  const int myRank = regProlong[1][0]->getRowMap()->getComm()->getRank();

  Teuchos::Array<LO> coarseCompositeToRegionLIDs;
  Teuchos::ArrayView<LO> compositeToRegionLIDs = compositeToRegionLIDsFinest;
  for(int currentLevel = 1; currentLevel < numLevels; ++currentLevel) {

    // Extracting some basic information about local mesh in composite/region format
    const size_t numFineRegionNodes    = regProlong[currentLevel][0]->getNodeNumRows();
    const size_t numFineCompositeNodes = compositeToRegionLIDs.size();
    const size_t numFineDuplicateNodes = numFineRegionNodes - numFineCompositeNodes;

    const size_t numCoarseRegionNodes  = regProlong[currentLevel][0]->getColMap()->getNodeNumElements();

    // Find the regionLIDs associated with local duplicated nodes
    // This will allow us to later loop only on duplicated nodes
    size_t countComposites = 0, countDuplicates = 0;
    Array<LO> fineDuplicateLIDs(numFineDuplicateNodes);
    for(size_t regionIdx = 0; regionIdx < numFineRegionNodes; ++regionIdx) {
      if(compositeToRegionLIDs[countComposites] == static_cast<LO>(regionIdx)) {
        ++countComposites;
      } else {
        fineDuplicateLIDs[countDuplicates] = regionIdx;
        ++countDuplicates;
      }
    }

    // We gather the coarse GIDs associated with each fine point in the local composite mesh part.
    RCP<Xpetra::Vector<MT,LO,GO,NO> > coarseCompositeGIDs
      = Xpetra::VectorFactory<MT,LO,GO,NO>::Build(regRowImporters[currentLevel - 1][0]->getSourceMap(), false);
    Teuchos::ArrayRCP<MT> coarseCompositeGIDsData = coarseCompositeGIDs->getDataNonConst(0);

    for(size_t compositeNodeIdx = 0; compositeNodeIdx < numFineCompositeNodes; ++compositeNodeIdx) {
      ArrayView<const LO> coarseRegionLID; // Should contain a single value
      ArrayView<const SC> dummyData; // Should contain a single value
      regProlong[currentLevel][0]->getLocalRowView(compositeToRegionLIDs[compositeNodeIdx],
                                                   coarseRegionLID,
                                                   dummyData);
      if(coarseRegionLID.size() == 1) {
        coarseCompositeGIDsData[compositeNodeIdx] = regProlong[currentLevel][0]->getColMap()->getGlobalElement(coarseRegionLID[0]);
      } else {
        coarseCompositeGIDsData[compositeNodeIdx] = -1;
      }
    }

    // We communicate the above GIDs to their duplicate so that we can replace GIDs of the region
    // column map and form the quasiRegion column map.
    Array<RCP<Xpetra::Vector<MT, LO, GO, NO> > > coarseQuasiregionGIDs(1);
    Array<RCP<Xpetra::Vector<MT, LO, GO, NO> > > coarseRegionGIDs(1);
    compositeToRegional(coarseCompositeGIDs,
                        coarseQuasiregionGIDs,
                        coarseRegionGIDs,
                        regRowMaps[currentLevel - 1],
                        regRowImporters[currentLevel - 1]);

    regionsPerGIDWithGhosts[currentLevel] =
      Xpetra::MultiVectorFactory<LO, LO, GO, NO>::Build(regRowMaps[currentLevel][0],
                                                        maxRegPerGID,
                                                        false);
    interfaceGIDs[currentLevel] =
      Xpetra::MultiVectorFactory<GO, LO, GO, NO>::Build(regRowMaps[currentLevel][0],
                                                        maxRegPerGID,
                                                        false);
    Array<ArrayRCP<const LO> > regionsPerGIDWithGhostsFine(maxRegPerGID);
    Array<ArrayRCP<LO> > regionsPerGIDWithGhostsCoarse(maxRegPerGID);
    Array<ArrayRCP<GO> > interfaceGIDsCoarse(maxRegPerGID);
    for(size_t idx = 0; idx < static_cast<size_t>(maxRegPerGID); ++idx) {
      regionsPerGIDWithGhostsFine[idx]   = regionsPerGIDWithGhosts[currentLevel - 1]->getData(idx);
      regionsPerGIDWithGhostsCoarse[idx] = regionsPerGIDWithGhosts[currentLevel]->getDataNonConst(idx);
      interfaceGIDsCoarse[idx]          = interfaceGIDs[currentLevel]->getDataNonConst(idx);
      for(size_t coarseIdx = 0;
          coarseIdx < regionsPerGIDWithGhosts[currentLevel]->getLocalLength(); ++coarseIdx) {
        regionsPerGIDWithGhostsCoarse[idx][coarseIdx] = -1;
        interfaceGIDsCoarse[idx][coarseIdx] = 0;
      }
    }

    for(size_t fineIdx = 0; fineIdx < numFineRegionNodes; ++fineIdx) {
      ArrayView<const LO> coarseRegionLID; // Should contain a single value
      ArrayView<const SC> dummyData;       // Should contain a single value
      regProlong[currentLevel][0]->getLocalRowView(fineIdx,
                                                   coarseRegionLID,
                                                   dummyData);
      const LO coarseIdx = coarseRegionLID[0];

      // Now fill regionPerGIDWithGhostsCoarse[:][coarseRegionLID]
      // with data from regionPerGIDWithGhostsFine[:][fineIdx].
      // The problem is we might have more then maxRegPerGID on currentLevel
      // then on currentLevel-1... we might need to do a union or something?
      // I guess technically using the restriction operator here would be more
      // helpful than the prolongator, this way we could find all the fine interface
      // points easily and compute the union of the PIDs they belong too.
      // This can actually be done extracting the local matrix and compute the transpose.
      // For now let us assume that maxRegPerGID is constant and hope for the best.
      LO countFinePIDs   = 0;
      LO countCoarsePIDs = 0;
      for(LO idx = 0; idx < maxRegPerGID; ++idx) {
        if(-1 < regionsPerGIDWithGhostsFine[idx][fineIdx]) {++countFinePIDs;}
        if(-1 < regionsPerGIDWithGhostsCoarse[idx][coarseIdx]) {++countCoarsePIDs;}
      }

      if(countCoarsePIDs < countFinePIDs) {
        for(LO idx = 0; idx < countFinePIDs; ++idx) {
          regionsPerGIDWithGhostsCoarse[idx][coarseIdx] = regionsPerGIDWithGhostsFine[idx][fineIdx];
          if(regionsPerGIDWithGhostsCoarse[idx][coarseIdx] == myRank) {
            interfaceGIDsCoarse[idx][coarseIdx] = regRowMaps[currentLevel][0]->getGlobalElement(coarseIdx);
          }
        }
      }
    }

    Array<GO> fineRegionDuplicateCoarseLIDs(numFineDuplicateNodes);
    Array<GO> fineRegionDuplicateCoarseGIDs(numFineDuplicateNodes);
    for(size_t duplicateIdx = 0; duplicateIdx < numFineDuplicateNodes; ++duplicateIdx) {
      ArrayView<const LO> coarseRegionLID; // Should contain a single value
      ArrayView<const SC> dummyData; // Should contain a single value
      regProlong[currentLevel][0]->getLocalRowView(fineDuplicateLIDs[duplicateIdx],
                                                   coarseRegionLID,
                                                   dummyData);
      fineRegionDuplicateCoarseLIDs[duplicateIdx] = regProlong[currentLevel][0]->getColMap()->getGlobalElement(coarseRegionLID[0]);
      fineRegionDuplicateCoarseGIDs[duplicateIdx] = (coarseQuasiregionGIDs[0]->getDataNonConst(0))[fineDuplicateLIDs[duplicateIdx]];
    }

    // Create the coarseQuasiregRowMap, it will be based on the coarseRegRowMap
    LO countCoarseComposites = 0;
    coarseCompositeToRegionLIDs.resize(numCoarseRegionNodes);
    Array<GO> coarseQuasiregRowMapData = regProlong[currentLevel][0]->getColMap()->getNodeElementList();
    Array<GO> coarseCompRowMapData(numCoarseRegionNodes, -1);
    for(size_t regionIdx = 0; regionIdx < numCoarseRegionNodes; ++regionIdx) {
      const GO initialValue = coarseQuasiregRowMapData[regionIdx];
      for(size_t duplicateIdx = 0; duplicateIdx < numFineDuplicateNodes; ++duplicateIdx) {
        if((initialValue == fineRegionDuplicateCoarseLIDs[duplicateIdx]) &&
           (fineRegionDuplicateCoarseGIDs[duplicateIdx] < coarseQuasiregRowMapData[regionIdx]) &&
           (-1 < fineRegionDuplicateCoarseGIDs[duplicateIdx])){
          coarseQuasiregRowMapData[regionIdx] = fineRegionDuplicateCoarseGIDs[duplicateIdx];
        }
      }
      if(initialValue == coarseQuasiregRowMapData[regionIdx]) {
        coarseCompRowMapData[countCoarseComposites] = coarseQuasiregRowMapData[regionIdx];
        coarseCompositeToRegionLIDs[countCoarseComposites] = regionIdx;
        ++countCoarseComposites;
      }
    }
    coarseCompRowMapData.resize(countCoarseComposites);
    coarseCompositeToRegionLIDs.resize(countCoarseComposites);

    // We are now ready to fill up the outputs
    regRowMaps[currentLevel][0] = Teuchos::rcp_const_cast<Map>(regProlong[currentLevel][0]->getColMap());
    regColMaps[currentLevel][0] = Teuchos::rcp_const_cast<Map>(regProlong[currentLevel][0]->getColMap());
    quasiRegRowMaps[currentLevel][0] = MapFactory::Build(regProlong[currentLevel][0]->getColMap()->lib(),
                                                         GO_INV,
                                                         coarseQuasiregRowMapData(),
                                                         regProlong[currentLevel][0]->getColMap()->getIndexBase(),
                                                         regProlong[currentLevel][0]->getColMap()->getComm());
    quasiRegColMaps[currentLevel][0] = quasiRegRowMaps[currentLevel][0];
    compRowMaps[currentLevel] = MapFactory::Build(regProlong[currentLevel][0]->getColMap()->lib(),
                                                  GO_INV,
                                                  coarseCompRowMapData(),
                                                  regProlong[currentLevel][0]->getColMap()->getIndexBase(),
                                                  regProlong[currentLevel][0]->getColMap()->getComm());
    regRowImporters[currentLevel][0] = ImportFactory::Build(compRowMaps[currentLevel], quasiRegRowMaps[currentLevel][0]);

    // Now generate matvec data
    SetupMatVec(interfaceGIDs[currentLevel], regionsPerGIDWithGhosts[currentLevel],
                regRowMaps[currentLevel], regRowImporters[currentLevel],
                regionMatVecLIDs[currentLevel], regionInterfaceImporter[currentLevel]);

    // Finally reset compositeToRegionLIDs
    compositeToRegionLIDs = coarseCompositeToRegionLIDs();
  } // Loop over numLevels
} // MakeCoarseLevelMaps


// Form the composite coarse level operator
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeCoarseCompositeOperator(const int maxRegPerProc,
                                 RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& compRowMap,
                                 std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& quasiRegRowMap,
                                 std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& quasiRegColMap,
                                 std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >& regRowImporter,
                                 std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regMatrix,
                                 Array<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > >& regCoarseCoordinates,
                                 RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& coarseCompOp,
                                 RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& compCoarseCoordinates,
                                 const bool makeCompCoords)
{
#include "Xpetra_UseShortNames.hpp"
  using CoordType = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  coarseCompOp = MatrixFactory::Build(compRowMap,
                                       // This estimate is very conservative and probably costs us lots of memory...
                                      8*regMatrix[0]->getCrsGraph()->getNodeMaxNumRowEntries());
  regionalToComposite(regMatrix,
                      quasiRegRowMap, quasiRegColMap,
                      regRowImporter, Xpetra::ADD,
                      coarseCompOp);

  const int dofsPerNode = regMatrix[0]->GetFixedBlockSize();
  coarseCompOp->SetFixedBlockSize(dofsPerNode);
  coarseCompOp->setObjectLabel("coarse composite operator");

  // Create coarse composite coordinates for repartitioning
  if(makeCompCoords) {
    for(int regIdx = 0; regIdx < maxRegPerProc; ++regIdx) {
      const int check       = regMatrix[regIdx]->getRowMap()->getNodeNumElements()
        % regCoarseCoordinates[regIdx]->getMap()->getNodeNumElements();
      TEUCHOS_ASSERT(check == 0);
    }

    RCP<Map> compCoordMap;
    std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > regCoordImporter(maxRegPerProc);
    if(dofsPerNode == 1) {
      compCoordMap = compRowMap;
      for(int regionIdx = 0; regionIdx < maxRegPerProc; ++regionIdx) {
        regCoordImporter[regionIdx] = regRowImporter[regionIdx];
      }
    } else {
      using size_type = typename Teuchos::Array<GlobalOrdinal>::size_type;
      Array<GlobalOrdinal> compCoordMapData(compRowMap->getNodeNumElements() / dofsPerNode);
      ArrayView<const GlobalOrdinal> compRowMapData = compRowMap->getNodeElementList();
      for(size_type nodeIdx = 0; nodeIdx < compCoordMapData.size(); ++nodeIdx) {
        compCoordMapData[nodeIdx] = compRowMapData[nodeIdx*dofsPerNode] / dofsPerNode;
      }
      compCoordMap = MapFactory::Build(compRowMap->lib(),
                                       compRowMap->getGlobalNumElements() / dofsPerNode,
                                       compCoordMapData(),
                                       compRowMap->getIndexBase(),
                                       compRowMap->getComm());

      std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > quasiRegCoordMap(maxRegPerProc);
      for(int regionIdx = 0; regionIdx < maxRegPerProc; ++regionIdx) {
        Array<GlobalOrdinal> quasiRegCoordMapData(quasiRegRowMap[regionIdx]->getNodeNumElements() / dofsPerNode);
        ArrayView<const GlobalOrdinal> quasiRegRowMapData = quasiRegRowMap[regionIdx]->getNodeElementList();
        for(size_type nodeIdx = 0; nodeIdx < quasiRegCoordMapData.size(); ++nodeIdx) {
          quasiRegCoordMapData[nodeIdx] = quasiRegRowMapData[nodeIdx*dofsPerNode] / dofsPerNode;
        }
        quasiRegCoordMap[regionIdx] = MapFactory::Build(quasiRegRowMap[regionIdx]->lib(),
                                                        quasiRegRowMap[regionIdx]->getGlobalNumElements() / dofsPerNode,
                                                        quasiRegCoordMapData(),
                                                        quasiRegRowMap[regionIdx]->getIndexBase(),
                                                        quasiRegRowMap[regionIdx]->getComm());
        regCoordImporter[regionIdx] = ImportFactory::Build(compCoordMap, quasiRegCoordMap[regionIdx]);
      }
    }
    compCoarseCoordinates = Xpetra::MultiVectorFactory<CoordType, LocalOrdinal, GlobalOrdinal, Node>
      ::Build(compCoordMap, regCoarseCoordinates[0]->getNumVectors());
    TEUCHOS_ASSERT(Teuchos::nonnull(compCoarseCoordinates));

    // The following looks like regionalToComposite for Vector
    // but it is a bit different since we do not want to add
    // entries in the coordinate MultiVector as we would for
    // a solution or residual vector.
    RCP<Xpetra::MultiVector<CoordType, LocalOrdinal, GlobalOrdinal, Node>>
      quasiRegCoarseCoordinates;
    for(int grpIdx = 0; grpIdx < maxRegPerProc; ++grpIdx) {
      quasiRegCoarseCoordinates = regCoarseCoordinates[grpIdx];
      TEUCHOS_ASSERT(Teuchos::nonnull(quasiRegCoarseCoordinates));
      quasiRegCoarseCoordinates->replaceMap(regCoordImporter[grpIdx]->getTargetMap());
      compCoarseCoordinates->doExport(*quasiRegCoarseCoordinates, *(regCoordImporter[grpIdx]), Xpetra::INSERT);
    }
  }
} // MakeCoarseCompositeOperator


/* Create a direct solver for a composite operator
 *
 * Create the solver object and compute symbolic and numeric factorization.
 * Finally, the solver object will be ready to be applied during the V-cycle call.
 *
 * \note For now, we're limited to Tpetra/Amesos2. From Amesos2, we use KLU as direct solver.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Amesos2::Solver<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >
MakeCompositeDirectSolver(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& compOp)
{
  using Tpetra_CrsMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Utilities = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Teuchos::TimeMonitor;

  RCP<Amesos2::Solver<Tpetra_CrsMatrix, Tpetra_MultiVector> > coarseSolver;
  {
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeCompositeDirectSolver: 1 - Setup")));

    // convert matrix to Tpetra
    RCP<Tpetra_CrsMatrix> tMat = Utilities::Op2NonConstTpetraCrs(compOp);

    // Amesos2-specific key phrase that denote smoother type
    std::string amesos2SolverName = "KLU2";
    TEUCHOS_ASSERT(Amesos2::query(amesos2SolverName));
    coarseSolver = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(amesos2SolverName, tMat);

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist(amesos2SolverName).set("IsContiguous", false, "Are GIDs Contiguous");
    coarseSolver->setParameters(Teuchos::rcpFromRef(amesos2_params));
  }

  {
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeCompositeDirectSolver: 2 - Factorization")));

    coarseSolver->symbolicFactorization();
    coarseSolver->numericFactorization();
  }

  return coarseSolver;
} // MakeCorseCompositeDirectSolver

/* Create an AMG hierarchy for a composite operator
 *
 * Create the hierarchy object and perform the multigrid setup.
 * Finally, the hierarhcy object will be ready to be applied during the region MG V-cycle call.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
MakeCompositeAMGHierarchy(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& compOp,
                          const std::string& xmlFileName,
                          RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > coordinates)
{
#include "MueLu_UseShortNames.hpp"
  using coordinates_type = typename Teuchos::ScalarTraits<Scalar>::coordinateType;

  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);
  *fos << "Attempting to setup AMG hierarchy for the composite coarse grid problem" << std::endl;

  // Get parameter list for AMG hierarchy
  RCP<ParameterList> mueluParams = Teuchos::rcp(new ParameterList());
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, mueluParams.ptr(),
      *compOp->getRowMap()->getComm());


  // Get the user data sublist
  const std::string userName = "user data";
  Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);

  // Add nullspace information
  {
    // Compute nullspace
    RCP<MultiVector> nullspace;
    if((compOp->GetFixedBlockSize() == 1) || Teuchos::is_null(coordinates)) { // Scalar problem, constant nullspace
      nullspace = MultiVectorFactory::Build(compOp->getRowMap(), 1);
      nullspace->putScalar(one);
    } else if(compOp->GetFixedBlockSize() == 2) { // 2D Elasticity
      nullspace = MultiVectorFactory::Build(compOp->getRowMap(), 3);
      Array<ArrayRCP<SC> > nullspaceData(3);
      Array<ArrayRCP<const coordinates_type> > coordinateData(2);

      // Calculate center
      const coordinates_type cx = coordinates->getVector(0)->meanValue();
      const coordinates_type cy = coordinates->getVector(1)->meanValue();

      coordinateData[0] = coordinates->getData(0);
      coordinateData[1] = coordinates->getData(1);

      for(int vecIdx = 0; vecIdx < 3; ++vecIdx) {
        nullspaceData[vecIdx] = nullspace->getDataNonConst(vecIdx);
      }

      for(size_t nodeIdx = 0; nodeIdx < coordinates->getLocalLength(); ++nodeIdx) {
        // translations
        nullspaceData[0][2*nodeIdx + 0] = one;
        nullspaceData[1][2*nodeIdx + 1] = one;

        // rotation about z axis
        nullspaceData[2][2*nodeIdx + 0] = -(coordinateData[1][nodeIdx] - cy);
        nullspaceData[2][2*nodeIdx + 1] =  (coordinateData[0][nodeIdx] - cx);
      }

    } else if(compOp->GetFixedBlockSize() == 3) { // 3D Elasticity
      nullspace = MultiVectorFactory::Build(compOp->getRowMap(), 6);
      Array<ArrayRCP<SC> > nullspaceData(6);
      Array<ArrayRCP<const coordinates_type> > coordinateData(3);

      // Calculate center
      const coordinates_type cx = coordinates->getVector(0)->meanValue();
      const coordinates_type cy = coordinates->getVector(1)->meanValue();
      const coordinates_type cz = coordinates->getVector(2)->meanValue();

      coordinateData[0] = coordinates->getData(0);
      coordinateData[1] = coordinates->getData(1);
      coordinateData[2] = coordinates->getData(2);

      for(int vecIdx = 0; vecIdx < 6; ++vecIdx) {
        nullspaceData[vecIdx] = nullspace->getDataNonConst(vecIdx);
      }

      for(size_t nodeIdx = 0; nodeIdx < coordinates->getLocalLength(); ++nodeIdx) {
        // translations
        nullspaceData[0][3*nodeIdx + 0] = one;
        nullspaceData[1][3*nodeIdx + 1] = one;
        nullspaceData[2][3*nodeIdx + 2] = one;

        // rotation about z axis
        nullspaceData[3][3*nodeIdx + 0] = -(coordinateData[1][nodeIdx] - cy);
        nullspaceData[3][3*nodeIdx + 1] =  (coordinateData[0][nodeIdx] - cx);

        // rotation about x axis
        nullspaceData[4][3*nodeIdx + 1] = -(coordinateData[2][nodeIdx] - cz);
        nullspaceData[4][3*nodeIdx + 2] =  (coordinateData[1][nodeIdx] - cy);

        // rotation about y axis
        nullspaceData[5][3*nodeIdx + 0] =  (coordinateData[2][nodeIdx] - cz);
        nullspaceData[5][3*nodeIdx + 2] = -(coordinateData[0][nodeIdx] - cx);
      }
    }

    // Equalize norms of all vectors to that of the first one
    // We do not normalize them as a vector of ones seems nice
    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms2(nullspace->getNumVectors());
    nullspace->norm2(norms2);
    Teuchos::Array<Scalar> norms2scalar(nullspace->getNumVectors());
    for (size_t vectorIdx = 0; vectorIdx < nullspace->getNumVectors(); ++vectorIdx) {
      norms2scalar[vectorIdx] = norms2[0] / norms2[vectorIdx];
    }
    nullspace->scale(norms2scalar);

    // Insert into parameter list
    userParamList.set("Nullspace", nullspace);
  }

  // Add coordinate information for rebalancing
  {
    if(Teuchos::nonnull(coordinates)) {
      userParamList.set("Coordinates", coordinates);
    }
  }

  // Create an AMG hierarchy based on the composite coarse level operator from the region MG scheme
  RCP<Hierarchy> compOpHiearchy = MueLu::CreateXpetraPreconditioner(compOp, *mueluParams);

  // We will use the hiearchy as a solver
  compOpHiearchy->IsPreconditioner(false);
  compOpHiearchy->SetVerbLevel(MueLu::VERB_NONE);

  return compOpHiearchy;
} // MakeCompositeAMGHierarchy

  // Make interface scaling factors recursively
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeInterfaceScalingFactors(const int maxRegPerProc,
                                 const int numLevels,
                                 Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                                 Array<Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regInterfaceScalings,
                                 Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                                 Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                                 Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps)
{
#include "Xpetra_UseShortNames.hpp"
  // std::cout << compRowMaps[0]->getComm()->getRank() << " | Computing interface scaling factors ..." << std::endl;

  const SC SC_ONE = Teuchos::ScalarTraits<SC>::one();

  TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

  for (int l = 0; l < numLevels; l++) {
    // initialize region vector with all ones.
    for (int j = 0; j < maxRegPerProc; j++) {
      regInterfaceScalings[l][j] = VectorFactory::Build(regRowMaps[l][j]);
      regInterfaceScalings[l][j]->putScalar(SC_ONE);
    }

    // transform to composite layout while adding interface values via the Export() combine mode
    RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(compRowMaps[l], true);
    regionalToComposite(regInterfaceScalings[l], compInterfaceScalingSum, regRowImporters[l]);

    /* transform composite layout back to regional layout. Now, GIDs associated
     * with region interface should carry a scaling factor (!= 1).
     */
    Array<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc); // Is that vector really needed?
    compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                        regInterfaceScalings[l],
                        regRowMaps[l], regRowImporters[l]);
  }
} // MakeInterfaceScalingFactors


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionHierarchy(const int maxRegPerProc,
                           const int numDimensions,
                           const Array<Array<int> > lNodesPerDim,
                           const Array<std::string> aggregationRegionType,
                           RCP<Teuchos::ParameterList>& interfaceParams,
                           const std::string xmlFileName,
                           Array<RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& nullspace,
                           Array<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > >& coordinates,
                           std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regionGrpMats,
                           const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > colMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedColMapPerGrp,
                           const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp,
                           Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                           Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compColMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regColMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegColMaps,
                           Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regMatrices,
                           Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regProlong,
                           Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                           Array<Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regInterfaceScalings,
                           Array<Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >& interfaceGIDs,
                           Array<Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >& regionsPerGIDWithGhosts,
                           Array<Teuchos::ArrayRCP<LocalOrdinal> >& regionMatVecLIDs,
                           Array<Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >& regionInterfaceImporter,
                           const int maxRegPerGID,
                           ArrayView<LocalOrdinal> compositeToRegionLIDs,
                           RCP<Teuchos::ParameterList>& coarseSolverData,
                           Array<RCP<Teuchos::ParameterList> >& smootherParams,
                           RCP<Teuchos::ParameterList> hierarchyData,
                           const bool keepCoarseCoords)
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  // This monitor times everything and gets the overall setting cost
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy")));

  using Hierarchy = MueLu::Hierarchy<SC, LO, GO, NO>;
  using Utilities = MueLu::Utilities<SC, LO, GO, NO>;
  using DirectCoarseSolver = Amesos2::Solver<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >;

  // std::cout << mapComp->getComm()->getRank() << " | Setting up MueLu hierarchies ..." << std::endl;
  int numLevels = 0;

  // A hierarchy for each group
  std::vector<RCP<Hierarchy> > regGrpHierarchy(maxRegPerProc);

  RCP<TimeMonitor> tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: GrpHierarchy")));
  for (int j = 0; j < maxRegPerProc; j++) {

    /* Set number of nodes per processor per dimension
     *
     * We don't use the number of owned nodes provided on input.
     * Use the region dimensions instead. This is the right thing to do
     * since duplication of interface nodes has added duplicated nodes to those regions
     * where inpData_ownedX/inpData_ownedY and inpData_regionX/inpData_regionY have been different on input.
     */

    // Read MueLu parameter list form xml file
    RCP<ParameterList> mueluParams = Teuchos::rcp(new ParameterList());
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, mueluParams.ptr(), *mapComp->getComm());

    // Insert region-specific data into parameter list
    const std::string userName = "user data";
    Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);
    userParamList.set<int>        ("int numDimensions", numDimensions);
    userParamList.set<Array<LO> > ("Array<LO> lNodesPerDim", lNodesPerDim[j]);
    userParamList.set<std::string>("string aggregationRegionType", aggregationRegionType[j]);
    userParamList.set<Array<LO> > ("Array<LO> nodeOnInterface", interfaceParams->get<Array<LO> >("interfaces: interface nodes"));
    userParamList.set<Array<LO> > ("Array<LO> interfacesDimensions", interfaceParams->get<Array<LO> >("interfaces: nodes per dimensions"));
    if(Teuchos::nonnull(coordinates[j])) {
      userParamList.set("Coordinates", coordinates[j]);
    }
    if(Teuchos::nonnull(nullspace[j])) {
      userParamList.set("Nullspace", nullspace[j]);
    }

    // Setup hierarchy
    regGrpHierarchy[j] = MueLu::CreateXpetraPreconditioner(regionGrpMats[j], *mueluParams);
  }

  // std::cout << mapComp->getComm()->getRank() << " | Resize containers..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: ExtractGrpData")));
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
    interfaceGIDs.resize(numLevels);
    regionsPerGIDWithGhosts.resize(numLevels);
    regionMatVecLIDs.resize(numLevels);
    regionInterfaceImporter.resize(numLevels);
    smootherParams.resize(numLevels);

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

      // Also doing some initialization in the smootherParams
      if(l > 0) {smootherParams[l] = rcp(new Teuchos::ParameterList(*smootherParams[0]));}
    }
  }

  // std::cout << mapComp->getComm()->getRank() << " | Fill fine level containers..." << std::endl;

  // Fill fine level with our data
  {
    compRowMaps[0]     = mapComp;
    quasiRegRowMaps[0] = rowMapPerGrp;
    quasiRegColMaps[0] = colMapPerGrp;
    regRowMaps[0][0]   = revisedRowMapPerGrp[0];
    regColMaps[0][0]   = revisedColMapPerGrp[0];
    regRowImporters[0] = rowImportPerGrp;
    regMatrices[0]     = regionGrpMats;

    /* MueLu stores prolongator on coarse level, so there is no prolongator
     * on the fine level. To have level containers of the same size, let's
     * just put in dummy data
     */
    std::vector<RCP<Matrix> > fineLevelProlong(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; ++j)
      fineLevelProlong[j] = Teuchos::null;
    regProlong[0] = fineLevelProlong;
  }

  // std::cout << mapComp->getComm()->getRank() << " | Fill coarser level containers..." << std::endl;

  /* Get coarse level matrices and prolongators from MueLu hierarchy
   * Note: fine level has been dealt with previously, so we start at level 1 here.
   */
  using real_type = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  using realvaluedmultivector_type = Xpetra::MultiVector<real_type, LocalOrdinal, GlobalOrdinal, Node>;
  Array<RCP<realvaluedmultivector_type> > regCoarseCoordinates(maxRegPerProc);
  for (int l = 1; l < numLevels; ++l) { // Note: we start at level 1 (which is the first coarse level)
    for (int j = 0; j < maxRegPerProc; ++j) {
      RCP<MueLu::Level> level = regGrpHierarchy[j]->GetLevel(l);

      if(keepCoarseCoords && (l == numLevels - 1)) {
        regCoarseCoordinates[j] = level->Get<RCP<realvaluedmultivector_type> >("Coordinates2", MueLu::NoFactory::get());
      }

      regProlong[l][j]  = level->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
      regMatrices[l][j] = level->Get<RCP<Matrix> >("A", MueLu::NoFactory::get());

      regRowMaps[l][j]  = Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> >(regMatrices[l][j]->getRowMap());
      regColMaps[l][j]  = Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> >(regMatrices[l][j]->getColMap());
    }

    // Create residual and solution vectors and cache them for vCycle apply
    std::string levelName("level");
    levelName += std::to_string(l);
    ParameterList& levelList = hierarchyData->sublist(levelName, false, "list of data on current level");
    Teuchos::Array<RCP<Vector> > regRes(maxRegPerProc), regSol(maxRegPerProc);
    createRegionalVector(regRes, revisedRowMapPerGrp);
    createRegionalVector(regSol, revisedRowMapPerGrp);
    levelList.set<Teuchos::Array<RCP<Vector> > >("residual", regRes, "Cached residual vector");
    levelList.set<Teuchos::Array<RCP<Vector> > >("solution", regSol, "Cached solution vector");
  }

  // std::cout << mapComp->getComm()->getRank() << " | MakeCoarseLevelMaps ..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: MakeCoarseLevel")));
  MakeCoarseLevelMaps(maxRegPerGID,
                      compositeToRegionLIDs,
                      regProlong,
                      regRowMaps,
                      regColMaps,
                      quasiRegRowMaps,
                      quasiRegColMaps,
                      compRowMaps,
                      regRowImporters,
                      interfaceGIDs,
                      regionsPerGIDWithGhosts,
                      regionMatVecLIDs,
                      regionInterfaceImporter);

  // std::cout << mapComp->getComm()->getRank() << " | MakeInterfaceScalingFactors ..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: MakeInterfaceScaling")));
  MakeInterfaceScalingFactors(maxRegPerProc,
                              numLevels,
                              compRowMaps,
                              regInterfaceScalings,
                              regRowMaps,
                              regRowImporters,
                              quasiRegRowMaps);

  // std::cout << mapComp->getComm()->getRank() << " | Setup smoothers ..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: SmootherSetup")));
  // Only set the smoother up to the last but one level
  // if we want to use a smoother on the coarse level
  // we will handle that separately with "coarse solver type"
  for(int levelIdx = 0; levelIdx < numLevels - 1; ++levelIdx) {
    smootherParams[levelIdx]->set("smoother: level", levelIdx);
    smootherSetup(smootherParams[levelIdx], regRowMaps[levelIdx],
                  regMatrices[levelIdx], regInterfaceScalings[levelIdx],
                  regRowImporters[levelIdx]);
  }

  // std::cout << mapComp->getComm()->getRank() << " | CreateCoarseSolver ..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: CreateCoarseSolver")));
  const std::string coarseSolverType = coarseSolverData->get<std::string>("coarse solver type");
  if (coarseSolverType == "smoother") {
    // Set the smoother on the coarsest level.
    const std::string smootherXMLFileName = coarseSolverData->get<std::string>("smoother xml file");
    RCP<ParameterList> coarseSmootherParams = smootherParams[numLevels - 1];
    Teuchos::updateParametersFromXmlFileAndBroadcast(smootherXMLFileName, coarseSmootherParams.ptr(), *mapComp->getComm());
    coarseSmootherParams->set("smoother: level", numLevels - 1);
    coarseSmootherParams->print();
    smootherSetup(smootherParams[numLevels - 1], regRowMaps[numLevels - 1],
                  regMatrices[numLevels - 1], regInterfaceScalings[numLevels - 1],
                  regRowImporters[numLevels - 1]);
  } else if( (coarseSolverType == "direct") || (coarseSolverType == "amg") ) {
    // A composite coarse matrix is needed

    // std::cout << mapComp->getComm()->getRank() << " | MakeCoarseCompositeOperator ..." << std::endl;

    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > coarseCompOp;
    RCP<realvaluedmultivector_type> compCoarseCoordinates;
    MakeCoarseCompositeOperator(maxRegPerProc,
                                compRowMaps[numLevels - 1],
                                quasiRegRowMaps[numLevels - 1],
                                quasiRegColMaps[numLevels - 1],
                                regRowImporters[numLevels - 1],
                                regMatrices[numLevels - 1],
                                regCoarseCoordinates,
                                coarseCompOp,
                                compCoarseCoordinates,
                                keepCoarseCoords);

    coarseSolverData->set<RCP<const Map> >("compCoarseRowMap", coarseCompOp->getRowMap());

    // std::cout << mapComp->getComm()->getRank() << " | MakeCoarseCompositeSolver ..." << std::endl;
    if (coarseSolverType == "direct") {
      RCP<DirectCoarseSolver> coarseDirectSolver = MakeCompositeDirectSolver(coarseCompOp);
      coarseSolverData->set<RCP<DirectCoarseSolver> >("direct solver object", coarseDirectSolver);
    } else if (coarseSolverType == "amg") {
      if(keepCoarseCoords == false) {
        std::cout << "WARNING: you requested a coarse AMG solver but you did not request coarse coordinates to be kept, repartitioning is not possible!" << std::endl;
      }
      std::string amgXmlFileName = coarseSolverData->get<std::string>("amg xml file");
      RCP<Hierarchy> coarseAMGHierarchy = MakeCompositeAMGHierarchy(coarseCompOp, amgXmlFileName, compCoarseCoordinates);
      coarseSolverData->set<RCP<Hierarchy> >("amg hierarchy object", coarseAMGHierarchy);
    }
  } else  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(false, "Unknown coarse solver type.");
  }

} // createRegionHierarchy


// Wrapper to be used from the Matlab-based driver
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionHierarchy(const int maxRegPerProc,
                           const int numDimensions,
                           const Array<Array<int> > lNodesPerDim,
                           const Array<std::string> aggregationRegionType,
                           const std::string xmlFileName,
                           Array<RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& nullspace,
                           Array<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > >& coordinates,
                           std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regionGrpMats,
                           const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > colMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedColMapPerGrp,
                           const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp,
                           Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                           Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compColMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regColMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegColMaps,
                           Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regMatrices,
                           Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regProlong,
                           Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                           Array<Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regInterfaceScalings,
                           Array<RCP<Teuchos::ParameterList> >& smootherParams
                           )
{
  // Define dummy values
  const int maxRegPerGID = 0;
  ArrayView<LocalOrdinal> compositeToRegionLIDs = {};
  RCP<ParameterList> coarseSolverParams = rcp(new ParameterList("Coarse solver parameters"));
  coarseSolverParams->set<bool>("use coarse solver", true);
  RCP<ParameterList> dummy = Teuchos::parameterList();
  Array<RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > > interfaceGIDsPerLevel;
  Array<RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > > regionsPerGIDWithGhostsPerLevel;
  Array<ArrayRCP<LocalOrdinal> > regionMatVecLIDsPerLevel;
  Array<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > regionInterfaceImporterPerLevel;

  // Call the actual routine
  createRegionHierarchy(maxRegPerProc, numDimensions, lNodesPerDim,
                        aggregationRegionType, dummy, xmlFileName, nullspace, coordinates,
                        regionGrpMats, mapComp, rowMapPerGrp, colMapPerGrp, revisedRowMapPerGrp,
                        revisedColMapPerGrp, rowImportPerGrp, compRowMaps, compColMaps, regRowMaps,
                        regColMaps, quasiRegRowMaps, quasiRegColMaps, regMatrices, regProlong,
                        regRowImporters, regInterfaceScalings, interfaceGIDsPerLevel,
                        regionsPerGIDWithGhostsPerLevel, regionMatVecLIDsPerLevel,
                        regionInterfaceImporterPerLevel, maxRegPerGID,
                        compositeToRegionLIDs, coarseSolverParams, smootherParams, dummy, false);
}


//! Recursive V-cycle in region fashion
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void vCycle(const int l, ///< ID of current level
            const int numLevels, ///< Total number of levels
            Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& fineRegX, ///< solution
            Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > fineRegB, ///< right hand side
            Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > > regMatrices, ///< Matrices in region layout
            Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > > regProlong, ///< Prolongators in region layout
            Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > compRowMaps, ///< composite maps
            Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > > quasiRegRowMaps, ///< quasiRegional row maps
            Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > > regRowMaps, ///< regional row maps
            Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > > regRowImporters, ///< regional row importers
            Array<Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > > regInterfaceScalings, ///< regional interface scaling factors
            Array<RCP<Teuchos::ParameterList> > smootherParams, ///< region smoother parameter list
	    bool& zeroInitGuess,
            RCP<ParameterList> coarseSolverData = Teuchos::null,
            RCP<ParameterList> hierarchyData = Teuchos::null)
{
#include "MueLu_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  // Get max number of regions per process
  const int maxRegPerProc = fineRegX.size();

  if (l < numLevels - 1) { // fine or intermediate levels

//    std::cout << "level: " << l << std::endl;

    // extract data from hierarchy parameterlist
    std::string levelName("level" + std::to_string(l));
    ParameterList levelList;
    bool useCachedVectors = false;
    // if(Teuchos::nonnull(hierarchyData) &&  hierarchyData->isSublist(levelName)) {
    //   levelList = hierarchyData->sublist(levelName);
    //   if(levelList.isParameter("residual") && levelList.isParameter("solution")) {
    //     useCachedVectors = true;
    //   }
    // }

    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 1 - pre-smoother")));

    // pre-smoothing
    smootherApply(smootherParams[l], fineRegX, fineRegB, regMatrices[l],
                  regRowMaps[l], regRowImporters[l], zeroInitGuess);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 2 - compute residual")));

    const bool useFastMatVec = smootherParams[l]->get<bool>("Use fast MatVec");

    Array<RCP<Vector> > regRes(maxRegPerProc);
    if(useCachedVectors) {
      regRes = levelList.get<Teuchos::Array<RCP<Vector> > >("residual");
    } else {
      createRegionalVector(regRes, regRowMaps[l]);
    }
    if (useFastMatVec)
      computeResidual(regRes, fineRegX, fineRegB, regMatrices[l], *smootherParams[l]);
    else
      computeResidual(regRes, fineRegX, fineRegB, regMatrices[l],
                      regRowMaps[l], regRowImporters[l]);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 3 - scale interface")));

    scaleInterfaceDOFs(regRes, regInterfaceScalings[l], true);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 4 - create coarse vectors")));

    // Transfer to coarse level
    Array<RCP<Vector> > coarseRegX(maxRegPerProc);
    Array<RCP<Vector> > coarseRegB(maxRegPerProc);
    if (useFastMatVec)
    {
      // Get pre-communicated communication patterns for the fast MatVec
      const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = smootherParams[l+1]->get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
      const RCP<Import> regionInterfaceImporter = smootherParams[l+1]->get<RCP<Import>>("Fast MatVec: interface importer");

      for (int j = 0; j < maxRegPerProc; j++) {
        coarseRegX[j] = VectorFactory::Build(regRowMaps[l+1][j], true);
        coarseRegB[j] = VectorFactory::Build(regRowMaps[l+1][j], true);

        ApplyMatVec(SC_ONE, regProlong[l+1][j], regRes[j], SC_ZERO, regionInterfaceImporter,
            regionInterfaceLIDs, coarseRegB[j], Teuchos::TRANS, true);
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegB[j]->getMap()));
      }
    }
    else
    {
      for (int j = 0; j < maxRegPerProc; j++) {
        coarseRegX[j] = VectorFactory::Build(regRowMaps[l+1][j], true);
        coarseRegB[j] = VectorFactory::Build(regRowMaps[l+1][j], true);

        regProlong[l+1][j]->apply(*regRes[j], *coarseRegB[j], Teuchos::TRANS);
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegB[j]->getMap()));
      }

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 5 - sum interface values")));

      sumInterfaceValues(coarseRegB, regRowMaps[l+1], regRowImporters[l+1]);
    }

    tm = Teuchos::null;
    bool coarseZeroInitGuess = true;

    // Call V-cycle recursively
    vCycle(l+1, numLevels,
           coarseRegX, coarseRegB, regMatrices, regProlong, compRowMaps,
           quasiRegRowMaps, regRowMaps, regRowImporters, regInterfaceScalings,
           smootherParams, coarseZeroInitGuess, coarseSolverData, hierarchyData);

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 6 - transfer coarse to fine")));

    // Transfer coarse level correction to fine level
    Array<RCP<Vector> > regCorrection(maxRegPerProc);
    if (useFastMatVec)
    {
      // Get pre-communicated communication patterns for the fast MatVec
      const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = smootherParams[l]->get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
      const RCP<Import> regionInterfaceImporter = smootherParams[l]->get<RCP<Import>>("Fast MatVec: interface importer");

      for (int j = 0; j < maxRegPerProc; j++) {
        regCorrection[j] = VectorFactory::Build(regRowMaps[l][j], true);
        ApplyMatVec(SC_ONE, regProlong[l+1][j], coarseRegX[j], SC_ZERO, regionInterfaceImporter,
            regionInterfaceLIDs, regCorrection[j], Teuchos::NO_TRANS, false);
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegX[j]->getMap()));
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regCorrection[j]->getMap()));
      }
    }
    else{
      for (int j = 0; j < maxRegPerProc; j++) {
        regCorrection[j] = VectorFactory::Build(regRowMaps[l][j], true);
        regProlong[l+1][j]->apply(*coarseRegX[j], *regCorrection[j]);
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegX[j]->getMap()));
        // TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regCorrection[j]->getMap()));
      }
    }

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 7 - add coarse grid correction")));

    // apply coarse grid correction
    for (int j = 0; j < maxRegPerProc; j++) {
      fineRegX[j]->update(SC_ONE, *regCorrection[j], SC_ONE);
    }
    if (coarseZeroInitGuess) zeroInitGuess = true;

//    std::cout << "level: " << l << std::endl;

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 8 - post-smoother")));

    // post-smoothing
    smootherApply(smootherParams[l], fineRegX, fineRegB, regMatrices[l],
                  regRowMaps[l], regRowImporters[l], zeroInitGuess);

    tm = Teuchos::null;

  } else {

    // Coarsest grid solve

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(0);

    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: * - coarsest grid solve")));

    // std::cout << "Applying coarse colver" << std::endl;

    const std::string coarseSolverType = coarseSolverData->get<std::string>("coarse solver type");
    if (coarseSolverType == "smoother") {
      smootherApply(smootherParams[l], fineRegX, fineRegB, regMatrices[l],
                  regRowMaps[l], regRowImporters[l], zeroInitGuess);
    } else {
      zeroInitGuess = false;
      // First get the Xpetra vectors from region to composite format
      RCP<const Map> coarseRowMap = coarseSolverData->get<RCP<const Map> >("compCoarseRowMap");
      RCP<Vector> compX = VectorFactory::Build(coarseRowMap, true);
      RCP<Vector> compRhs = VectorFactory::Build(coarseRowMap, true);
      {
        for (int j = 0; j < maxRegPerProc; j++) {
          RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(regInterfaceScalings[l][j]->getMap());
          inverseInterfaceScaling->reciprocal(*regInterfaceScalings[l][j]);
          fineRegB[j]->elementWiseMultiply(SC_ONE, *fineRegB[j], *inverseInterfaceScaling, SC_ZERO);
        }

        regionalToComposite(fineRegB, compRhs, regRowImporters[l]);
      }

      if (coarseSolverType == "direct")
      {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)

        using DirectCoarseSolver = Amesos2::Solver<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >;
        RCP<DirectCoarseSolver> coarseSolver = coarseSolverData->get<RCP<DirectCoarseSolver> >("direct solver object");

        TEUCHOS_TEST_FOR_EXCEPT_MSG(coarseRowMap->lib()!=Xpetra::UseTpetra,
            "Coarse solver requires Tpetra/Amesos2 stack.");
        TEUCHOS_ASSERT(!coarseSolver.is_null());

        // using Utilities = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

        // From here on we switch to Tpetra for simplicity
        // we could also implement a similar Epetra branch
        using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

      //    *fos << "Attempting to use Amesos2 to solve the coarse grid problem" << std::endl;
        RCP<Tpetra_MultiVector> tX = Utilities::MV2NonConstTpetraMV2(*compX);
        RCP<const Tpetra_MultiVector> tB = Utilities::MV2TpetraMV(compRhs);

        /* Solve!
         *
         * Calling solve() on the coarseSolver should just do a triangular solve, since symbolic
         * and numeric factorization are supposed to have happened during hierarchy setup.
         * Here, we just check if they're done and print message if not.
         *
         * We don't have to change the map of tX and tB since we have configured the Amesos2 solver
         * during its construction to work with non-continuous maps.
         */
        if (not coarseSolver->getStatus().symbolicFactorizationDone())
          *fos << "Symbolic factorization should have been done during hierarchy setup, "
              "but actually is missing. Anyway ... just do it right now." << std::endl;
        if (not coarseSolver->getStatus().numericFactorizationDone())
          *fos << "Numeric factorization should have been done during hierarchy setup, "
              "but actually is missing. Anyway ... just do it right now." << std::endl;
        coarseSolver->solve(tX.ptr(), tB.ptr());
#else
        *fos << "+++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++\n"
             << "+ Coarse level direct solver requires Tpetra and Amesos2.   +\n"
             << "+ Skipping the coarse level solve.                          +\n"
             << "+++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++"
             << std::endl;
#endif
      }
      else if (coarseSolverType == "amg") // use AMG as coarse level solver
      {

        // Extract the hierarchy from the coarseSolverData
        RCP<Hierarchy> amgHierarchy = coarseSolverData->get<RCP<Hierarchy>>("amg hierarchy object");

        // Run a single V-cycle

        amgHierarchy->Iterate(*compRhs, *compX, true, 1);
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(false, "Unknown coarse solver type.");
      }

      // Transform back to region format
      Array<RCP<Vector> > quasiRegX(maxRegPerProc);
      compositeToRegional(compX, quasiRegX, fineRegX,
                          regRowMaps[l],
                          regRowImporters[l]);

      tm = Teuchos::null;
    }
  }

  return;
} // vCycle

#endif // MUELU_SETUPREGIONHIERARCHY_DEF_HPP
