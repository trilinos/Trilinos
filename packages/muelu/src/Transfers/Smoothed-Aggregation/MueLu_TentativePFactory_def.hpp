// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TENTATIVEPFACTORY_DEF_HPP
#define MUELU_TENTATIVEPFACTORY_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_IO.hpp>

#include "Xpetra_TpetraBlockCrsMatrix.hpp"

#include "MueLu_TentativePFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("tentative: build coarse coordinates");
  SET_VALID_ENTRY("tentative: constant column sums");
#undef SET_VALID_ENTRY
  validParamList->set<std::string>("Nullspace name", "Nullspace", "Name for the input nullspace");

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Generating factory of the aggregates");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<RCP<const FactoryBase> >("Scaled Nullspace", Teuchos::null, "Generating factory of the scaled nullspace");
  validParamList->set<RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase> >("CoarseMap", Teuchos::null, "Generating factory of the coarse map");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Generating factory of the coordinates");
  validParamList->set<RCP<const FactoryBase> >("Node Comm", Teuchos::null, "Generating factory of the node level communicator");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  const ParameterList& pL = GetParameterList();
  // NOTE: This guy can only either be 'Nullspace' or 'Scaled Nullspace' or else the validator above will cause issues
  std::string nspName = "Nullspace";
  if (pL.isParameter("Nullspace name")) nspName = pL.get<std::string>("Nullspace name");

  Input(fineLevel, "A");
  Input(fineLevel, "Aggregates");
  Input(fineLevel, nspName);
  Input(fineLevel, "UnAmalgamationInfo");
  Input(fineLevel, "CoarseMap");
  if (fineLevel.GetLevelID() == 0 &&
      fineLevel.IsAvailable("Coordinates", NoFactory::get()) &&  // we have coordinates (provided by user app)
      pL.get<bool>("tentative: build coarse coordinates")) {     // and we want coordinates on other levels
    bTransferCoordinates_ = true;                                // then set the transfer coordinates flag to true
    Input(fineLevel, "Coordinates");
  } else if (bTransferCoordinates_) {
    Input(fineLevel, "Coordinates");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);
  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinate_type;
  typedef Xpetra::MultiVector<coordinate_type, LO, GO, NO> RealValuedMultiVector;
  typedef Xpetra::MultiVectorFactory<coordinate_type, LO, GO, NO> RealValuedMultiVectorFactory;

  const ParameterList& pL = GetParameterList();
  std::string nspName     = "Nullspace";
  if (pL.isParameter("Nullspace name")) nspName = pL.get<std::string>("Nullspace name");

  RCP<Matrix> Ptentative;
  RCP<Matrix> A              = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<Aggregates> aggregates = Get<RCP<Aggregates> >(fineLevel, "Aggregates");
  // No coarse DoFs so we need to bail by setting Ptentattive to null and returning
  //  This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
  if (aggregates->GetNumGlobalAggregatesComputeIfNeeded() == 0) {
    Ptentative = Teuchos::null;
    Set(coarseLevel, "P", Ptentative);
    return;
  }

  RCP<AmalgamationInfo> amalgInfo = Get<RCP<AmalgamationInfo> >(fineLevel, "UnAmalgamationInfo");
  RCP<MultiVector> fineNullspace  = Get<RCP<MultiVector> >(fineLevel, nspName);
  RCP<const Map> coarseMap        = Get<RCP<const Map> >(fineLevel, "CoarseMap");
  RCP<RealValuedMultiVector> fineCoords;
  if (bTransferCoordinates_) {
    fineCoords = Get<RCP<RealValuedMultiVector> >(fineLevel, "Coordinates");
  }

  // FIXME: We should remove the NodeComm on levels past the threshold
  if (fineLevel.IsAvailable("Node Comm")) {
    RCP<const Teuchos::Comm<int> > nodeComm = Get<RCP<const Teuchos::Comm<int> > >(fineLevel, "Node Comm");
    Set<RCP<const Teuchos::Comm<int> > >(coarseLevel, "Node Comm", nodeComm);
  }

  // NOTE:  We check DomainMap here rather than RowMap because those are different for BlockCrs matrices
  TEUCHOS_TEST_FOR_EXCEPTION(A->getDomainMap()->getLocalNumElements() != fineNullspace->getMap()->getLocalNumElements(),
                             Exceptions::RuntimeError, "MueLu::TentativePFactory::MakeTentative: Size mismatch between A and Nullspace");

  RCP<MultiVector> coarseNullspace;
  RCP<RealValuedMultiVector> coarseCoords;

  if (bTransferCoordinates_) {
    //*** Create the coarse coordinates ***
    // First create the coarse map and coarse multivector
    ArrayView<const GO> elementAList = coarseMap->getLocalElementList();
    LO blkSize                       = 1;
    if (rcp_dynamic_cast<const StridedMap>(coarseMap) != Teuchos::null) {
      blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap)->getFixedBlockSize();
    }
    GO indexBase      = coarseMap->getIndexBase();
    LO numCoarseNodes = Teuchos::as<LO>(elementAList.size() / blkSize);
    Array<GO> nodeList(numCoarseNodes);
    const int numDimensions = fineCoords->getNumVectors();

    for (LO i = 0; i < numCoarseNodes; i++) {
      nodeList[i] = (elementAList[i * blkSize] - indexBase) / blkSize + indexBase;
    }
    RCP<const Map> coarseCoordsMap = MapFactory::Build(fineCoords->getMap()->lib(),
                                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                       nodeList,
                                                       indexBase,
                                                       fineCoords->getMap()->getComm());
    coarseCoords                   = RealValuedMultiVectorFactory::Build(coarseCoordsMap, numDimensions);

    // Create overlapped fine coordinates to reduce global communication
    RCP<RealValuedMultiVector> ghostedCoords;
    if (aggregates->AggregatesCrossProcessors()) {
      RCP<const Map> aggMap      = aggregates->GetMap();
      RCP<const Import> importer = ImportFactory::Build(fineCoords->getMap(), aggMap);

      ghostedCoords = RealValuedMultiVectorFactory::Build(aggMap, numDimensions);
      ghostedCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);
    } else {
      ghostedCoords = fineCoords;
    }

    // Get some info about aggregates
    int myPID                             = coarseCoordsMap->getComm()->getRank();
    LO numAggs                            = aggregates->GetNumAggregates();
    ArrayRCP<LO> aggSizes                 = aggregates->ComputeAggregateSizesArrayRCP();
    const ArrayRCP<const LO> vertex2AggID = aggregates->GetVertex2AggId()->getData(0);
    const ArrayRCP<const LO> procWinner   = aggregates->GetProcWinner()->getData(0);

    // Fill in coarse coordinates
    for (int dim = 0; dim < numDimensions; ++dim) {
      ArrayRCP<const coordinate_type> fineCoordsData = ghostedCoords->getData(dim);
      ArrayRCP<coordinate_type> coarseCoordsData     = coarseCoords->getDataNonConst(dim);

      for (LO lnode = 0; lnode < Teuchos::as<LO>(vertex2AggID.size()); lnode++) {
        if (procWinner[lnode] == myPID &&
            lnode < fineCoordsData.size() &&
            vertex2AggID[lnode] < coarseCoordsData.size() &&
            Teuchos::ScalarTraits<coordinate_type>::isnaninf(fineCoordsData[lnode]) == false) {
          coarseCoordsData[vertex2AggID[lnode]] += fineCoordsData[lnode];
        }
      }
      for (LO agg = 0; agg < numAggs; agg++) {
        coarseCoordsData[agg] /= aggSizes[agg];
      }
    }
  }

  if (!aggregates->AggregatesCrossProcessors()) {
    if (Xpetra::Helpers<SC, LO, GO, NO>::isTpetraBlockCrs(A)) {
      BuildPuncoupledBlockCrs(A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace, coarseLevel.GetLevelID());
    } else {
      BuildPuncoupled(A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace, coarseLevel.GetLevelID());
    }
  } else
    BuildPcoupled(A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace);

  // If available, use striding information of fine level matrix A for range
  // map and coarseMap as domain map; otherwise use plain range map of
  // Ptent = plain range map of A for range map and coarseMap as domain map.
  // NOTE:
  // The latter is not really safe, since there is no striding information
  // for the range map. This is not really a problem, since striding
  // information is always available on the intermedium levels and the
  // coarsest levels.
  if (A->IsView("stridedMaps") == true)
    Ptentative->CreateView("stridedMaps", A->getRowMap("stridedMaps"), coarseMap);
  else
    Ptentative->CreateView("stridedMaps", Ptentative->getRangeMap(), coarseMap);

  if (bTransferCoordinates_) {
    Set(coarseLevel, "Coordinates", coarseCoords);
  }
  Set(coarseLevel, "Nullspace", coarseNullspace);
  Set(coarseLevel, "P", Ptentative);

  if (IsPrint(Statistics2)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Ptentative, "Ptent", params);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildPuncoupledBlockCrs(RCP<Matrix> A, RCP<Aggregates> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                            RCP<const Map> coarsePointMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace, const int levelID) const {
  /* This routine generates a BlockCrs P for a BlockCrs A.  There are a few assumptions here, which meet the use cases we care about, but could
     be generalized later, if we ever need to do so:
     1) Null space dimension === block size of matrix:  So no elasticity right now
     2) QR is not supported:  Under assumption #1, this shouldn't cause problems.
     3) Maps are "good": Aka the first chunk of the ColMap is the RowMap.

     These assumptions keep our code way simpler and still support the use cases we actually care about.
   */

  RCP<const Map> rowMap   = A->getRowMap();
  RCP<const Map> rangeMap = A->getRangeMap();
  RCP<const Map> colMap   = A->getColMap();
  //    const size_t numFinePointRows = rangeMap->getLocalNumElements();
  const size_t numFineBlockRows = rowMap->getLocalNumElements();

  typedef Teuchos::ScalarTraits<SC> STS;
  // typedef typename STS::magnitudeType Magnitude;
  const SC zero    = STS::zero();
  const SC one     = STS::one();
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  const GO numAggs      = aggregates->GetNumAggregates();
  const size_t NSDim    = fineNullspace->getNumVectors();
  ArrayRCP<LO> aggSizes = aggregates->ComputeAggregateSizesArrayRCP();

  // Need to generate the coarse block map
  // NOTE: We assume NSDim == block size here
  // NOTE: We also assume that coarseMap has contiguous GIDs
  // const size_t numCoarsePointRows = coarsePointMap->getLocalNumElements();
  const size_t numCoarseBlockRows = coarsePointMap->getLocalNumElements() / NSDim;
  RCP<const Map> coarseBlockMap   = MapFactory::Build(coarsePointMap->lib(),
                                                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                      numCoarseBlockRows,
                                                      coarsePointMap->getIndexBase(),
                                                      coarsePointMap->getComm());
  // Sanity checking
  const ParameterList& pL     = GetParameterList();
  const bool& doQRStep        = pL.get<bool>("tentative: calculate qr");
  const bool& constantColSums = pL.get<bool>("tentative: constant column sums");

  TEUCHOS_TEST_FOR_EXCEPTION(doQRStep && constantColSums, Exceptions::RuntimeError,
                             "MueLu::TentativePFactory::MakeTentative: cannot use 'constant column sums' and 'calculate qr' at the same time");

  // The aggregates use the amalgamated column map, which in this case is what we want

  // Aggregates map is based on the amalgamated column map
  // We can skip global-to-local conversion if LIDs in row map are
  // same as LIDs in column map
  bool goodMap = MueLu::Utilities<SC, LO, GO, NO>::MapsAreNested(*rowMap, *colMap);

  // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
  // aggStart is a pointer into aggToRowMapLO
  // aggStart[i]..aggStart[i+1] are indices into aggToRowMapLO
  // aggToRowMapLO[aggStart[i]]..aggToRowMapLO[aggStart[i+1]-1] are the DOFs in aggregate i
  ArrayRCP<LO> aggStart;
  ArrayRCP<LO> aggToRowMapLO;
  ArrayRCP<GO> aggToRowMapGO;
  if (goodMap) {
    amalgInfo->UnamalgamateAggregatesLO(*aggregates, aggStart, aggToRowMapLO);
    GetOStream(Runtime1) << "Column map is consistent with the row map, good." << std::endl;
  } else {
    throw std::runtime_error("TentativePFactory::PuncoupledBlockCrs: Inconsistent maps not currently supported");
  }

  coarseNullspace = MultiVectorFactory::Build(coarsePointMap, NSDim);

  // Pull out the nullspace vectors so that we can have random access.
  ArrayRCP<ArrayRCP<const SC> > fineNS(NSDim);
  ArrayRCP<ArrayRCP<SC> > coarseNS(NSDim);
  for (size_t i = 0; i < NSDim; i++) {
    fineNS[i] = fineNullspace->getData(i);
    if (coarsePointMap->getLocalNumElements() > 0)
      coarseNS[i] = coarseNullspace->getDataNonConst(i);
  }

  // BlockCrs requires that we build the (block) graph first, so let's do that...
  // NOTE: Because we're assuming that the NSDim == BlockSize, we only have one
  // block non-zero per row in the matrix;
  RCP<CrsGraph> BlockGraph = CrsGraphFactory::Build(rowMap, coarseBlockMap, 0);
  ArrayRCP<size_t> iaPtent;
  ArrayRCP<LO> jaPtent;
  BlockGraph->allocateAllIndices(numFineBlockRows, iaPtent, jaPtent);
  ArrayView<size_t> ia = iaPtent();
  ArrayView<LO> ja     = jaPtent();

  for (size_t i = 0; i < numFineBlockRows; i++) {
    ia[i] = i;
    ja[i] = INVALID;
  }
  ia[numCoarseBlockRows] = numCoarseBlockRows;

  for (GO agg = 0; agg < numAggs; agg++) {
    LO aggSize                   = aggStart[agg + 1] - aggStart[agg];
    Xpetra::global_size_t offset = agg;

    for (LO j = 0; j < aggSize; j++) {
      // FIXME: Allow for bad maps
      const LO localRow     = aggToRowMapLO[aggStart[agg] + j];
      const size_t rowStart = ia[localRow];
      ja[rowStart]          = offset;
    }
  }

  // Compress storage (remove all INVALID, which happen when we skip zeros)
  // We do that in-place
  size_t ia_tmp = 0, nnz = 0;
  for (size_t i = 0; i < numFineBlockRows; i++) {
    for (size_t j = ia_tmp; j < ia[i + 1]; j++)
      if (ja[j] != INVALID) {
        ja[nnz] = ja[j];
        nnz++;
      }
    ia_tmp    = ia[i + 1];
    ia[i + 1] = nnz;
  }

  if (rowMap->lib() == Xpetra::UseTpetra) {
    // - Cannot resize for Epetra, as it checks for same pointers
    // - Need to resize for Tpetra, as it check ().size() == ia[numRows]
    // NOTE: these invalidate ja and val views
    jaPtent.resize(nnz);
  }

  GetOStream(Runtime1) << "TentativePFactory : generating block graph" << std::endl;
  BlockGraph->setAllIndices(iaPtent, jaPtent);

  // Managing labels & constants for ESFC
  {
    RCP<ParameterList> FCparams;
    if (pL.isSublist("matrixmatrix: kernel params"))
      FCparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
    else
      FCparams = rcp(new ParameterList);
    // By default, we don't need global constants for TentativeP, but we do want it for the graph
    // if we're printing statistics, so let's leave it on for now.
    FCparams->set("compute global constants", FCparams->get("compute global constants", true));
    std::string levelIDs = toString(levelID);
    FCparams->set("Timer Label", std::string("MueLu::TentativeP-") + levelIDs);
    RCP<const Export> dummy_e;
    RCP<const Import> dummy_i;
    BlockGraph->expertStaticFillComplete(coarseBlockMap, rowMap, dummy_i, dummy_e, FCparams);
  }

  // Now let's make a BlockCrs Matrix
  // NOTE: Assumes block size== NSDim
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO> > P_xpetra            = Xpetra::CrsMatrixFactory<SC, LO, GO, NO>::BuildBlock(BlockGraph, coarsePointMap, rangeMap, NSDim);
  RCP<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO> > P_tpetra = rcp_dynamic_cast<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO> >(P_xpetra);
  if (P_tpetra.is_null()) throw std::runtime_error("BuildPUncoupled: Matrix factory did not return a Tpetra::BlockCrsMatrix");
  RCP<CrsMatrixWrap> P_wrap = rcp(new CrsMatrixWrap(P_xpetra));

  /////////////////////////////
  //      "no-QR" option     //
  /////////////////////////////
  // Local Q factor is just the fine nullspace support over the current aggregate.
  // Local R factor is the identity.
  // NOTE: We're not going to do a QR here as we're assuming that blocksize == NSDim
  // NOTE: "goodMap" case only
  Teuchos::Array<Scalar> block(NSDim * NSDim, zero);
  Teuchos::Array<LO> bcol(1);

  GetOStream(Runtime1) << "TentativePFactory : bypassing local QR phase" << std::endl;
  for (LO agg = 0; agg < numAggs; agg++) {
    bcol[0]                      = agg;
    const LO aggSize             = aggStart[agg + 1] - aggStart[agg];
    Xpetra::global_size_t offset = agg * NSDim;

    // Process each row in the local Q factor
    // NOTE: Blocks are in row-major order
    for (LO j = 0; j < aggSize; j++) {
      const LO localBlockRow = aggToRowMapLO[aggStart[agg] + j];

      for (size_t r = 0; r < NSDim; r++) {
        LO localPointRow = localBlockRow * NSDim + r;
        for (size_t c = 0; c < NSDim; c++)
          block[r * NSDim + c] = fineNS[c][localPointRow];
      }
      // NOTE: Assumes columns==aggs and are ordered sequentially
      P_tpetra->replaceLocalValues(localBlockRow, bcol(), block());

    }  // end aggSize

    for (size_t j = 0; j < NSDim; j++)
      coarseNS[j][offset + j] = one;

  }  // for (GO agg = 0; agg < numAggs; agg++)

  Ptentative = P_wrap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildPcoupled(RCP<Matrix> A, RCP<Aggregates> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                  RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const {
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType Magnitude;
  const SC zero = STS::zero();
  const SC one  = STS::one();

  // number of aggregates
  GO numAggs = aggregates->GetNumAggregates();

  // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
  // aggStart is a pointer into aggToRowMap
  // aggStart[i]..aggStart[i+1] are indices into aggToRowMap
  // aggToRowMap[aggStart[i]]..aggToRowMap[aggStart[i+1]-1] are the DOFs in aggregate i
  ArrayRCP<LO> aggStart;
  ArrayRCP<GO> aggToRowMap;
  amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);

  // find size of the largest aggregate
  LO maxAggSize = 0;
  for (GO i = 0; i < numAggs; ++i) {
    LO sizeOfThisAgg = aggStart[i + 1] - aggStart[i];
    if (sizeOfThisAgg > maxAggSize) maxAggSize = sizeOfThisAgg;
  }

  // dimension of fine level nullspace
  const size_t NSDim = fineNullspace->getNumVectors();

  // index base for coarse Dof map (usually 0)
  GO indexBase = A->getRowMap()->getIndexBase();

  const RCP<const Map> nonUniqueMap         = amalgInfo->ComputeUnamalgamatedImportDofMap(*aggregates);
  const RCP<const Map> uniqueMap            = A->getDomainMap();
  RCP<const Import> importer                = ImportFactory::Build(uniqueMap, nonUniqueMap);
  RCP<MultiVector> fineNullspaceWithOverlap = MultiVectorFactory::Build(nonUniqueMap, NSDim);
  fineNullspaceWithOverlap->doImport(*fineNullspace, *importer, Xpetra::INSERT);

  // Pull out the nullspace vectors so that we can have random access.
  ArrayRCP<ArrayRCP<const SC> > fineNS(NSDim);
  for (size_t i = 0; i < NSDim; ++i)
    fineNS[i] = fineNullspaceWithOverlap->getData(i);

  // Allocate storage for the coarse nullspace.
  coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

  ArrayRCP<ArrayRCP<SC> > coarseNS(NSDim);
  for (size_t i = 0; i < NSDim; ++i)
    if (coarseMap->getLocalNumElements() > 0) coarseNS[i] = coarseNullspace->getDataNonConst(i);

  // This makes the rowmap of Ptent the same as that of A->
  // This requires moving some parts of some local Q's to other processors
  // because aggregates can span processors.
  RCP<const Map> rowMapForPtent = A->getRowMap();
  const Map& rowMapForPtentRef  = *rowMapForPtent;

  // Set up storage for the rows of the local Qs that belong to other processors.
  // FIXME This is inefficient and could be done within the main loop below with std::vector's.
  RCP<const Map> colMap = A->getColMap();

  RCP<const Map> ghostQMap;
  RCP<MultiVector> ghostQvalues;
  Array<RCP<Xpetra::Vector<GO, LO, GO, Node> > > ghostQcolumns;
  RCP<Xpetra::Vector<GO, LO, GO, Node> > ghostQrowNums;
  ArrayRCP<ArrayRCP<SC> > ghostQvals;
  ArrayRCP<ArrayRCP<GO> > ghostQcols;
  ArrayRCP<GO> ghostQrows;

  Array<GO> ghostGIDs;
  for (LO j = 0; j < numAggs; ++j) {
    for (LO k = aggStart[j]; k < aggStart[j + 1]; ++k) {
      if (rowMapForPtentRef.isNodeGlobalElement(aggToRowMap[k]) == false) {
        ghostGIDs.push_back(aggToRowMap[k]);
      }
    }
  }
  ghostQMap = MapFactory::Build(A->getRowMap()->lib(),
                                Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                ghostGIDs,
                                indexBase, A->getRowMap()->getComm());  // JG:Xpetra::global_size_t>?
  // Vector to hold bits of Q that go to other processors.
  ghostQvalues = MultiVectorFactory::Build(ghostQMap, NSDim);
  // Note that Epetra does not support MultiVectors templated on Scalar != double.
  // So to work around this, we allocate an array of Vectors.  This shouldn't be too
  // expensive, as the number of Vectors is NSDim.
  ghostQcolumns.resize(NSDim);
  for (size_t i = 0; i < NSDim; ++i)
    ghostQcolumns[i] = Xpetra::VectorFactory<GO, LO, GO, Node>::Build(ghostQMap);
  ghostQrowNums = Xpetra::VectorFactory<GO, LO, GO, Node>::Build(ghostQMap);
  if (ghostQvalues->getLocalLength() > 0) {
    ghostQvals.resize(NSDim);
    ghostQcols.resize(NSDim);
    for (size_t i = 0; i < NSDim; ++i) {
      ghostQvals[i] = ghostQvalues->getDataNonConst(i);
      ghostQcols[i] = ghostQcolumns[i]->getDataNonConst(0);
    }
    ghostQrows = ghostQrowNums->getDataNonConst(0);
  }

  // importer to handle moving Q
  importer = ImportFactory::Build(ghostQMap, A->getRowMap());

  // Dense QR solver
  Teuchos::SerialQRDenseSolver<LO, SC> qrSolver;

  // Allocate temporary storage for the tentative prolongator.
  Array<GO> globalColPtr(maxAggSize * NSDim, 0);
  Array<LO> localColPtr(maxAggSize * NSDim, 0);
  Array<SC> valPtr(maxAggSize * NSDim, 0.);

  // Create column map for Ptent, estimate local #nonzeros in Ptent,  and create Ptent itself.
  const Map& coarseMapRef = *coarseMap;

  // For the 3-arrays constructor
  ArrayRCP<size_t> ptent_rowptr;
  ArrayRCP<LO> ptent_colind;
  ArrayRCP<Scalar> ptent_values;

  // Because ArrayRCPs are slow...
  ArrayView<size_t> rowptr_v;
  ArrayView<LO> colind_v;
  ArrayView<Scalar> values_v;

  // For temporary usage
  Array<size_t> rowptr_temp;
  Array<LO> colind_temp;
  Array<Scalar> values_temp;

  RCP<CrsMatrix> PtentCrs;

  RCP<CrsMatrixWrap> PtentCrsWrap = rcp(new CrsMatrixWrap(rowMapForPtent, NSDim));
  PtentCrs                        = PtentCrsWrap->getCrsMatrix();
  Ptentative                      = PtentCrsWrap;

  //*****************************************************************
  // Loop over all aggregates and calculate local QR decompositions.
  //*****************************************************************
  GO qctr                    = 0;  // for indexing into Ptent data vectors
  const Map& nonUniqueMapRef = *nonUniqueMap;

  size_t total_nnz_count = 0;

  for (GO agg = 0; agg < numAggs; ++agg) {
    LO myAggSize = aggStart[agg + 1] - aggStart[agg];
    // For each aggregate, extract the corresponding piece of the nullspace and put it in the flat array,
    // "localQR" (in column major format) for the QR routine.
    Teuchos::SerialDenseMatrix<LO, SC> localQR(myAggSize, NSDim);
    for (size_t j = 0; j < NSDim; ++j) {
      bool bIsZeroNSColumn = true;
      for (LO k = 0; k < myAggSize; ++k) {
        // aggToRowMap[aggPtr[i]+k] is the kth DOF in the ith aggregate
        // fineNS[j][n] is the nth entry in the jth NS vector
        try {
          SC nsVal      = fineNS[j][nonUniqueMapRef.getLocalElement(aggToRowMap[aggStart[agg] + k])];  // extract information from fine level NS
          localQR(k, j) = nsVal;
          if (nsVal != zero) bIsZeroNSColumn = false;
        } catch (...) {
          GetOStream(Runtime1, -1) << "length of fine level nsp: " << fineNullspace->getGlobalLength() << std::endl;
          GetOStream(Runtime1, -1) << "length of fine level nsp w overlap: " << fineNullspaceWithOverlap->getGlobalLength() << std::endl;
          GetOStream(Runtime1, -1) << "(local?) aggId=" << agg << std::endl;
          GetOStream(Runtime1, -1) << "aggSize=" << myAggSize << std::endl;
          GetOStream(Runtime1, -1) << "agg DOF=" << k << std::endl;
          GetOStream(Runtime1, -1) << "NS vector j=" << j << std::endl;
          GetOStream(Runtime1, -1) << "j*myAggSize + k = " << j * myAggSize + k << std::endl;
          GetOStream(Runtime1, -1) << "aggToRowMap[" << agg << "][" << k << "] = " << aggToRowMap[aggStart[agg] + k] << std::endl;
          GetOStream(Runtime1, -1) << "id aggToRowMap[agg][k]=" << aggToRowMap[aggStart[agg] + k] << " is global element in nonUniqueMap = " << nonUniqueMapRef.isNodeGlobalElement(aggToRowMap[aggStart[agg] + k]) << std::endl;
          GetOStream(Runtime1, -1) << "colMap local id aggToRowMap[agg][k]=" << nonUniqueMapRef.getLocalElement(aggToRowMap[aggStart[agg] + k]) << std::endl;
          GetOStream(Runtime1, -1) << "fineNS...=" << fineNS[j][nonUniqueMapRef.getLocalElement(aggToRowMap[aggStart[agg] + k])] << std::endl;
          GetOStream(Errors, -1) << "caught an error!" << std::endl;
        }
      }  // for (LO k=0 ...
      TEUCHOS_TEST_FOR_EXCEPTION(bIsZeroNSColumn == true, Exceptions::RuntimeError, "MueLu::TentativePFactory::MakeTentative: fine level NS part has a zero column. Error.");
    }  // for (LO j=0 ...

    Xpetra::global_size_t offset = agg * NSDim;

    if (myAggSize >= Teuchos::as<LocalOrdinal>(NSDim)) {
      // calculate QR decomposition (standard)
      // R is stored in localQR (size: myAggSize x NSDim)

      // Householder multiplier
      SC tau = localQR(0, 0);

      if (NSDim == 1) {
        // Only one nullspace vector, so normalize by hand
        Magnitude dtemp = 0;
        for (size_t k = 0; k < Teuchos::as<size_t>(myAggSize); ++k) {
          Magnitude tmag = STS::magnitude(localQR(k, 0));
          dtemp += tmag * tmag;
        }
        dtemp         = Teuchos::ScalarTraits<Magnitude>::squareroot(dtemp);
        tau           = localQR(0, 0);
        localQR(0, 0) = dtemp;
      } else {
        qrSolver.setMatrix(Teuchos::rcp(&localQR, false));
        qrSolver.factor();
      }

      // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
      // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
      // This stores the (offset+k)th entry only if it is local according to the coarseMap.
      for (size_t j = 0; j < NSDim; ++j) {
        for (size_t k = 0; k <= j; ++k) {
          try {
            if (coarseMapRef.isNodeLocalElement(offset + k)) {
              coarseNS[j][offset + k] = localQR(k, j);  // TODO is offset+k the correct local ID?!
            }
          } catch (...) {
            GetOStream(Errors, -1) << "caught error in coarseNS insert, j=" << j << ", offset+k = " << offset + k << std::endl;
          }
        }
      }

      // Calculate Q, the tentative prolongator.
      // The Lapack GEQRF call only works for myAggsize >= NSDim

      if (NSDim == 1) {
        // Only one nullspace vector, so calculate Q by hand
        Magnitude dtemp = Teuchos::ScalarTraits<SC>::magnitude(localQR(0, 0));
        localQR(0, 0)   = tau;
        dtemp           = 1 / dtemp;
        for (LocalOrdinal i = 0; i < myAggSize; ++i) {
          localQR(i, 0) *= dtemp;
        }
      } else {
        qrSolver.formQ();
        Teuchos::RCP<Teuchos::SerialDenseMatrix<LO, SC> > qFactor = qrSolver.getQ();
        for (size_t j = 0; j < NSDim; j++) {
          for (size_t i = 0; i < Teuchos::as<size_t>(myAggSize); i++) {
            localQR(i, j) = (*qFactor)(i, j);
          }
        }
      }

      // end default case (myAggSize >= NSDim)
    } else {  // special handling for myAggSize < NSDim (i.e. 1pt nodes)
      // See comments for the uncoupled case

      // R = extended (by adding identity rows) localQR
      for (size_t j = 0; j < NSDim; j++)
        for (size_t k = 0; k < NSDim; k++) {
          TEUCHOS_TEST_FOR_EXCEPTION(!coarseMapRef.isNodeLocalElement(offset + k), Exceptions::RuntimeError,
                                     "Caught error in coarseNS insert, j=" << j << ", offset+k = " << offset + k);

          if (k < as<size_t>(myAggSize))
            coarseNS[j][offset + k] = localQR(k, j);
          else
            coarseNS[j][offset + k] = (k == j ? one : zero);
        }

      // Q = I (rectangular)
      for (size_t i = 0; i < as<size_t>(myAggSize); i++)
        for (size_t j = 0; j < NSDim; j++)
          localQR(i, j) = (j == i ? one : zero);
    }  // end else (special handling for 1pt aggregates)

    // Process each row in the local Q factor.  If the row is local to the current processor
    // according to the rowmap, insert it into Ptentative.  Otherwise, save it in ghostQ
    // to be communicated later to the owning processor.
    // FIXME -- what happens if maps are blocked?
    for (GO j = 0; j < myAggSize; ++j) {
      // This loop checks whether row associated with current DOF is local, according to rowMapForPtent.
      // If it is, the row is inserted.  If not, the row number, columns, and values are saved in
      // MultiVectors that will be sent to other processors.
      GO globalRow = aggToRowMap[aggStart[agg] + j];

      // TODO is the use of Xpetra::global_size_t below correct?
      if (rowMapForPtentRef.isNodeGlobalElement(globalRow) == false) {
        ghostQrows[qctr] = globalRow;
        for (size_t k = 0; k < NSDim; ++k) {
          ghostQcols[k][qctr] = coarseMapRef.getGlobalElement(agg * NSDim + k);
          ghostQvals[k][qctr] = localQR(j, k);
        }
        ++qctr;
      } else {
        size_t nnz = 0;
        for (size_t k = 0; k < NSDim; ++k) {
          try {
            if (localQR(j, k) != Teuchos::ScalarTraits<SC>::zero()) {
              localColPtr[nnz]  = agg * NSDim + k;
              globalColPtr[nnz] = coarseMapRef.getGlobalElement(localColPtr[nnz]);
              valPtr[nnz]       = localQR(j, k);
              ++total_nnz_count;
              ++nnz;
            }
          } catch (...) {
            GetOStream(Errors, -1) << "caught error in colPtr/valPtr insert, current index=" << nnz << std::endl;
          }
        }  // for (size_t k=0; k<NSDim; ++k)

        try {
          Ptentative->insertGlobalValues(globalRow, globalColPtr.view(0, nnz), valPtr.view(0, nnz));
        } catch (...) {
          GetOStream(Errors, -1) << "pid " << A->getRowMap()->getComm()->getRank()
                                 << "caught error during Ptent row insertion, global row "
                                 << globalRow << std::endl;
        }
      }
    }  // for (GO j=0; j<myAggSize; ++j)

  }  // for (LO agg=0; agg<numAggs; ++agg)

  // ***********************************************************
  // ************* end of aggregate-wise QR ********************
  // ***********************************************************
  GetOStream(Runtime1) << "TentativePFactory : aggregates may cross process boundaries" << std::endl;
  // Import ghost parts of Q factors and insert into Ptentative.
  // First import just the global row numbers.
  RCP<Xpetra::Vector<GO, LO, GO, Node> > targetQrowNums = Xpetra::VectorFactory<GO, LO, GO, Node>::Build(rowMapForPtent);
  targetQrowNums->putScalar(-1);
  targetQrowNums->doImport(*ghostQrowNums, *importer, Xpetra::INSERT);
  ArrayRCP<GO> targetQrows = targetQrowNums->getDataNonConst(0);

  // Now create map based on just the row numbers imported.
  Array<GO> gidsToImport;
  gidsToImport.reserve(targetQrows.size());
  for (typename ArrayRCP<GO>::iterator r = targetQrows.begin(); r != targetQrows.end(); ++r) {
    if (*r > -1) {
      gidsToImport.push_back(*r);
    }
  }
  RCP<const Map> reducedMap = MapFactory::Build(A->getRowMap()->lib(),
                                                Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                gidsToImport, indexBase, A->getRowMap()->getComm());

  // Import using the row numbers that this processor will receive.
  importer = ImportFactory::Build(ghostQMap, reducedMap);

  Array<RCP<Xpetra::Vector<GO, LO, GO, Node> > > targetQcolumns(NSDim);
  for (size_t i = 0; i < NSDim; ++i) {
    targetQcolumns[i] = Xpetra::VectorFactory<GO, LO, GO, Node>::Build(reducedMap);
    targetQcolumns[i]->doImport(*(ghostQcolumns[i]), *importer, Xpetra::INSERT);
  }
  RCP<MultiVector> targetQvalues = MultiVectorFactory::Build(reducedMap, NSDim);
  targetQvalues->doImport(*ghostQvalues, *importer, Xpetra::INSERT);

  ArrayRCP<ArrayRCP<SC> > targetQvals;
  ArrayRCP<ArrayRCP<GO> > targetQcols;
  if (targetQvalues->getLocalLength() > 0) {
    targetQvals.resize(NSDim);
    targetQcols.resize(NSDim);
    for (size_t i = 0; i < NSDim; ++i) {
      targetQvals[i] = targetQvalues->getDataNonConst(i);
      targetQcols[i] = targetQcolumns[i]->getDataNonConst(0);
    }
  }

  valPtr       = Array<SC>(NSDim, 0.);
  globalColPtr = Array<GO>(NSDim, 0);
  for (typename Array<GO>::iterator r = gidsToImport.begin(); r != gidsToImport.end(); ++r) {
    if (targetQvalues->getLocalLength() > 0) {
      for (size_t j = 0; j < NSDim; ++j) {
        valPtr[j]       = targetQvals[j][reducedMap->getLocalElement(*r)];
        globalColPtr[j] = targetQcols[j][reducedMap->getLocalElement(*r)];
      }
      Ptentative->insertGlobalValues(*r, globalColPtr.view(0, NSDim), valPtr.view(0, NSDim));
    }  // if (targetQvalues->getLocalLength() > 0)
  }

  Ptentative->fillComplete(coarseMap, A->getDomainMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildPuncoupled(RCP<Matrix> A, RCP<Aggregates> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                    RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace, const int levelID) const {
  RCP<const Map> rowMap = A->getRowMap();
  RCP<const Map> colMap = A->getColMap();
  const size_t numRows  = rowMap->getLocalNumElements();

  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType Magnitude;
  const SC zero    = STS::zero();
  const SC one     = STS::one();
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  const GO numAggs      = aggregates->GetNumAggregates();
  const size_t NSDim    = fineNullspace->getNumVectors();
  ArrayRCP<LO> aggSizes = aggregates->ComputeAggregateSizesArrayRCP();

  // Sanity checking
  const ParameterList& pL     = GetParameterList();
  const bool& doQRStep        = pL.get<bool>("tentative: calculate qr");
  const bool& constantColSums = pL.get<bool>("tentative: constant column sums");

  TEUCHOS_TEST_FOR_EXCEPTION(doQRStep && constantColSums, Exceptions::RuntimeError,
                             "MueLu::TentativePFactory::MakeTentative: cannot use 'constant column sums' and 'calculate qr' at the same time");

  // Aggregates map is based on the amalgamated column map
  // We can skip global-to-local conversion if LIDs in row map are
  // same as LIDs in column map
  bool goodMap = MueLu::Utilities<SC, LO, GO, NO>::MapsAreNested(*rowMap, *colMap);

  // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
  // aggStart is a pointer into aggToRowMapLO
  // aggStart[i]..aggStart[i+1] are indices into aggToRowMapLO
  // aggToRowMapLO[aggStart[i]]..aggToRowMapLO[aggStart[i+1]-1] are the DOFs in aggregate i
  ArrayRCP<LO> aggStart;
  ArrayRCP<LO> aggToRowMapLO;
  ArrayRCP<GO> aggToRowMapGO;
  if (goodMap) {
    amalgInfo->UnamalgamateAggregatesLO(*aggregates, aggStart, aggToRowMapLO);
    GetOStream(Runtime1) << "Column map is consistent with the row map, good." << std::endl;

  } else {
    amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMapGO);
    GetOStream(Warnings0) << "Column map is not consistent with the row map\n"
                          << "using GO->LO conversion with performance penalty" << std::endl;
  }
  coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

  // Pull out the nullspace vectors so that we can have random access.
  ArrayRCP<ArrayRCP<const SC> > fineNS(NSDim);
  ArrayRCP<ArrayRCP<SC> > coarseNS(NSDim);
  for (size_t i = 0; i < NSDim; i++) {
    fineNS[i] = fineNullspace->getData(i);
    if (coarseMap->getLocalNumElements() > 0)
      coarseNS[i] = coarseNullspace->getDataNonConst(i);
  }

  size_t nnzEstimate = numRows * NSDim;

  // Time to construct the matrix and fill in the values
  Ptentative              = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0));
  RCP<CrsMatrix> PtentCrs = rcp_dynamic_cast<CrsMatrixWrap>(Ptentative)->getCrsMatrix();

  ArrayRCP<size_t> iaPtent;
  ArrayRCP<LO> jaPtent;
  ArrayRCP<SC> valPtent;

  PtentCrs->allocateAllValues(nnzEstimate, iaPtent, jaPtent, valPtent);

  ArrayView<size_t> ia = iaPtent();
  ArrayView<LO> ja     = jaPtent();
  ArrayView<SC> val    = valPtent();

  ia[0] = 0;
  for (size_t i = 1; i <= numRows; i++)
    ia[i] = ia[i - 1] + NSDim;

  for (size_t j = 0; j < nnzEstimate; j++) {
    ja[j]  = INVALID;
    val[j] = zero;
  }

  if (doQRStep) {
    ////////////////////////////////
    // Standard aggregate-wise QR //
    ////////////////////////////////
    for (GO agg = 0; agg < numAggs; agg++) {
      LO aggSize = aggStart[agg + 1] - aggStart[agg];

      Xpetra::global_size_t offset = agg * NSDim;

      // Extract the piece of the nullspace corresponding to the aggregate, and
      // put it in the flat array, "localQR" (in column major format) for the
      // QR routine.
      Teuchos::SerialDenseMatrix<LO, SC> localQR(aggSize, NSDim);
      if (goodMap) {
        for (size_t j = 0; j < NSDim; j++)
          for (LO k = 0; k < aggSize; k++)
            localQR(k, j) = fineNS[j][aggToRowMapLO[aggStart[agg] + k]];
      } else {
        for (size_t j = 0; j < NSDim; j++)
          for (LO k = 0; k < aggSize; k++)
            localQR(k, j) = fineNS[j][rowMap->getLocalElement(aggToRowMapGO[aggStart[agg] + k])];
      }

      // Test for zero columns
      for (size_t j = 0; j < NSDim; j++) {
        bool bIsZeroNSColumn = true;

        for (LO k = 0; k < aggSize; k++)
          if (localQR(k, j) != zero)
            bIsZeroNSColumn = false;

        TEUCHOS_TEST_FOR_EXCEPTION(bIsZeroNSColumn == true, Exceptions::RuntimeError,
                                   "MueLu::TentativePFactory::MakeTentative: fine level NS part has a zero column in NS column " << j);
      }

      // Calculate QR decomposition (standard)
      // NOTE: Q is stored in localQR and R is stored in coarseNS
      if (aggSize >= Teuchos::as<LO>(NSDim)) {
        if (NSDim == 1) {
          // Only one nullspace vector, calculate Q and R by hand
          Magnitude norm = STS::magnitude(zero);
          for (size_t k = 0; k < Teuchos::as<size_t>(aggSize); k++)
            norm += STS::magnitude(localQR(k, 0) * localQR(k, 0));
          norm = Teuchos::ScalarTraits<Magnitude>::squareroot(norm);

          // R = norm
          coarseNS[0][offset] = norm;

          // Q = localQR(:,0)/norm
          for (LO i = 0; i < aggSize; i++)
            localQR(i, 0) /= norm;

        } else {
          Teuchos::SerialQRDenseSolver<LO, SC> qrSolver;
          qrSolver.setMatrix(Teuchos::rcp(&localQR, false));
          qrSolver.factor();

          // R = upper triangular part of localQR
          for (size_t j = 0; j < NSDim; j++)
            for (size_t k = 0; k <= j; k++)
              coarseNS[j][offset + k] = localQR(k, j);  // TODO is offset+k the correct local ID?!

          // Calculate Q, the tentative prolongator.
          // The Lapack GEQRF call only works for myAggsize >= NSDim
          qrSolver.formQ();
          Teuchos::RCP<Teuchos::SerialDenseMatrix<LO, SC> > qFactor = qrSolver.getQ();
          for (size_t j = 0; j < NSDim; j++)
            for (size_t i = 0; i < Teuchos::as<size_t>(aggSize); i++)
              localQR(i, j) = (*qFactor)(i, j);
        }

      } else {
        // Special handling for aggSize < NSDim (i.e. single node aggregates in structural mechanics)

        // The local QR decomposition is not possible in the "overconstrained"
        // case (i.e. number of columns in localQR > number of rows), which
        // corresponds to #DOFs in Aggregate < NSDim. For usual problems this
        // is only possible for single node aggregates in structural mechanics.
        // (Similar problems may arise in discontinuous Galerkin problems...)
        // We bypass the QR decomposition and use an identity block in the
        // tentative prolongator for the single node aggregate and transfer the
        // corresponding fine level null space information 1-to-1 to the coarse
        // level null space part.

        // NOTE: The resulting tentative prolongation operator has
        // (aggSize*DofsPerNode-NSDim) zero columns leading to a singular
        // coarse level operator A.  To deal with that one has the following
        // options:
        // - Use the "RepairMainDiagonal" flag in the RAPFactory (default:
        //   false) to add some identity block to the diagonal of the zero rows
        //   in the coarse level operator A, such that standard level smoothers
        //   can be used again.
        // - Use special (projection-based) level smoothers, which can deal
        //   with singular matrices (very application specific)
        // - Adapt the code below to avoid zero columns. However, we do not
        //   support a variable number of DOFs per node in MueLu/Xpetra which
        //   makes the implementation really hard.

        // R = extended (by adding identity rows) localQR
        for (size_t j = 0; j < NSDim; j++)
          for (size_t k = 0; k < NSDim; k++)
            if (k < as<size_t>(aggSize))
              coarseNS[j][offset + k] = localQR(k, j);
            else
              coarseNS[j][offset + k] = (k == j ? one : zero);

        // Q = I (rectangular)
        for (size_t i = 0; i < as<size_t>(aggSize); i++)
          for (size_t j = 0; j < NSDim; j++)
            localQR(i, j) = (j == i ? one : zero);
      }

      // Process each row in the local Q factor
      // FIXME: What happens if maps are blocked?
      for (LO j = 0; j < aggSize; j++) {
        LO localRow = (goodMap ? aggToRowMapLO[aggStart[agg] + j] : rowMap->getLocalElement(aggToRowMapGO[aggStart[agg] + j]));

        size_t rowStart = ia[localRow];
        for (size_t k = 0, lnnz = 0; k < NSDim; k++) {
          // Skip zeros (there may be plenty of them, i.e., NSDim > 1 or boundary conditions)
          if (localQR(j, k) != zero) {
            ja[rowStart + lnnz]  = offset + k;
            val[rowStart + lnnz] = localQR(j, k);
            lnnz++;
          }
        }
      }
    }

  } else {
    GetOStream(Runtime1) << "TentativePFactory : bypassing local QR phase" << std::endl;
    if (NSDim > 1)
      GetOStream(Warnings0) << "TentativePFactory : for nontrivial nullspace, this may degrade performance" << std::endl;
    /////////////////////////////
    //      "no-QR" option     //
    /////////////////////////////
    // Local Q factor is just the fine nullspace support over the current aggregate.
    // Local R factor is the identity.
    // TODO I have not implemented any special handling for aggregates that are too
    // TODO small to locally support the nullspace, as is done in the standard QR
    // TODO case above.
    if (goodMap) {
      for (GO agg = 0; agg < numAggs; agg++) {
        const LO aggSize             = aggStart[agg + 1] - aggStart[agg];
        Xpetra::global_size_t offset = agg * NSDim;

        // Process each row in the local Q factor
        // FIXME: What happens if maps are blocked?
        for (LO j = 0; j < aggSize; j++) {
          // TODO Here I do not check for a zero nullspace column on the aggregate.
          //      as is done in the standard QR case.

          const LO localRow = aggToRowMapLO[aggStart[agg] + j];

          const size_t rowStart = ia[localRow];

          for (size_t k = 0, lnnz = 0; k < NSDim; k++) {
            // Skip zeros (there may be plenty of them, i.e., NSDim > 1 or boundary conditions)
            SC qr_jk = fineNS[k][aggToRowMapLO[aggStart[agg] + j]];
            if (constantColSums) qr_jk = qr_jk / (Magnitude)aggSizes[agg];
            if (qr_jk != zero) {
              ja[rowStart + lnnz]  = offset + k;
              val[rowStart + lnnz] = qr_jk;
              lnnz++;
            }
          }
        }
        for (size_t j = 0; j < NSDim; j++)
          coarseNS[j][offset + j] = one;
      }  // for (GO agg = 0; agg < numAggs; agg++)

    } else {
      for (GO agg = 0; agg < numAggs; agg++) {
        const LO aggSize             = aggStart[agg + 1] - aggStart[agg];
        Xpetra::global_size_t offset = agg * NSDim;
        for (LO j = 0; j < aggSize; j++) {
          const LO localRow = rowMap->getLocalElement(aggToRowMapGO[aggStart[agg] + j]);

          const size_t rowStart = ia[localRow];

          for (size_t k = 0, lnnz = 0; k < NSDim; ++k) {
            // Skip zeros (there may be plenty of them, i.e., NSDim > 1 or boundary conditions)
            SC qr_jk = fineNS[k][rowMap->getLocalElement(aggToRowMapGO[aggStart[agg] + j])];
            if (constantColSums) qr_jk = qr_jk / (Magnitude)aggSizes[agg];
            if (qr_jk != zero) {
              ja[rowStart + lnnz]  = offset + k;
              val[rowStart + lnnz] = qr_jk;
              lnnz++;
            }
          }
        }
        for (size_t j = 0; j < NSDim; j++)
          coarseNS[j][offset + j] = one;
      }  // for (GO agg = 0; agg < numAggs; agg++)

    }  // if (goodmap) else ...

  }  // if doQRStep ... else

  // Compress storage (remove all INVALID, which happen when we skip zeros)
  // We do that in-place
  size_t ia_tmp = 0, nnz = 0;
  for (size_t i = 0; i < numRows; i++) {
    for (size_t j = ia_tmp; j < ia[i + 1]; j++)
      if (ja[j] != INVALID) {
        ja[nnz]  = ja[j];
        val[nnz] = val[j];
        nnz++;
      }
    ia_tmp    = ia[i + 1];
    ia[i + 1] = nnz;
  }
  if (rowMap->lib() == Xpetra::UseTpetra) {
    // - Cannot resize for Epetra, as it checks for same pointers
    // - Need to resize for Tpetra, as it check ().size() == ia[numRows]
    // NOTE: these invalidate ja and val views
    jaPtent.resize(nnz);
    valPtent.resize(nnz);
  }

  GetOStream(Runtime1) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;

  PtentCrs->setAllValues(iaPtent, jaPtent, valPtent);

  // Managing labels & constants for ESFC
  RCP<ParameterList> FCparams;
  if (pL.isSublist("matrixmatrix: kernel params"))
    FCparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
  else
    FCparams = rcp(new ParameterList);
  // By default, we don't need global constants for TentativeP
  FCparams->set("compute global constants", FCparams->get("compute global constants", false));
  std::string levelIDs = toString(levelID);
  FCparams->set("Timer Label", std::string("MueLu::TentativeP-") + levelIDs);
  RCP<const Export> dummy_e;
  RCP<const Import> dummy_i;

  PtentCrs->expertStaticFillComplete(coarseMap, A->getDomainMap(), dummy_i, dummy_e, FCparams);
}

}  // namespace MueLu

// TODO ReUse: If only P or Nullspace is missing, TentativePFactory can be smart and skip part of the computation.

#define MUELU_TENTATIVEPFACTORY_SHORT
#endif  // MUELU_TENTATIVEPFACTORY_DEF_HPP
