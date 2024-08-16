// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP
#define MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP

#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_TentativePFactory_kokkos_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_LocalQR.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TentativePFactory_kokkos() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TentativePFactory_kokkos() = default;


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("tentative: build coarse coordinates");
  // SET_VALID_ENTRY("tentative: constant column sums");
#undef SET_VALID_ENTRY
  validParamList->set<std::string>("Nullspace name", "Nullspace", "Name for the input nullspace");

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase>>("Aggregates", Teuchos::null, "Generating factory of the aggregates");
  validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<RCP<const FactoryBase>>("Scaled Nullspace", Teuchos::null, "Generating factory of the scaled nullspace");
  validParamList->set<RCP<const FactoryBase>>("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase>>("CoarseMap", Teuchos::null, "Generating factory of the coarse map");
  validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Generating factory of the coordinates");
  validParamList->set<RCP<const FactoryBase>>("Node Comm", Teuchos::null, "Generating factory of the node level communicator");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
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
void TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinate_type;
  typedef Xpetra::MultiVector<coordinate_type, LO, GO, NO> RealValuedMultiVector;
  typedef Xpetra::MultiVectorFactory<coordinate_type, LO, GO, NO> RealValuedMultiVectorFactory;

  const ParameterList& pL = GetParameterList();
  std::string nspName     = "Nullspace";
  if (pL.isParameter("Nullspace name")) nspName = pL.get<std::string>("Nullspace name");

  RCP<Matrix> Ptentative;
  auto A          = Get<RCP<Matrix>>(fineLevel, "A");
  auto aggregates = Get<RCP<Aggregates>>(fineLevel, "Aggregates");
  // No coarse DoFs so we need to bail by setting Ptentattive to null and returning
  //  This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
  if (aggregates->GetNumGlobalAggregatesComputeIfNeeded() == 0) {
    Ptentative = Teuchos::null;
    Set(coarseLevel, "P", Ptentative);
    return;
  }

  auto amalgInfo     = Get<RCP<AmalgamationInfo>>(fineLevel, "UnAmalgamationInfo");
  auto fineNullspace = Get<RCP<MultiVector>>(fineLevel, nspName);
  auto coarseMap     = Get<RCP<const Map>>(fineLevel, "CoarseMap");
  RCP<RealValuedMultiVector> fineCoords;
  if (bTransferCoordinates_) {
    fineCoords = Get<RCP<RealValuedMultiVector>>(fineLevel, "Coordinates");
  }

  // FIXME: We should remove the NodeComm on levels past the threshold
  if (fineLevel.IsAvailable("Node Comm")) {
    RCP<const Teuchos::Comm<int>> nodeComm = Get<RCP<const Teuchos::Comm<int>>>(fineLevel, "Node Comm");
    Set<RCP<const Teuchos::Comm<int>>>(coarseLevel, "Node Comm", nodeComm);
  }

  // NOTE:  We check DomainMap here rather than RowMap because those are different for BlockCrs matrices
  TEUCHOS_TEST_FOR_EXCEPTION(A->getDomainMap()->getLocalNumElements() != fineNullspace->getMap()->getLocalNumElements(),
                             Exceptions::RuntimeError, "MueLu::TentativePFactory::MakeTentative: Size mismatch between A and Nullspace");

  RCP<MultiVector> coarseNullspace;
  RCP<RealValuedMultiVector> coarseCoords;

  if (bTransferCoordinates_) {
    RCP<const Map> coarseCoordMap;

    LO blkSize = 1;
    if (rcp_dynamic_cast<const StridedMap>(coarseMap) != Teuchos::null)
      blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap)->getFixedBlockSize();

    if (blkSize == 1) {
      // Scalar system
      // No amalgamation required, we can use the coarseMap
      coarseCoordMap = coarseMap;
    } else {
      // Vector system
      AmalgamationFactory<SC, LO, GO, NO>::AmalgamateMap(rcp_dynamic_cast<const StridedMap>(coarseMap), coarseCoordMap);
    }

    coarseCoords = RealValuedMultiVectorFactory::Build(coarseCoordMap, fineCoords->getNumVectors());

    // Create overlapped fine coordinates to reduce global communication
    auto uniqueMap                           = fineCoords->getMap();
    RCP<RealValuedMultiVector> ghostedCoords = fineCoords;
    if (aggregates->AggregatesCrossProcessors()) {
      auto nonUniqueMap = aggregates->GetMap();
      auto importer     = ImportFactory::Build(uniqueMap, nonUniqueMap);

      ghostedCoords = RealValuedMultiVectorFactory::Build(nonUniqueMap, fineCoords->getNumVectors());
      ghostedCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);
    }

    // The good new is that his graph has already been constructed for the
    // TentativePFactory and was cached in Aggregates. So this is a no-op.
    auto aggGraph = aggregates->GetGraph();
    auto numAggs  = aggGraph.numRows();

    auto fineCoordsView   = fineCoords->getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto coarseCoordsView = coarseCoords->getLocalViewDevice(Tpetra::Access::OverwriteAll);

    // Fill in coarse coordinates
    {
      SubFactoryMonitor m2(*this, "AverageCoords", coarseLevel);

      const auto dim = fineCoords->getNumVectors();

      typename AppendTrait<decltype(fineCoordsView), Kokkos::RandomAccess>::type fineCoordsRandomView = fineCoordsView;
      for (size_t j = 0; j < dim; j++) {
        Kokkos::parallel_for(
            "MueLu::TentativeP::BuildCoords", Kokkos::RangePolicy<LocalOrdinal, execution_space>(0, numAggs),
            KOKKOS_LAMBDA(const LO i) {
              // A row in this graph represents all node ids in the aggregate
              // Therefore, averaging is very easy

              auto aggregate = aggGraph.rowConst(i);

              coordinate_type sum = 0.0;  // do not use Scalar here (Stokhos)
              for (size_t colID = 0; colID < static_cast<size_t>(aggregate.length); colID++)
                sum += fineCoordsRandomView(aggregate(colID), j);

              coarseCoordsView(i, j) = sum / aggregate.length;
            });
      }
    }
  }

  if (!aggregates->AggregatesCrossProcessors()) {
    if (Xpetra::Helpers<SC, LO, GO, NO>::isTpetraBlockCrs(A)) {
      BuildPuncoupledBlockCrs(coarseLevel, A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace,
                              coarseLevel.GetLevelID());
    } else {
      BuildPuncoupled(coarseLevel, A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace, coarseLevel.GetLevelID());
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
void TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildPuncoupled(Level& coarseLevel, RCP<Matrix> A, RCP<Aggregates> aggregates,
                    RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                    RCP<const Map> coarseMap, RCP<Matrix>& Ptentative,
                    RCP<MultiVector>& coarseNullspace, const int levelID) const {
  auto rowMap = A->getRowMap();
  auto colMap = A->getColMap();

  const size_t numRows = rowMap->getLocalNumElements();
  const size_t NSDim   = fineNullspace->getNumVectors();

#if KOKKOS_VERSION >= 40799
  typedef KokkosKernels::ArithTraits<SC> ATS;
#else
  typedef Kokkos::ArithTraits<SC> ATS;
#endif
  using impl_SC      = typename ATS::val_type;
#if KOKKOS_VERSION >= 40799
  using impl_ATS = KokkosKernels::ArithTraits<impl_SC>;
#else
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;
#endif
  const impl_SC zero = impl_ATS::zero();
  const impl_SC one  = impl_ATS::one();

  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  typename Aggregates::local_graph_type aggGraph;
  {
    SubFactoryMonitor m2(*this, "Get Aggregates graph", coarseLevel);
    aggGraph = aggregates->GetGraph();
  }
  auto aggRows = aggGraph.row_map;
  auto aggCols = aggGraph.entries;

  // Aggregates map is based on the amalgamated column map
  // We can skip global-to-local conversion if LIDs in row map are
  // same as LIDs in column map
  bool goodMap;
  {
    SubFactoryMonitor m2(*this, "Check good map", coarseLevel);
    goodMap = Utilities::MapsAreNested(*rowMap, *colMap);
  }
  // FIXME_KOKKOS: need to proofread later code for bad maps
  TEUCHOS_TEST_FOR_EXCEPTION(!goodMap, Exceptions::RuntimeError,
                             "MueLu: TentativePFactory_kokkos: for now works only with good maps "
                             "(i.e. \"matching\" row and column maps)");

  // STEP 1: do unamalgamation
  // The non-kokkos version uses member functions from the AmalgamationInfo
  // container class to unamalgamate the data. In contrast, the kokkos
  // version of TentativePFactory does the unamalgamation here and only uses
  // the data of the AmalgamationInfo container class

  // Extract information for unamalgamation
  LO fullBlockSize;
  LO blockID;
  LO stridingOffset;
  LO stridedBlockSize;
  GO indexBase;
  amalgInfo->GetStridingInformation(fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase);
  GO globalOffset = amalgInfo->GlobalOffset();

  // Extract aggregation info (already in Kokkos host views)
  auto procWinner            = aggregates->GetProcWinner()->getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto vertex2AggId          = aggregates->GetVertex2AggId()->getLocalViewDevice(Tpetra::Access::ReadOnly);
  const size_t numAggregates = aggregates->GetNumAggregates();

  int myPID = aggregates->GetMap()->getComm()->getRank();

  // Create Kokkos::View (on the device) to store the aggreate dof sizes
  // Later used to get aggregate dof offsets
  // NOTE: This zeros itself on construction
  typedef typename Aggregates::aggregates_sizes_type::non_const_type AggSizeType;
  AggSizeType aggDofSizes;

  if (stridedBlockSize == 1) {
    SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

    // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
    aggDofSizes = AggSizeType("agg_dof_sizes", numAggregates + 1);

    auto sizesConst = aggregates->ComputeAggregateSizes();
    Kokkos::deep_copy(Kokkos::subview(aggDofSizes, Kokkos::make_pair(static_cast<size_t>(1), numAggregates + 1)), sizesConst);

  } else {
    SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

    // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
    aggDofSizes = AggSizeType("agg_dof_sizes", numAggregates + 1);

    auto nodeMap = aggregates->GetMap()->getLocalMap();
    auto dofMap  = colMap->getLocalMap();

    Kokkos::parallel_for(
        "MueLu:TentativePF:Build:compute_agg_sizes", range_type(0, numAggregates),
        KOKKOS_LAMBDA(const LO agg) {
          auto aggRowView = aggGraph.rowConst(agg);

          size_t size = 0;
          for (LO colID = 0; colID < aggRowView.length; colID++) {
            GO nodeGID = nodeMap.getGlobalElement(aggRowView(colID));

            for (LO k = 0; k < stridedBlockSize; k++) {
              GO dofGID = (nodeGID - indexBase) * fullBlockSize + k + indexBase + globalOffset + stridingOffset;

              if (dofMap.getLocalElement(dofGID) != INVALID)
                size++;
            }
          }
          aggDofSizes(agg + 1) = size;
        });
  }

  // Find maximum dof size for aggregates
  // Later used to reserve enough scratch space for local QR decompositions
  LO maxAggSize = 0;
  LocalQR::ReduceMaxFunctor<LO, decltype(aggDofSizes)> reduceMax(aggDofSizes);
  Kokkos::parallel_reduce("MueLu:TentativePF:Build:max_agg_size", range_type(0, aggDofSizes.extent(0)), reduceMax, maxAggSize);

  // parallel_scan (exclusive)
  // The aggDofSizes View then contains the aggregate dof offsets
  Kokkos::parallel_scan(
      "MueLu:TentativePF:Build:aggregate_sizes:stage1_scan", range_type(0, numAggregates + 1),
      KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
        update += aggDofSizes(i);
        if (final_pass)
          aggDofSizes(i) = update;
      });

  // Create Kokkos::View on the device to store mapping
  // between (local) aggregate id and row map ids (LIDs)
  Kokkos::View<LO*, DeviceType> agg2RowMapLO(Kokkos::ViewAllocateWithoutInitializing("agg2row_map_LO"), numRows);
  {
    SubFactoryMonitor m2(*this, "Create Agg2RowMap", coarseLevel);

    AggSizeType aggOffsets(Kokkos::ViewAllocateWithoutInitializing("aggOffsets"), numAggregates);
    Kokkos::deep_copy(aggOffsets, Kokkos::subview(aggDofSizes, Kokkos::make_pair(static_cast<size_t>(0), numAggregates)));

    Kokkos::parallel_for(
        "MueLu:TentativePF:Build:createAgg2RowMap", range_type(0, vertex2AggId.extent(0)),
        KOKKOS_LAMBDA(const LO lnode) {
          if (procWinner(lnode, 0) == myPID) {
            // No need for atomics, it's one-to-one
            auto aggID = vertex2AggId(lnode, 0);

            auto offset = Kokkos::atomic_fetch_add(&aggOffsets(aggID), stridedBlockSize);
            // FIXME: I think this may be wrong
            // We unconditionally add the whole block here. When we calculated
            // aggDofSizes, we did the isLocalElement check. Something's fishy.
            for (LO k = 0; k < stridedBlockSize; k++)
              agg2RowMapLO(offset + k) = lnode * stridedBlockSize + k;
          }
        });
  }

  // STEP 2: prepare local QR decomposition
  // Reserve memory for tentative prolongation operator
  coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

  // Pull out the nullspace vectors so that we can have random access (on the device)
  auto fineNS   = fineNullspace->getLocalViewDevice(Tpetra::Access::ReadWrite);
  auto coarseNS = coarseNullspace->getLocalViewDevice(Tpetra::Access::OverwriteAll);

  size_t nnz = 0;  // actual number of nnz

  typedef typename Xpetra::Matrix<SC, LO, GO, NO>::local_matrix_device_type local_matrix_type;
  typedef typename local_matrix_type::row_map_type::non_const_type rows_type;
  typedef typename local_matrix_type::index_type::non_const_type cols_type;
  typedef typename local_matrix_type::values_type::non_const_type vals_type;

  // Device View for status (error messages...)
  typedef Kokkos::View<int[10], DeviceType> status_type;
  status_type status("status");

  typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
  typename AppendTrait<status_type, Kokkos::Atomic>::type statusAtomic            = status;

  const ParameterList& pL = GetParameterList();
  const bool& doQRStep    = pL.get<bool>("tentative: calculate qr");
  if (!doQRStep) {
    GetOStream(Runtime1) << "TentativePFactory : bypassing local QR phase" << std::endl;
    if (NSDim > 1)
      GetOStream(Warnings0) << "TentativePFactor : for nontrivial nullspace, this may degrade performance" << std::endl;
  }

  size_t nnzEstimate = numRows * NSDim;
  rows_type rowsAux(Kokkos::ViewAllocateWithoutInitializing("Ptent_aux_rows"), numRows + 1);
  cols_type colsAux(Kokkos::ViewAllocateWithoutInitializing("Ptent_aux_cols"), nnzEstimate);
  vals_type valsAux("Ptent_aux_vals", nnzEstimate);
  rows_type rows("Ptent_rows", numRows + 1);
  {
    // Stage 0: fill in views.
    SubFactoryMonitor m2(*this, "Stage 0 (InitViews)", coarseLevel);

    // The main thing to notice is initialization of vals with INVALID. These
    // values will later be used to compress the arrays
    Kokkos::parallel_for(
        "MueLu:TentativePF:BuildPuncoupled:for1", range_type(0, numRows + 1),
        KOKKOS_LAMBDA(const LO row) {
          rowsAux(row) = row * NSDim;
        });
    Kokkos::parallel_for(
        "MueLu:TentativePF:BuildUncoupled:for2", range_type(0, nnzEstimate),
        KOKKOS_LAMBDA(const LO j) {
          colsAux(j) = INVALID;
        });
  }

  if (NSDim == 1) {
    // 1D is special, as it is the easiest. We don't even need to the QR,
    // just normalize an array. Plus, no worries abot small aggregates.  In
    // addition, we do not worry about compression. It is unlikely that
    // nullspace will have zeros. If it does, a prolongator row would be
    // zero and we'll get singularity anyway.
    SubFactoryMonitor m2(*this, "Stage 1 (LocalQR)", coarseLevel);

    // Set up team policy with numAggregates teams and one thread per team.
    // Each team handles a slice of the data associated with one aggregate
    // and performs a local QR decomposition (in this case real QR is
    // unnecessary).
    const Kokkos::TeamPolicy<execution_space> policy(numAggregates, 1);

    if (doQRStep) {
      Kokkos::parallel_for(
          "MueLu:TentativePF:BuildUncoupled:main_loop", policy,
          KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread) {
            auto agg = thread.league_rank();

            // size of the aggregate (number of DOFs in aggregate)
            LO aggSize = aggRows(agg + 1) - aggRows(agg);

            // Extract the piece of the nullspace corresponding to the aggregate, and
            // put it in the flat array, "localQR" (in column major format) for the
            // QR routine. Trivial in 1D.
            auto norm = impl_ATS::magnitude(zero);

            // Calculate QR by hand
            // FIXME: shouldn't there be stridedblock here?
            // FIXME_KOKKOS: shouldn't there be stridedblock here?
            for (decltype(aggSize) k = 0; k < aggSize; k++) {
              auto dnorm = impl_ATS::magnitude(fineNSRandom(agg2RowMapLO(aggRows(agg) + k), 0));
              norm += dnorm * dnorm;
            }
            norm = sqrt(norm);

            if (norm == zero) {
              // zero column; terminate the execution
              statusAtomic(1) = true;
              return;
            }

            // R = norm
            coarseNS(agg, 0) = norm;

            // Q = localQR(:,0)/norm
            for (decltype(aggSize) k = 0; k < aggSize; k++) {
              LO localRow      = agg2RowMapLO(aggRows(agg) + k);
              impl_SC localVal = fineNSRandom(agg2RowMapLO(aggRows(agg) + k), 0) / norm;

              rows(localRow + 1) = 1;
              colsAux(localRow)  = agg;
              valsAux(localRow)  = localVal;
            }
          });

      typename status_type::host_mirror_type statusHost = Kokkos::create_mirror_view(status);
      Kokkos::deep_copy(statusHost, status);
      for (decltype(statusHost.size()) i = 0; i < statusHost.size(); i++)
        if (statusHost(i)) {
          std::ostringstream oss;
          oss << "MueLu::TentativePFactory::MakeTentative: ";
          switch (i) {
            case 0: oss << "!goodMap is not implemented"; break;
            case 1: oss << "fine level NS part has a zero column"; break;
          }
          throw Exceptions::RuntimeError(oss.str());
        }

    } else {
      Kokkos::parallel_for(
          "MueLu:TentativePF:BuildUncoupled:main_loop_noqr", policy,
          KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread) {
            auto agg = thread.league_rank();

            // size of the aggregate (number of DOFs in aggregate)
            LO aggSize = aggRows(agg + 1) - aggRows(agg);

            // R = norm
            coarseNS(agg, 0) = one;

            // Q = localQR(:,0)/norm
            for (decltype(aggSize) k = 0; k < aggSize; k++) {
              LO localRow      = agg2RowMapLO(aggRows(agg) + k);
              impl_SC localVal = fineNSRandom(agg2RowMapLO(aggRows(agg) + k), 0);

              rows(localRow + 1) = 1;
              colsAux(localRow)  = agg;
              valsAux(localRow)  = localVal;
            }
          });
    }

    Kokkos::parallel_reduce(
        "MueLu:TentativeP:CountNNZ", range_type(0, numRows + 1),
        KOKKOS_LAMBDA(const LO i, size_t& nnz_count) {
          nnz_count += rows(i);
        },
        nnz);

  } else {  // NSdim > 1
    // FIXME_KOKKOS: This code branch is completely unoptimized.
    // Work to do:
    //   - Optimize QR decomposition
    //   - Remove INVALID usage similarly to CoalesceDropFactory_kokkos by
    //     packing new values in the beginning of each row
    // We do use auxilary view in this case, so keep a second rows view for
    // counting nonzeros in rows

    {
      SubFactoryMonitor m2 = SubFactoryMonitor(*this, doQRStep ? "Stage 1 (LocalQR)" : "Stage 1 (Fill coarse nullspace and tentative P)", coarseLevel);
      // Set up team policy with numAggregates teams and one thread per team.
      // Each team handles a slice of the data associated with one aggregate
      // and performs a local QR decomposition
      const Kokkos::TeamPolicy<execution_space> policy(numAggregates, 1);  // numAggregates teams a 1 thread
      LocalQR::LocalQRDecompFunctor<LocalOrdinal, GlobalOrdinal, Scalar, DeviceType, decltype(fineNSRandom),
                                    decltype(aggDofSizes /*aggregate sizes in dofs*/), decltype(maxAggSize), decltype(agg2RowMapLO),
                                    decltype(statusAtomic), decltype(rows), decltype(rowsAux), decltype(colsAux),
                                    decltype(valsAux)>
          localQRFunctor(fineNSRandom, coarseNS, aggDofSizes, maxAggSize, agg2RowMapLO, statusAtomic,
                         rows, rowsAux, colsAux, valsAux, doQRStep);
      Kokkos::parallel_reduce("MueLu:TentativePF:BuildUncoupled:main_qr_loop", policy, localQRFunctor, nnz);
    }

    typename status_type::host_mirror_type statusHost = Kokkos::create_mirror_view(status);
    Kokkos::deep_copy(statusHost, status);
    for (decltype(statusHost.size()) i = 0; i < statusHost.size(); i++)
      if (statusHost(i)) {
        std::ostringstream oss;
        oss << "MueLu::TentativePFactory::MakeTentative: ";
        switch (i) {
          case 0: oss << "!goodMap is not implemented"; break;
          case 1: oss << "fine level NS part has a zero column"; break;
        }
        throw Exceptions::RuntimeError(oss.str());
      }
  }

  // Compress the cols and vals by ignoring INVALID column entries that correspond
  // to 0 in QR.

  // The real cols and vals are constructed using calculated (not estimated) nnz
  cols_type cols;
  vals_type vals;

  if (nnz != nnzEstimate) {
    {
      // Stage 2: compress the arrays
      SubFactoryMonitor m2(*this, "Stage 2 (CompressRows)", coarseLevel);

      Kokkos::parallel_scan(
          "MueLu:TentativePF:Build:compress_rows", range_type(0, numRows + 1),
          KOKKOS_LAMBDA(const LO i, LO& upd, const bool& final) {
            upd += rows(i);
            if (final)
              rows(i) = upd;
          });
    }

    {
      SubFactoryMonitor m2(*this, "Stage 2 (CompressCols)", coarseLevel);

      cols = cols_type("Ptent_cols", nnz);
      vals = vals_type("Ptent_vals", nnz);

      // FIXME_KOKKOS: this can be spedup by moving correct cols and vals values
      // to the beginning of rows. See CoalesceDropFactory_kokkos for
      // example.
      Kokkos::parallel_for(
          "MueLu:TentativePF:Build:compress_cols_vals", range_type(0, numRows),
          KOKKOS_LAMBDA(const LO i) {
            LO rowStart = rows(i);

            size_t lnnz = 0;
            for (auto j = rowsAux(i); j < rowsAux(i + 1); j++)
              if (colsAux(j) != INVALID) {
                cols(rowStart + lnnz) = colsAux(j);
                vals(rowStart + lnnz) = valsAux(j);
                lnnz++;
              }
          });
    }

  } else {
    rows = rowsAux;
    cols = colsAux;
    vals = valsAux;
  }

  GetOStream(Runtime1) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;

  {
    // Stage 3: construct Xpetra::Matrix
    SubFactoryMonitor m2(*this, "Stage 3 (LocalMatrix+FillComplete)", coarseLevel);

    local_matrix_type lclMatrix = local_matrix_type("A", numRows, coarseMap->getLocalNumElements(), nnz, vals, rows, cols);

    // Managing labels & constants for ESFC
    RCP<ParameterList> FCparams;
    if (pL.isSublist("matrixmatrix: kernel params"))
      FCparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
    else
      FCparams = rcp(new ParameterList);

    // By default, we don't need global constants for TentativeP
    FCparams->set("compute global constants", FCparams->get("compute global constants", false));
    FCparams->set("Timer Label", std::string("MueLu::TentativeP-") + toString(levelID));

    auto PtentCrs = CrsMatrixFactory::Build(lclMatrix, rowMap, coarseMap, coarseMap, A->getDomainMap());
    Ptentative    = rcp(new CrsMatrixWrap(PtentCrs));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildPuncoupledBlockCrs(Level& coarseLevel, RCP<Matrix> A, RCP<Aggregates> aggregates,
                            RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                            RCP<const Map> coarsePointMap, RCP<Matrix>& Ptentative,
                            RCP<MultiVector>& coarseNullspace, const int levelID) const {
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

  // typedef Teuchos::ScalarTraits<SC> STS;
  // typedef typename STS::magnitudeType Magnitude;
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

#if KOKKOS_VERSION >= 40799
  typedef KokkosKernels::ArithTraits<SC> ATS;
#else
  typedef Kokkos::ArithTraits<SC> ATS;
#endif
  using impl_SC = typename ATS::val_type;
#if KOKKOS_VERSION >= 40799
  using impl_ATS = KokkosKernels::ArithTraits<impl_SC>;
#else
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;
#endif
  const impl_SC one = impl_ATS::one();

  //    const GO     numAggs   = aggregates->GetNumAggregates();
  const size_t NSDim = fineNullspace->getNumVectors();
  auto aggSizes      = aggregates->ComputeAggregateSizes();

  typename Aggregates::local_graph_type aggGraph;
  {
    SubFactoryMonitor m2(*this, "Get Aggregates graph", coarseLevel);
    aggGraph = aggregates->GetGraph();
  }
  auto aggRows = aggGraph.row_map;
  auto aggCols = aggGraph.entries;

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
  const ParameterList& pL = GetParameterList();
  //    const bool &doQRStep = pL.get<bool>("tentative: calculate qr");

  // The aggregates use the amalgamated column map, which in this case is what we want

  // Aggregates map is based on the amalgamated column map
  // We can skip global-to-local conversion if LIDs in row map are
  // same as LIDs in column map
  bool goodMap = Utilities::MapsAreNested(*rowMap, *colMap);
  TEUCHOS_TEST_FOR_EXCEPTION(!goodMap, Exceptions::RuntimeError,
                             "MueLu: TentativePFactory_kokkos: for now works only with good maps "
                             "(i.e. \"matching\" row and column maps)");

  // STEP 1: do unamalgamation
  // The non-kokkos version uses member functions from the AmalgamationInfo
  // container class to unamalgamate the data. In contrast, the kokkos
  // version of TentativePFactory does the unamalgamation here and only uses
  // the data of the AmalgamationInfo container class

  // Extract information for unamalgamation
  LO fullBlockSize;
  LO blockID;
  LO stridingOffset;
  LO stridedBlockSize;
  GO indexBase;
  amalgInfo->GetStridingInformation(fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase);
  // GO globalOffset = amalgInfo->GlobalOffset();

  // Extract aggregation info (already in Kokkos host views)
  auto procWinner            = aggregates->GetProcWinner()->getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto vertex2AggId          = aggregates->GetVertex2AggId()->getLocalViewDevice(Tpetra::Access::ReadOnly);
  const size_t numAggregates = aggregates->GetNumAggregates();

  int myPID = aggregates->GetMap()->getComm()->getRank();

  // Create Kokkos::View (on the device) to store the aggreate dof sizes
  // Later used to get aggregate dof offsets
  // NOTE: This zeros itself on construction
  typedef typename Aggregates::aggregates_sizes_type::non_const_type AggSizeType;
  AggSizeType aggDofSizes;  // This turns into "starts" after the parallel_scan

  {
    SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

    // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
    aggDofSizes = AggSizeType("agg_dof_sizes", numAggregates + 1);

    Kokkos::deep_copy(Kokkos::subview(aggDofSizes, Kokkos::make_pair(static_cast<size_t>(1), numAggregates + 1)), aggSizes);
  }

  // Find maximum dof size for aggregates
  // Later used to reserve enough scratch space for local QR decompositions
  LO maxAggSize = 0;
  LocalQR::ReduceMaxFunctor<LO, decltype(aggDofSizes)> reduceMax(aggDofSizes);
  Kokkos::parallel_reduce("MueLu:TentativePF:Build:max_agg_size", range_type(0, aggDofSizes.extent(0)), reduceMax, maxAggSize);

  // parallel_scan (exclusive)
  // The aggDofSizes View then contains the aggregate dof offsets
  Kokkos::parallel_scan(
      "MueLu:TentativePF:Build:aggregate_sizes:stage1_scan", range_type(0, numAggregates + 1),
      KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
        update += aggDofSizes(i);
        if (final_pass)
          aggDofSizes(i) = update;
      });

  // Create Kokkos::View on the device to store mapping
  // between (local) aggregate id and row map ids (LIDs)
  Kokkos::View<LO*, DeviceType> aggToRowMapLO(Kokkos::ViewAllocateWithoutInitializing("aggtorow_map_LO"), numFineBlockRows);
  {
    SubFactoryMonitor m2(*this, "Create AggToRowMap", coarseLevel);

    AggSizeType aggOffsets(Kokkos::ViewAllocateWithoutInitializing("aggOffsets"), numAggregates);
    Kokkos::deep_copy(aggOffsets, Kokkos::subview(aggDofSizes, Kokkos::make_pair(static_cast<size_t>(0), numAggregates)));

    Kokkos::parallel_for(
        "MueLu:TentativePF:Build:createAgg2RowMap", range_type(0, vertex2AggId.extent(0)),
        KOKKOS_LAMBDA(const LO lnode) {
          if (procWinner(lnode, 0) == myPID) {
            // No need for atomics, it's one-to-one
            auto aggID = vertex2AggId(lnode, 0);

            auto offset = Kokkos::atomic_fetch_add(&aggOffsets(aggID), stridedBlockSize);
            // FIXME: I think this may be wrong
            // We unconditionally add the whole block here. When we calculated
            // aggDofSizes, we did the isLocalElement check. Something's fishy.
            for (LO k = 0; k < stridedBlockSize; k++)
              aggToRowMapLO(offset + k) = lnode * stridedBlockSize + k;
          }
        });
  }

  // STEP 2: prepare local QR decomposition
  // Reserve memory for tentative prolongation operator
  coarseNullspace = MultiVectorFactory::Build(coarsePointMap, NSDim);

  // Pull out the nullspace vectors so that we can have random access (on the device)
  auto fineNS   = fineNullspace->getLocalViewDevice(Tpetra::Access::ReadWrite);
  auto coarseNS = coarseNullspace->getLocalViewDevice(Tpetra::Access::OverwriteAll);

  typedef typename Xpetra::Matrix<SC, LO, GO, NO>::local_matrix_device_type local_matrix_type;
  typedef typename local_matrix_type::row_map_type::non_const_type rows_type;
  typedef typename local_matrix_type::index_type::non_const_type cols_type;
  // typedef typename local_matrix_type::values_type::non_const_type    vals_type;

  // Device View for status (error messages...)
  typedef Kokkos::View<int[10], DeviceType> status_type;
  status_type status("status");

  typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
  typename AppendTrait<status_type, Kokkos::Atomic>::type statusAtomic            = status;

  // We're going to bypass QR in the BlockCrs version of the code regardless of what the user asks for
  GetOStream(Runtime1) << "TentativePFactory : bypassing local QR phase" << std::endl;

  // BlockCrs requires that we build the (block) graph first, so let's do that...

  // NOTE: Because we're assuming that the NSDim == BlockSize, we only have one
  // block non-zero per row in the matrix;
  rows_type ia(Kokkos::ViewAllocateWithoutInitializing("BlockGraph_rowptr"), numFineBlockRows + 1);
  cols_type ja(Kokkos::ViewAllocateWithoutInitializing("BlockGraph_colind"), numFineBlockRows);

  Kokkos::parallel_for(
      "MueLu:TentativePF:BlockCrs:graph_init", range_type(0, numFineBlockRows),
      KOKKOS_LAMBDA(const LO j) {
        ia[j] = j;
        ja[j] = INVALID;

        if (j == (LO)numFineBlockRows - 1)
          ia[numFineBlockRows] = numFineBlockRows;
      });

  // Fill Graph
  const Kokkos::TeamPolicy<execution_space> policy(numAggregates, 1);
  Kokkos::parallel_for(
      "MueLu:TentativePF:BlockCrs:fillGraph", policy,
      KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread) {
        auto agg                     = thread.league_rank();
        Xpetra::global_size_t offset = agg;

        // size of the aggregate (number of DOFs in aggregate)
        LO aggSize = aggRows(agg + 1) - aggRows(agg);

        for (LO j = 0; j < aggSize; j++) {
          // FIXME: Allow for bad maps
          const LO localRow     = aggToRowMapLO[aggDofSizes[agg] + j];
          const size_t rowStart = ia[localRow];
          ja[rowStart]          = offset;
        }
      });

  // Compress storage (remove all INVALID, which happen when we skip zeros)
  // We do that in-place
  {
    // Stage 2: compress the arrays
    SubFactoryMonitor m2(*this, "Stage 2 (CompressData)", coarseLevel);
    // Fill i_temp with the correct row starts
    rows_type i_temp(Kokkos::ViewAllocateWithoutInitializing("BlockGraph_rowptr"), numFineBlockRows + 1);
    LO nnz = 0;
    Kokkos::parallel_scan(
        "MueLu:TentativePF:BlockCrs:compress_rows", range_type(0, numFineBlockRows),
        KOKKOS_LAMBDA(const LO i, LO& upd, const bool& final) {
          if (final)
            i_temp[i] = upd;
          for (auto j = ia[i]; j < ia[i + 1]; j++)
            if (ja[j] != INVALID)
              upd++;
          if (final && i == (LO)numFineBlockRows - 1)
            i_temp[numFineBlockRows] = upd;
        },
        nnz);

    cols_type j_temp(Kokkos::ViewAllocateWithoutInitializing("BlockGraph_colind"), nnz);

    Kokkos::parallel_for(
        "MueLu:TentativePF:BlockCrs:compress_cols", range_type(0, numFineBlockRows),
        KOKKOS_LAMBDA(const LO i) {
          size_t rowStart = i_temp[i];
          size_t lnnz     = 0;
          for (auto j = ia[i]; j < ia[i + 1]; j++)
            if (ja[j] != INVALID) {
              j_temp[rowStart + lnnz] = ja[j];
              lnnz++;
            }
        });

    ia = i_temp;
    ja = j_temp;
  }

  RCP<CrsGraph> BlockGraph = CrsGraphFactory::Build(rowMap, coarseBlockMap, ia, ja);

  // Managing labels & constants for ESFC
  {
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
    BlockGraph->expertStaticFillComplete(coarseBlockMap, rowMap, dummy_i, dummy_e, FCparams);
  }

  // We can't leave the ia/ja pointers floating around, because of host/device view counting, so
  // we clear them here
  ia = rows_type();
  ja = cols_type();

  // Now let's make a BlockCrs Matrix
  // NOTE: Assumes block size== NSDim
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> P_xpetra            = Xpetra::CrsMatrixFactory<SC, LO, GO, NO>::BuildBlock(BlockGraph, coarsePointMap, rangeMap, NSDim);
  RCP<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>> P_tpetra = rcp_dynamic_cast<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>>(P_xpetra);
  if (P_tpetra.is_null()) throw std::runtime_error("BuildPUncoupled: Matrix factory did not return a Tpetra::BlockCrsMatrix");
  RCP<CrsMatrixWrap> P_wrap = rcp(new CrsMatrixWrap(P_xpetra));

  auto values     = P_tpetra->getTpetra_BlockCrsMatrix()->getValuesDeviceNonConst();
  const LO stride = NSDim * NSDim;

  Kokkos::parallel_for(
      "MueLu:TentativePF:BlockCrs:main_loop_noqr", policy,
      KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread) {
        auto agg = thread.league_rank();

        // size of the aggregate (number of DOFs in aggregate)
        LO aggSize                   = aggRows(agg + 1) - aggRows(agg);
        Xpetra::global_size_t offset = agg * NSDim;

        // Q = localQR(:,0)/norm
        for (LO j = 0; j < aggSize; j++) {
          LO localBlockRow = aggToRowMapLO(aggRows(agg) + j);
          LO rowStart      = localBlockRow * stride;
          for (LO r = 0; r < (LO)NSDim; r++) {
            LO localPointRow = localBlockRow * NSDim + r;
            for (LO c = 0; c < (LO)NSDim; c++) {
              values[rowStart + r * NSDim + c] = fineNSRandom(localPointRow, c);
            }
          }
        }

        // R = norm
        for (LO j = 0; j < (LO)NSDim; j++)
          coarseNS(offset + j, j) = one;
      });

  Ptentative = P_wrap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildPcoupled(RCP<Matrix> /* A */, RCP<Aggregates> /* aggregates */,
                  RCP<AmalgamationInfo> /* amalgInfo */, RCP<MultiVector> /* fineNullspace */,
                  RCP<const Map> /* coarseMap */, RCP<Matrix>& /* Ptentative */,
                  RCP<MultiVector>& /* coarseNullspace */) const {
  throw Exceptions::RuntimeError("MueLu: Construction of coupled tentative P is not implemented");
}

}  // namespace MueLu

#define MUELU_TENTATIVEPFACTORY_KOKKOS_SHORT
#endif  // MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP
