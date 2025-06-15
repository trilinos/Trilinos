// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <sstream>
#include <tuple>

#include "Xpetra_Matrix.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_decl.hpp"

#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

// #define MUELU_COALESCE_DROP_DEBUG 1

#include "MueLu_BoundaryDetection.hpp"
#include "MueLu_ClassicalDropping.hpp"
#include "MueLu_CutDrop.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_DistanceLaplacianDropping.hpp"
#include "MueLu_MatrixConstruction.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: drop tol");
  SET_VALID_ENTRY("aggregation: use ml scaling of drop tol");
  SET_VALID_ENTRY("aggregation: Dirichlet threshold");
  SET_VALID_ENTRY("aggregation: greedy Dirichlet");
  SET_VALID_ENTRY("aggregation: row sum drop tol");
  SET_VALID_ENTRY("aggregation: drop scheme");
  SET_VALID_ENTRY("aggregation: block diagonal: interleaved blocksize");
  SET_VALID_ENTRY("aggregation: distance laplacian metric");
  SET_VALID_ENTRY("aggregation: distance laplacian directional weights");
  SET_VALID_ENTRY("aggregation: dropping may create Dirichlet");
  SET_VALID_ENTRY("aggregation: distance laplacian algo");
  SET_VALID_ENTRY("aggregation: classical algo");
  SET_VALID_ENTRY("aggregation: coloring: localize color graph");

  SET_VALID_ENTRY("filtered matrix: use lumping");
  SET_VALID_ENTRY("filtered matrix: reuse graph");
  SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");

  SET_VALID_ENTRY("filtered matrix: use root stencil");
  SET_VALID_ENTRY("filtered matrix: use spread lumping");
  SET_VALID_ENTRY("filtered matrix: spread lumping diag dom growth factor");
  SET_VALID_ENTRY("filtered matrix: spread lumping diag dom cap");
  SET_VALID_ENTRY("filtered matrix: Dirichlet threshold");

#undef SET_VALID_ENTRY
  validParamList->set<bool>("lightweight wrap", true, "Experimental option for lightweight graph access");

  // "signed classical" is the Ruge-Stuben style (relative to max off-diagonal), "sign classical sa" is the signed version of the sa criterion (relative to the diagonal values)
  validParamList->getEntry("aggregation: drop scheme").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("signed classical sa", "classical", "distance laplacian", "signed classical", "block diagonal", "block diagonal classical", "block diagonal distance laplacian", "block diagonal signed classical", "block diagonal colored signed classical", "signed classical distance laplacian", "signed classical sa distance laplacian"))));
  validParamList->getEntry("aggregation: classical algo").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("default", "unscaled cut", "scaled cut", "scaled cut symmetric"))));
  validParamList->getEntry("aggregation: distance laplacian algo").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("default", "unscaled cut", "scaled cut", "scaled cut symmetric"))));
  validParamList->getEntry("aggregation: distance laplacian metric").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("unweighted", "material"))));

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase>>("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Generating factory for Coordinates");
  validParamList->set<RCP<const FactoryBase>>("BlockNumber", Teuchos::null, "Generating factory for BlockNumber");
  validParamList->set<RCP<const FactoryBase>>("Material", Teuchos::null, "Generating factory for Material");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "UnAmalgamationInfo");

  const ParameterList& pL    = GetParameterList();
  std::string algo           = pL.get<std::string>("aggregation: drop scheme");
  std::string distLaplMetric = pL.get<std::string>("aggregation: distance laplacian metric");
  if (algo == "distance laplacian" || algo == "block diagonal distance laplacian") {
    Input(currentLevel, "Coordinates");
    if (distLaplMetric == "material")
      Input(currentLevel, "Material");
  }
  if (algo == "signed classical sa")
    ;
  else if (algo.find("block diagonal") != std::string::npos) {
    Input(currentLevel, "BlockNumber");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& currentLevel) const {
  auto A = Get<RCP<Matrix>>(currentLevel, "A");
  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() % A->GetStorageBlockSize() != 0, Exceptions::RuntimeError, "A->GetFixedBlockSize() needs to be a multiple of A->GetStorageBlockSize()");
  LO blkSize = A->GetFixedBlockSize() / A->GetStorageBlockSize();

  std::tuple<GlobalOrdinal, boundary_nodes_type> results;
  if (blkSize == 1)
    results = BuildScalar(currentLevel);
  else
    results = BuildVector(currentLevel);

  if (GetVerbLevel() & Statistics1) {
    GlobalOrdinal numDropped = std::get<0>(results);
    auto boundaryNodes       = std::get<1>(results);

    GO numLocalBoundaryNodes  = 0;
    GO numGlobalBoundaryNodes = 0;

    Kokkos::parallel_reduce(
        "MueLu:CoalesceDropF:Build:bnd", range_type(0, boundaryNodes.extent(0)),
        KOKKOS_LAMBDA(const LO i, GO& n) {
          if (boundaryNodes(i))
            n++;
        },
        numLocalBoundaryNodes);

    auto comm = A->getRowMap()->getComm();
    MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);

    GO numGlobalTotal = A->getGlobalNumEntries();
    GO numGlobalDropped;
    MueLu_sumAll(comm, numDropped, numGlobalDropped);

    GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    if (numGlobalTotal != 0) {
      GetOStream(Statistics1) << "Number of dropped entries: "
                              << numGlobalDropped << "/" << numGlobalTotal
                              << " (" << 100 * Teuchos::as<double>(numGlobalDropped) / Teuchos::as<double>(numGlobalTotal) << "%)" << std::endl;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::tuple<Teuchos::RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>, Teuchos::RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetBlockNumberMVs(Level& currentLevel) const {
  RCP<LocalOrdinalVector> BlockNumber = Get<RCP<LocalOrdinalVector>>(currentLevel, "BlockNumber");
  RCP<LocalOrdinalVector> ghostedBlockNumber;
  GetOStream(Statistics1) << "Using BlockDiagonal Graph before dropping (with provided blocking)" << std::endl;

  // Ghost the column block numbers if we need to
  auto A                     = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<const Import> importer = A->getCrsGraph()->getImporter();
  if (!importer.is_null()) {
    SubFactoryMonitor m1(*this, "Block Number import", currentLevel);
    ghostedBlockNumber = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(importer->getTargetMap());
    ghostedBlockNumber->doImport(*BlockNumber, *importer, Xpetra::INSERT);
  } else {
    ghostedBlockNumber = BlockNumber;
  }
  return std::make_tuple(BlockNumber, ghostedBlockNumber);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::tuple<GlobalOrdinal, typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildScalar(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  using MatrixType        = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType         = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type = typename MatrixType::local_matrix_type;
  using local_graph_type  = typename GraphType::local_graph_type;
  using rowptr_type       = typename local_graph_type::row_map_type::non_const_type;
  using entries_type      = typename local_graph_type::entries_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using device_type       = typename Node::device_type;
  using memory_space      = typename device_type::memory_space;

  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType MT;
  const MT zero = Teuchos::ScalarTraits<MT>::zero();

  auto A = Get<RCP<Matrix>>(currentLevel, "A");

  //////////////////////////////////////////////////////////////////////
  // Process parameterlist
  const ParameterList& pL = GetParameterList();

  // Boundary detection
  const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));
  const typename STS::magnitudeType rowSumTol          = as<typename STS::magnitudeType>(pL.get<double>("aggregation: row sum drop tol"));
  const LocalOrdinal dirichletNonzeroThreshold         = 1;

  // Dropping
  const std::string algo               = pL.get<std::string>("aggregation: drop scheme");
  std::string classicalAlgoStr         = pL.get<std::string>("aggregation: classical algo");
  std::string distanceLaplacianAlgoStr = pL.get<std::string>("aggregation: distance laplacian algo");
  std::string distanceLaplacianMetric  = pL.get<std::string>("aggregation: distance laplacian metric");
  MT threshold;
  // If we're doing the ML-style halving of the drop tol at each level, we do that here.
  if (pL.get<bool>("aggregation: use ml scaling of drop tol"))
    threshold = pL.get<double>("aggregation: drop tol") / pow(2.0, currentLevel.GetLevelID());
  else
    threshold = as<MT>(pL.get<double>("aggregation: drop tol"));
  bool aggregationMayCreateDirichlet = pL.get<bool>("aggregation: dropping may create Dirichlet");

  // Fill
  const bool lumping         = pL.get<bool>("filtered matrix: use lumping");
  const bool reuseGraph      = pL.get<bool>("filtered matrix: reuse graph");
  const bool reuseEigenvalue = pL.get<bool>("filtered matrix: reuse eigenvalue");

  const bool useRootStencil   = pL.get<bool>("filtered matrix: use root stencil");
  const bool useSpreadLumping = pL.get<bool>("filtered matrix: use spread lumping");

  const MT filteringDirichletThreshold = as<MT>(pL.get<double>("filtered matrix: Dirichlet threshold"));
  TEUCHOS_ASSERT(!useRootStencil);
  TEUCHOS_ASSERT(!useSpreadLumping);

  if (algo == "classical")
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\" classical algorithm = \"" << classicalAlgoStr << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
  else if (algo == "distance laplacian")
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\" distance laplacian algorithm = \"" << distanceLaplacianAlgoStr << "\" distance laplacian metric = \"" << distanceLaplacianMetric << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
  else
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

  if (((algo == "classical") && (classicalAlgoStr.find("scaled") != std::string::npos)) || ((algo == "distance laplacian") && (distanceLaplacianAlgoStr.find("scaled") != std::string::npos)))
    TEUCHOS_TEST_FOR_EXCEPTION(threshold > 1.0, Exceptions::RuntimeError, "For cut-drop algorithms, \"aggregation: drop tol\" = " << threshold << ", needs to be <= 1.0");

  //////////////////////////////////////////////////////////////////////
  // We perform four sweeps over the rows of A:
  // Pass 1: detection of boundary nodes
  // Pass 2: diagonal extraction
  // Pass 3: drop decision for each entry and construction of the rowptr of the filtered matrix
  // Pass 4: fill of the filtered matrix
  //
  // Pass 1 and 3 apply a sequence of criteria to each row of the matrix.

  // TODO: We could merge pass 1 and 2.

  auto crsA  = toCrsMatrix(A);
  auto lclA  = crsA->getLocalMatrixDevice();
  auto range = range_type(0, lclA.numRows());

  //////////////////////////////////////////////////////////////////////
  // Pass 1: Detect boundary nodes
  //
  // The following criteria are available:
  // - BoundaryDetection::PointDirichletFunctor
  //   Marks rows as Dirichlet based on value threshold and number of off-diagonal entries
  // - BoundaryDetection::RowSumFunctor
  //   Marks rows as Dirichlet bases on row-sum criterion

  // Dirichlet nodes
  auto boundaryNodes = boundary_nodes_type("boundaryNodes", lclA.numRows());  // initialized to false
  {
    SubFactoryMonitor mBoundary(*this, "Boundary detection", currentLevel);

    // macro that applies boundary detection functors
#define MueLu_runBoundaryFunctors(...)                                          \
  {                                                                             \
    auto boundaries = BoundaryDetection::BoundaryFunctor(lclA, __VA_ARGS__);    \
    Kokkos::parallel_for("CoalesceDrop::BoundaryDetection", range, boundaries); \
  }

    auto dirichlet_detection = BoundaryDetection::PointDirichletFunctor(lclA, boundaryNodes, dirichletThreshold, dirichletNonzeroThreshold);

    if (rowSumTol <= 0.) {
      MueLu_runBoundaryFunctors(dirichlet_detection);
    } else {
      auto apply_rowsum = BoundaryDetection::RowSumFunctor(lclA, boundaryNodes, rowSumTol);
      MueLu_runBoundaryFunctors(dirichlet_detection,
                                apply_rowsum);
    }
#undef MueLu_runBoundaryFunctors
  }
  // In what follows, boundaryNodes can still still get modified if aggregationMayCreateDirichlet == true.
  // Otherwise we're now done with it now.

  //////////////////////////////////////////////////////////////////////
  // Pass 2 & 3: Diagonal extraction and determine dropping and construct
  //             rowptr of filtered matrix
  //
  // The following criteria are available:
  // - Misc::PointwiseDropBoundaryFunctor
  //   Drop all rows that have been marked as Dirichlet
  // - Misc::DropOffRankFunctor
  //   Drop all entries that are off-rank
  // - ClassicalDropping::DropFunctor
  //   Classical dropping
  // - DistanceLaplacian::DropFunctor
  //   Distance Laplacian dropping
  // - Misc::KeepDiagonalFunctor
  //   Mark diagonal as KEEP
  // - Misc::MarkSingletonFunctor
  //   Mark singletons after dropping as Dirichlet
  // - Misc::BlockDiagonalizeFunctor
  //   Drop coupling between blocks
  //
  // For the block diagonal variants we first block diagonalized and then apply "blocksize = 1" algorithms.

  // rowptr of filtered A
  auto filtered_rowptr = rowptr_type("filtered_rowptr", lclA.numRows() + 1);
  // Number of nonzeros of filtered A
  LocalOrdinal nnz_filtered = 0;
  // dropping decisions for each entry
  auto results = Kokkos::View<DecisionType*, memory_space>("results", lclA.nnz());  // initialized to UNDECIDED
  {
    SubFactoryMonitor mDropping(*this, "Dropping decisions", currentLevel);

    std::string functorLabel = "MueLu::CoalesceDrop::CountEntries";

    // macro that applied dropping functors
#if !defined(HAVE_MUELU_DEBUG)
#define MueLu_runDroppingFunctors(...)                                                                                \
  {                                                                                                                   \
    auto countingFunctor = MatrixConstruction::PointwiseCountingFunctor(lclA, results, filtered_rowptr, __VA_ARGS__); \
    Kokkos::parallel_scan(functorLabel, range, countingFunctor, nnz_filtered);                                        \
  }
#else
#define MueLu_runDroppingFunctors(...)                                                                                       \
  {                                                                                                                          \
    auto debug           = Misc::DebugFunctor(lclA, results);                                                                \
    auto countingFunctor = MatrixConstruction::PointwiseCountingFunctor(lclA, results, filtered_rowptr, __VA_ARGS__, debug); \
    Kokkos::parallel_scan(functorLabel, range, countingFunctor, nnz_filtered);                                               \
  }
#endif

    auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);

    if (threshold != zero) {
      auto preserve_diagonals          = Misc::KeepDiagonalFunctor(lclA, results);
      auto mark_singletons_as_boundary = Misc::MarkSingletonFunctor(lclA, boundaryNodes, results);

      if (algo == "classical" || algo == "block diagonal classical") {
        const auto SoC = Misc::SmoothedAggregationMeasure;

        if (algo == "block diagonal classical") {
          auto BlockNumbers      = GetBlockNumberMVs(currentLevel);
          auto block_diagonalize = Misc::BlockDiagonalizeFunctor(*A, *std::get<0>(BlockNumbers), *std::get<1>(BlockNumbers), results);

          if (classicalAlgoStr == "default") {
            auto classical_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        classical_dropping,
                                        drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);
            } else {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        classical_dropping,
                                        drop_boundaries,
                                        preserve_diagonals);
            }
          } else if (classicalAlgoStr == "unscaled cut") {
            auto comparison = CutDrop::UnscaledComparison(*A, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(block_diagonalize,
                                      drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);
          } else if (classicalAlgoStr == "scaled cut") {
            auto comparison = CutDrop::make_scaled_comparison_functor<SoC>(*A, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(block_diagonalize,
                                      drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);
          } else if (classicalAlgoStr == "scaled cut symmetric") {
            auto comparison = CutDrop::make_scaled_comparison_functor<SoC>(*A, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(block_diagonalize,
                                      drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);

            auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);

            MueLu_runDroppingFunctors(symmetrize);

          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << classicalAlgoStr << "\"");
          }
        } else {
          if (classicalAlgoStr == "default") {
            auto classical_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(classical_dropping,
                                        drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);
            } else {
              MueLu_runDroppingFunctors(classical_dropping,
                                        drop_boundaries,
                                        preserve_diagonals);
            }
          } else if (classicalAlgoStr == "unscaled cut") {
            auto comparison = CutDrop::UnscaledComparison(*A, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);
          } else if (classicalAlgoStr == "scaled cut") {
            auto comparison = CutDrop::make_scaled_comparison_functor<SoC>(*A, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);
          } else if (classicalAlgoStr == "scaled cut symmetric") {
            auto comparison = CutDrop::make_scaled_comparison_functor<SoC>(*A, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);

            auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);

            MueLu_runDroppingFunctors(symmetrize);

          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << classicalAlgoStr << "\"");
          }
        }
      } else if (algo == "signed classical" || algo == "block diagonal signed classical" || algo == "block diagonal colored signed classical") {
        const auto SoC = Misc::SignedRugeStuebenMeasure;

        auto signed_classical_rs_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

        if (algo == "block diagonal signed classical" || algo == "block diagonal colored signed classical") {
          auto BlockNumbers      = GetBlockNumberMVs(currentLevel);
          auto block_diagonalize = Misc::BlockDiagonalizeFunctor(*A, *std::get<0>(BlockNumbers), *std::get<1>(BlockNumbers), results);

          if (classicalAlgoStr == "default") {
            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        signed_classical_rs_dropping,
                                        drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);

            } else {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        signed_classical_rs_dropping,
                                        drop_boundaries,
                                        preserve_diagonals);
            }
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be default, not \"" << classicalAlgoStr << "\"");
          }
        } else {
          if (classicalAlgoStr == "default") {
            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(signed_classical_rs_dropping,
                                        drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);

            } else {
              MueLu_runDroppingFunctors(signed_classical_rs_dropping,
                                        drop_boundaries,
                                        preserve_diagonals);
            }
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be default, not \"" << classicalAlgoStr << "\"");
          }
        }
      } else if (algo == "signed classical sa") {
        const auto SoC = Misc::SignedSmoothedAggregationMeasure;
        if (classicalAlgoStr == "default") {
          auto signed_classical_sa_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

          if (aggregationMayCreateDirichlet) {
            MueLu_runDroppingFunctors(signed_classical_sa_dropping,
                                      drop_boundaries,
                                      preserve_diagonals,
                                      mark_singletons_as_boundary);

          } else {
            MueLu_runDroppingFunctors(signed_classical_sa_dropping,
                                      drop_boundaries,
                                      preserve_diagonals);
          }
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be default, not \"" << classicalAlgoStr << "\"");
        }
      } else if (algo == "distance laplacian" || algo == "block diagonal distance laplacian" || algo == "signed classical distance laplacian" || algo == "signed classical sa distance laplacian") {
        using doubleMultiVector = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>;
        auto coords             = Get<RCP<doubleMultiVector>>(currentLevel, "Coordinates");

        if (algo == "block diagonal distance laplacian") {
          const auto SoC = Misc::SmoothedAggregationMeasure;

          auto BlockNumbers      = GetBlockNumberMVs(currentLevel);
          auto block_diagonalize = Misc::BlockDiagonalizeFunctor(*A, *std::get<0>(BlockNumbers), *std::get<1>(BlockNumbers), results);

          auto dist2 = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);

          if (distanceLaplacianAlgoStr == "default") {
            auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        dist_laplacian_dropping,
                                        drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);
            } else {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        dist_laplacian_dropping,
                                        drop_boundaries,
                                        preserve_diagonals);
            }
          } else if (distanceLaplacianAlgoStr == "unscaled cut") {
            auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(block_diagonalize,
                                      drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);
          } else if (distanceLaplacianAlgoStr == "scaled cut") {
            auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(block_diagonalize,
                                      drop_boundaries,
                                      preserve_diagonals,
                                      cut_drop);
          } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
            auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
            auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

            MueLu_runDroppingFunctors(block_diagonalize,
                                      drop_boundaries,
                                      cut_drop,
                                      preserve_diagonals);

            auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);

            MueLu_runDroppingFunctors(symmetrize);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << distanceLaplacianAlgoStr << "\"");
          }
        } else if (algo == "distance laplacian") {
          const auto SoC = Misc::SmoothedAggregationMeasure;

          if (distanceLaplacianAlgoStr == "default") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2                   = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

              if (aggregationMayCreateDirichlet) {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          drop_boundaries,
                                          preserve_diagonals,
                                          mark_singletons_as_boundary);
              } else {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          drop_boundaries,
                                          preserve_diagonals);
              }
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2                   = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2                   = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              }
            }

          } else if (distanceLaplacianAlgoStr == "unscaled cut") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }
          } else if (distanceLaplacianAlgoStr == "scaled cut") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }
          } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }

            auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);

            MueLu_runDroppingFunctors(symmetrize);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << distanceLaplacianAlgoStr << "\"");
          }
        } else if (algo == "signed classical distance laplacian") {
          const auto SoC = Misc::SignedRugeStuebenMeasure;
          if (distanceLaplacianAlgoStr == "default") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2                   = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

              if (aggregationMayCreateDirichlet) {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          drop_boundaries,
                                          preserve_diagonals,
                                          mark_singletons_as_boundary);
              } else {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          drop_boundaries,
                                          preserve_diagonals);
              }
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2                   = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2                   = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              }
            }

          } else if (distanceLaplacianAlgoStr == "unscaled cut") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }
          } else if (distanceLaplacianAlgoStr == "scaled cut") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }
          } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }

            auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);

            MueLu_runDroppingFunctors(symmetrize);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << distanceLaplacianAlgoStr << "\"");
          }
        } else if (algo == "signed classical sa distance laplacian") {
          const auto SoC = Misc::SignedSmoothedAggregationMeasure;

          if (distanceLaplacianAlgoStr == "default") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2                   = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

              if (aggregationMayCreateDirichlet) {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          drop_boundaries,
                                          preserve_diagonals,
                                          mark_singletons_as_boundary);
              } else {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          drop_boundaries,
                                          preserve_diagonals);
              }
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2                   = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2                   = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(*A, threshold, dist2, results);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              }
            }

          } else if (distanceLaplacianAlgoStr == "unscaled cut") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::UnscaledDistanceLaplacianComparison(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }
          } else if (distanceLaplacianAlgoStr == "scaled cut") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }
          } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
            if (distanceLaplacianMetric == "unweighted") {
              auto dist2      = DistanceLaplacian::UnweightedDistanceFunctor(*A, coords);
              auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
              auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

              MueLu_runDroppingFunctors(drop_boundaries,
                                        preserve_diagonals,
                                        cut_drop);
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2      = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2      = DistanceLaplacian::TensorMaterialDistanceFunctor(*A, coords, material);
                auto comparison = CutDrop::make_scaled_dlap_comparison_functor<SoC>(*A, dist2, results);
                auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

                MueLu_runDroppingFunctors(drop_boundaries,
                                          preserve_diagonals,
                                          cut_drop);
              }
            }

            auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);

            MueLu_runDroppingFunctors(symmetrize);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << distanceLaplacianAlgoStr << "\"");
          }
        }
      } else if (algo == "block diagonal") {
        auto BlockNumbers      = GetBlockNumberMVs(currentLevel);
        auto block_diagonalize = Misc::BlockDiagonalizeFunctor(*A, *std::get<0>(BlockNumbers), *std::get<1>(BlockNumbers), results);

        MueLu_runDroppingFunctors(block_diagonalize);
      } else {
        TEUCHOS_ASSERT(false);
      }
    } else {
      Kokkos::deep_copy(results, KEEP);
      // FIXME: This seems inconsistent
      // MueLu_runDroppingFunctors(drop_boundaries);
      auto no_op = Misc::NoOpFunctor<LocalOrdinal>();
      MueLu_runDroppingFunctors(no_op);
    }
#undef MueLu_runDroppingFunctors
  }
  GO numDropped = lclA.nnz() - nnz_filtered;
  // We now know the number of entries of filtered A and have the final rowptr.

  //////////////////////////////////////////////////////////////////////
  // Pass 4: Create local matrix for filtered A
  //
  // Dropped entries are optionally lumped to the diagonal.

  RCP<Matrix> filteredA;
  RCP<LWGraph_kokkos> graph;
  {
    SubFactoryMonitor mFill(*this, "Filtered matrix fill", currentLevel);

    local_matrix_type lclFilteredA;
    local_graph_type lclGraph;
    if (reuseGraph) {
      filteredA    = MatrixFactory::BuildCopy(A);
      lclFilteredA = filteredA->getLocalMatrixDevice();

      auto colidx = entries_type("entries", nnz_filtered);
      lclGraph    = local_graph_type(colidx, filtered_rowptr);
    } else {
      auto colidx  = entries_type("entries", nnz_filtered);
      auto values  = values_type("values", nnz_filtered);
      lclFilteredA = local_matrix_type("filteredA",
                                       lclA.numRows(), lclA.numCols(),
                                       nnz_filtered,
                                       values, filtered_rowptr, colidx);
    }

    if (lumping) {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::PointwiseFillReuseFunctor<local_matrix_type, local_graph_type, true>(lclA, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_reuse", range, fillFunctor);
      } else {
        auto fillFunctor = MatrixConstruction::PointwiseFillNoReuseFunctor<local_matrix_type, true>(lclA, results, lclFilteredA, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_noreuse", range, fillFunctor);
      }
    } else {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::PointwiseFillReuseFunctor<local_matrix_type, local_graph_type, false>(lclA, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_reuse", range, fillFunctor);
      } else {
        auto fillFunctor = MatrixConstruction::PointwiseFillNoReuseFunctor<local_matrix_type, false>(lclA, results, lclFilteredA, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_noreuse", range, fillFunctor);
      }
    }

    if (!reuseGraph)
      filteredA = MatrixFactory::Build(lclFilteredA, A->getRowMap(), A->getColMap(), A->getDomainMap(), A->getRangeMap());
    filteredA->SetFixedBlockSize(A->GetFixedBlockSize());

    if (reuseEigenvalue) {
      // Reuse max eigenvalue from A
      // It is unclear what eigenvalue is the best for the smoothing, but we already may have
      // the D^{-1}A estimate in A, may as well use it.
      // NOTE: ML does that too
      filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
    } else {
      filteredA->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
    }

    if (!reuseGraph) {
      // Use graph of filteredA as graph.
      lclGraph = filteredA->getCrsGraph()->getLocalGraphDevice();
    }
    graph = rcp(new LWGraph_kokkos(lclGraph, filteredA->getRowMap(), filteredA->getColMap(), "amalgamated graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
  }

  LO dofsPerNode = 1;
  Set(currentLevel, "DofsPerNode", dofsPerNode);
  Set(currentLevel, "Graph", graph);
  Set(currentLevel, "A", filteredA);

  return std::make_tuple(numDropped, boundaryNodes);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::tuple<GlobalOrdinal, typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildVector(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  using MatrixType        = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType         = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type = typename MatrixType::local_matrix_type;
  using local_graph_type  = typename GraphType::local_graph_type;
  using rowptr_type       = typename local_graph_type::row_map_type::non_const_type;
  using entries_type      = typename local_graph_type::entries_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using device_type       = typename Node::device_type;
  using memory_space      = typename device_type::memory_space;

  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType MT;
  const MT zero = Teuchos::ScalarTraits<MT>::zero();

  auto A = Get<RCP<Matrix>>(currentLevel, "A");

  /* NOTE: storageblocksize (from GetStorageBlockSize()) is the size of a block in the chosen storage scheme.
     blkSize is the number of storage blocks that must kept together during the amalgamation process.

     Both of these quantities may be different than numPDEs (from GetFixedBlockSize()), but the following must always hold:

     numPDEs = blkSize * storageblocksize.

     If numPDEs==1
       Matrix is point storage (classical CRS storage).  storageblocksize=1 and  blkSize=1
       No other values makes sense.

     If numPDEs>1
       If matrix uses point storage, then storageblocksize=1  and blkSize=numPDEs.
       If matrix uses block storage, with block size of n, then storageblocksize=n, and blkSize=numPDEs/n.
       Thus far, only storageblocksize=numPDEs and blkSize=1 has been tested.
    */

  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() % A->GetStorageBlockSize() != 0, Exceptions::RuntimeError, "A->GetFixedBlockSize() needs to be a multiple of A->GetStorageBlockSize()");
  LO blkSize = A->GetFixedBlockSize() / A->GetStorageBlockSize();

  auto amalInfo = Get<RCP<AmalgamationInfo>>(currentLevel, "UnAmalgamationInfo");

  const RCP<const Map> rowMap = A->getRowMap();
  const RCP<const Map> colMap = A->getColMap();

  // build a node row map (uniqueMap = non-overlapping) and a node column map
  // (nonUniqueMap = overlapping). The arrays rowTranslation and colTranslation
  // stored in the AmalgamationInfo class container contain the local node id
  // given a local dof id. The data is calculated in the AmalgamationFactory and
  // stored in the variable "UnAmalgamationInfo" (which is of type AmalagamationInfo)
  const RCP<const Map> uniqueMap    = amalInfo->getNodeRowMap();
  const RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
  Array<LO> rowTranslationArray     = *(amalInfo->getRowTranslation());  // TAW should be transform that into a View?
  Array<LO> colTranslationArray     = *(amalInfo->getColTranslation());

  Kokkos::View<LO*, Kokkos::MemoryUnmanaged>
      rowTranslationView(rowTranslationArray.getRawPtr(), rowTranslationArray.size());
  Kokkos::View<LO*, Kokkos::MemoryUnmanaged>
      colTranslationView(colTranslationArray.getRawPtr(), colTranslationArray.size());

  // get number of local nodes
  LO numNodes = Teuchos::as<LocalOrdinal>(uniqueMap->getLocalNumElements());
  typedef typename Kokkos::View<LocalOrdinal*, typename Node::device_type> id_translation_type;
  id_translation_type rowTranslation("dofId2nodeId", rowTranslationArray.size());
  id_translation_type colTranslation("ov_dofId2nodeId", colTranslationArray.size());
  Kokkos::deep_copy(rowTranslation, rowTranslationView);
  Kokkos::deep_copy(colTranslation, colTranslationView);

  // extract striding information
  blkSize                  = A->GetFixedBlockSize();  //< the full block size (number of dofs per node in strided map)
  LocalOrdinal blkId       = -1;                      //< the block id within a strided map or -1 if it is a full block map
  LocalOrdinal blkPartSize = A->GetFixedBlockSize();  //< stores block size of part blkId (or the full block size)
  if (A->IsView("stridedMaps") == true) {
    const RCP<const Map> myMap         = A->getRowMap("stridedMaps");
    const RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
    TEUCHOS_TEST_FOR_EXCEPTION(strMap.is_null() == true, Exceptions::RuntimeError, "Map is not of type stridedMap");
    blkSize = Teuchos::as<const LocalOrdinal>(strMap->getFixedBlockSize());
    blkId   = strMap->getStridedBlockId();
    if (blkId > -1)
      blkPartSize = Teuchos::as<LocalOrdinal>(strMap->getStridingData()[blkId]);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getLocalNumElements() % blkPartSize != 0, MueLu::Exceptions::RuntimeError, "MueLu::CoalesceDropFactory: Number of local elements is " << A->getRowMap()->getLocalNumElements() << " but should be a multiple of " << blkPartSize);

  //////////////////////////////////////////////////////////////////////
  // Process parameterlist
  const ParameterList& pL = GetParameterList();

  // Boundary detection
  const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));
  const typename STS::magnitudeType rowSumTol          = as<typename STS::magnitudeType>(pL.get<double>("aggregation: row sum drop tol"));
  const LocalOrdinal dirichletNonzeroThreshold         = 1;
  const bool useGreedyDirichlet                        = pL.get<bool>("aggregation: greedy Dirichlet");
  TEUCHOS_TEST_FOR_EXCEPTION(rowSumTol > zero, MueLu::Exceptions::RuntimeError, "MueLu::CoalesceDropFactory: RowSum is not implemented for vectorial problems.");

  // Dropping
  const std::string algo               = pL.get<std::string>("aggregation: drop scheme");
  std::string classicalAlgoStr         = pL.get<std::string>("aggregation: classical algo");
  std::string distanceLaplacianAlgoStr = pL.get<std::string>("aggregation: distance laplacian algo");
  std::string distanceLaplacianMetric  = pL.get<std::string>("aggregation: distance laplacian metric");
  MT threshold;
  // If we're doing the ML-style halving of the drop tol at each level, we do that here.
  if (pL.get<bool>("aggregation: use ml scaling of drop tol"))
    threshold = pL.get<double>("aggregation: drop tol") / pow(2.0, currentLevel.GetLevelID());
  else
    threshold = as<MT>(pL.get<double>("aggregation: drop tol"));
  bool aggregationMayCreateDirichlet = pL.get<bool>("aggregation: dropping may create Dirichlet");

  // Fill
  const bool lumping         = pL.get<bool>("filtered matrix: use lumping");
  const bool reuseGraph      = pL.get<bool>("filtered matrix: reuse graph");
  const bool reuseEigenvalue = pL.get<bool>("filtered matrix: reuse eigenvalue");

  const bool useRootStencil   = pL.get<bool>("filtered matrix: use root stencil");
  const bool useSpreadLumping = pL.get<bool>("filtered matrix: use spread lumping");

  const MT filteringDirichletThreshold = as<MT>(pL.get<double>("filtered matrix: Dirichlet threshold"));

  TEUCHOS_ASSERT(!useRootStencil);
  TEUCHOS_ASSERT(!useSpreadLumping);

  if (algo == "classical") {
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\" classical algorithm = \"" << classicalAlgoStr << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
  } else if (algo == "distance laplacian") {
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\" distance laplacian algorithm = \"" << distanceLaplacianAlgoStr << "\" distance laplacian metric = \"" << distanceLaplacianMetric << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
  } else
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

  //////////////////////////////////////////////////////////////////////
  // We perform four sweeps over the rows of A:
  // Pass 1: detection of boundary nodes
  // Pass 2: diagonal extraction
  // Pass 3: drop decision for each entry and construction of the rowptr of the filtered matrix
  // Pass 4: fill of the filtered matrix
  //
  // Pass 1 and 3 apply a sequence of criteria to each row of the matrix.

  // TODO: We could merge pass 1 and 2.

  auto crsA  = toCrsMatrix(A);
  auto lclA  = crsA->getLocalMatrixDevice();
  auto range = range_type(0, numNodes);

  //////////////////////////////////////////////////////////////////////
  // Pass 1: Detect boundary nodes
  //
  // The following criteria are available:
  // - BoundaryDetection::VectorDirichletFunctor
  //   Marks rows as Dirichlet based on value threshold and number of off-diagonal entries

  // Dirichlet nodes
  auto boundaryNodes = boundary_nodes_type("boundaryNodes", numNodes);  // initialized to false
  {
    SubFactoryMonitor mBoundary(*this, "Boundary detection", currentLevel);

#define MueLu_runBoundaryFunctors(...)                                          \
  {                                                                             \
    auto boundaries = BoundaryDetection::BoundaryFunctor(lclA, __VA_ARGS__);    \
    Kokkos::parallel_for("CoalesceDrop::BoundaryDetection", range, boundaries); \
  }

    if (useGreedyDirichlet) {
      auto dirichlet_detection = BoundaryDetection::VectorDirichletFunctor<local_matrix_type, true>(lclA, blkPartSize, boundaryNodes, dirichletThreshold, dirichletNonzeroThreshold);
      MueLu_runBoundaryFunctors(dirichlet_detection);
    } else {
      auto dirichlet_detection = BoundaryDetection::VectorDirichletFunctor<local_matrix_type, false>(lclA, blkPartSize, boundaryNodes, dirichletThreshold, dirichletNonzeroThreshold);
      MueLu_runBoundaryFunctors(dirichlet_detection);
    }
#undef MueLu_runBoundaryFunctors
  }
  // In what follows, boundaryNodes can still still get modified if aggregationMayCreateDirichlet == true.
  // Otherwise we're now done with it now.

  //////////////////////////////////////////////////////////////////////
  // Pass 2 & 3: Diagonal extraction and determine dropping and construct
  //             rowptr of filtered matrix
  //
  // The following criteria are available:
  // - Misc::VectorDropBoundaryFunctor
  //   Drop all rows that have been marked as Dirichlet
  // - Misc::DropOffRankFunctor
  //   Drop all entries that are off-rank
  // - ClassicalDropping::DropFunctor
  //   Classical dropping
  // - DistanceLaplacian::VectorDropFunctor
  //   Distance Laplacian dropping
  // - Misc::KeepDiagonalFunctor
  //   Mark diagonal as KEEP
  // - Misc::MarkSingletonFunctor
  //   Mark singletons after dropping as Dirichlet

  // rowptr of filtered A
  auto filtered_rowptr = rowptr_type("rowptr", lclA.numRows() + 1);
  auto graph_rowptr    = rowptr_type("rowptr", numNodes + 1);
  // Number of nonzeros of filtered A and graph
  Kokkos::pair<LocalOrdinal, LocalOrdinal> nnz = {0, 0};

  // dropping decisions for each entry
  auto results = Kokkos::View<DecisionType*, memory_space>("results", lclA.nnz());  // initialized to UNDECIDED
  {
    SubFactoryMonitor mDropping(*this, "Dropping decisions", currentLevel);

    std::string functorLabel = "MueLu::CoalesceDrop::CountEntries";

#if !defined(HAVE_MUELU_DEBUG)
#define MueLu_runDroppingFunctors(...)                                                                                                                        \
  {                                                                                                                                                           \
    auto countingFunctor = MatrixConstruction::VectorCountingFunctor(lclA, blkPartSize, colTranslation, results, filtered_rowptr, graph_rowptr, __VA_ARGS__); \
    Kokkos::parallel_scan(functorLabel, range, countingFunctor, nnz);                                                                                         \
  }
#else
#define MueLu_runDroppingFunctors(...)                                                                                                                               \
  {                                                                                                                                                                  \
    auto debug           = Misc::DebugFunctor(lclA, results);                                                                                                        \
    auto countingFunctor = MatrixConstruction::VectorCountingFunctor(lclA, blkPartSize, colTranslation, results, filtered_rowptr, graph_rowptr, __VA_ARGS__, debug); \
    Kokkos::parallel_scan(functorLabel, range, countingFunctor, nnz);                                                                                                \
  }
#endif

    auto drop_boundaries = Misc::VectorDropBoundaryFunctor(lclA, rowTranslation, boundaryNodes, results);

    if (threshold != zero) {
      auto preserve_diagonals          = Misc::KeepDiagonalFunctor(lclA, results);
      auto mark_singletons_as_boundary = Misc::MarkSingletonVectorFunctor(lclA, rowTranslation, boundaryNodes, results);

      if (algo == "classical" || algo == "block diagonal classical") {
        const auto SoC = Misc::SmoothedAggregationMeasure;

        if (classicalAlgoStr == "default") {
          auto classical_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

          if (algo == "block diagonal classical") {
            RCP<LocalOrdinalVector> BlockNumber = Get<RCP<LocalOrdinalVector>>(currentLevel, "BlockNumber");
            auto block_diagonalize              = Misc::BlockDiagonalizeVectorFunctor(*A, *BlockNumber, nonUniqueMap, results, rowTranslation, colTranslation);

            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        classical_dropping,
                                        // drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);
            } else {
              MueLu_runDroppingFunctors(block_diagonalize,
                                        classical_dropping,
                                        // drop_boundaries,
                                        preserve_diagonals);
            }
          } else {
            if (aggregationMayCreateDirichlet) {
              MueLu_runDroppingFunctors(classical_dropping,
                                        // drop_boundaries,
                                        preserve_diagonals,
                                        mark_singletons_as_boundary);
            } else {
              MueLu_runDroppingFunctors(classical_dropping,
                                        // drop_boundaries,
                                        preserve_diagonals);
            }
          }
        } else if (classicalAlgoStr == "unscaled cut") {
          TEUCHOS_ASSERT(false);
        } else if (classicalAlgoStr == "scaled cut") {
          TEUCHOS_ASSERT(false);
        } else if (classicalAlgoStr == "scaled cut symmetric") {
          TEUCHOS_ASSERT(false);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << classicalAlgoStr << "\"");
        }
      } else if (algo == "signed classical" || algo == "block diagonal colored signed classical" || algo == "block diagonal signed classical") {
        const auto SoC = Misc::SignedRugeStuebenMeasure;

        auto signed_classical_rs_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

        if (algo == "block diagonal signed classical" || algo == "block diagonal colored signed classical") {
          auto BlockNumbers      = GetBlockNumberMVs(currentLevel);
          auto block_diagonalize = Misc::BlockDiagonalizeFunctor(*A, *std::get<0>(BlockNumbers), *std::get<1>(BlockNumbers), results);

          if (aggregationMayCreateDirichlet) {
            MueLu_runDroppingFunctors(block_diagonalize,
                                      signed_classical_rs_dropping,
                                      // drop_boundaries,
                                      preserve_diagonals,
                                      mark_singletons_as_boundary);

          } else {
            MueLu_runDroppingFunctors(block_diagonalize,
                                      signed_classical_rs_dropping,
                                      // drop_boundaries,
                                      preserve_diagonals);
          }
        } else {
          if (aggregationMayCreateDirichlet) {
            MueLu_runDroppingFunctors(signed_classical_rs_dropping,
                                      // drop_boundaries,
                                      preserve_diagonals,
                                      mark_singletons_as_boundary);

          } else {
            MueLu_runDroppingFunctors(signed_classical_rs_dropping,
                                      // drop_boundaries,
                                      preserve_diagonals);
          }
        }
      } else if (algo == "signed classical sa") {
        const auto SoC = Misc::SignedSmoothedAggregationMeasure;

        auto signed_classical_sa_dropping = ClassicalDropping::make_drop_functor<SoC>(*A, threshold, results);

        if (aggregationMayCreateDirichlet) {
          MueLu_runDroppingFunctors(signed_classical_sa_dropping,
                                    // drop_boundaries,
                                    preserve_diagonals,
                                    mark_singletons_as_boundary);

        } else {
          MueLu_runDroppingFunctors(signed_classical_sa_dropping,
                                    // drop_boundaries,
                                    preserve_diagonals);
        }
      } else if (algo == "distance laplacian" || algo == "block diagonal distance laplacian") {
        using doubleMultiVector = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>;
        auto coords             = Get<RCP<doubleMultiVector>>(currentLevel, "Coordinates");

        bool use_block_algorithm   = (algo == "block diagonal classical" || algo == "block diagonal distance laplacian");
        Array<double> dlap_weights = pL.get<Array<double>>("aggregation: distance laplacian directional weights");
        enum { NO_WEIGHTS = 0,
               SINGLE_WEIGHTS,
               BLOCK_WEIGHTS };
        int use_dlap_weights = NO_WEIGHTS;
        if (algo == "distance laplacian") {
          LO dim = (LO)coords->getNumVectors();
          // If anything isn't 1.0 we need to turn on the weighting
          bool non_unity = false;
          for (LO i = 0; !non_unity && i < (LO)dlap_weights.size(); i++) {
            if (dlap_weights[i] != 1.0) {
              non_unity = true;
            }
          }
          if (non_unity) {
            LO blocksize = use_block_algorithm ? as<LO>(pL.get<int>("aggregation: block diagonal: interleaved blocksize")) : 1;
            if ((LO)dlap_weights.size() == dim)
              use_dlap_weights = SINGLE_WEIGHTS;
            else if ((LO)dlap_weights.size() == blocksize * dim)
              use_dlap_weights = BLOCK_WEIGHTS;
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError,
                                         "length of 'aggregation: distance laplacian directional weights' must equal the coordinate dimension OR the coordinate dimension times the blocksize");
            }
            if (GetVerbLevel() & Statistics1)
              GetOStream(Statistics1) << "Using distance laplacian weights: " << dlap_weights << std::endl;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(use_dlap_weights != NO_WEIGHTS, Exceptions::RuntimeError, "Only the NO_WEIGHTS option is implemented for distance laplacian ");

        RCP<Matrix> mergedA;
        {
          // Construct merged A.

          auto merged_rowptr      = rowptr_type("rowptr", numNodes + 1);
          LocalOrdinal nnz_merged = 0;

          auto functor = MatrixConstruction::MergeCountFunctor(lclA, blkPartSize, colTranslation, merged_rowptr);
          Kokkos::parallel_scan("MergeCount", range, functor, nnz_merged);

          local_graph_type lclMergedGraph;
          auto colidx_merged = entries_type("entries", nnz_merged);
          auto values_merged = values_type("values", nnz_merged);

          local_matrix_type lclMergedA = local_matrix_type("mergedA",
                                                           numNodes, nonUniqueMap->getLocalNumElements(),
                                                           nnz_merged,
                                                           values_merged, merged_rowptr, colidx_merged);

          auto fillFunctor = MatrixConstruction::MergeFillFunctor<local_matrix_type>(lclA, blkSize, colTranslation, lclMergedA);
          Kokkos::parallel_for("MueLu::CoalesceDrop::MergeFill", range, fillFunctor);

          mergedA = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclMergedA, uniqueMap, nonUniqueMap, uniqueMap, uniqueMap);
        }

        if (algo == "block diagonal distance laplacian") {
          RCP<LocalOrdinalVector> BlockNumber = Get<RCP<LocalOrdinalVector>>(currentLevel, "BlockNumber");
          auto block_diagonalize              = Misc::BlockDiagonalizeVectorFunctor(*A, *BlockNumber, nonUniqueMap, results, rowTranslation, colTranslation);

          if (distanceLaplacianAlgoStr == "default") {
            const auto SoC = Misc::SmoothedAggregationMeasure;

            if (distanceLaplacianMetric == "unweighted") {
              auto dist2                   = DistanceLaplacian::UnweightedDistanceFunctor(*mergedA, coords);
              auto dist_laplacian_dropping = DistanceLaplacian::make_vector_drop_functor<SoC>(*A, *mergedA, threshold, dist2, results, rowTranslation, colTranslation);

              if (aggregationMayCreateDirichlet) {
                MueLu_runDroppingFunctors(block_diagonalize,
                                          dist_laplacian_dropping,
                                          // drop_boundaries,
                                          preserve_diagonals,
                                          mark_singletons_as_boundary);
              } else {
                MueLu_runDroppingFunctors(block_diagonalize,
                                          dist_laplacian_dropping,
                                          // drop_boundaries,
                                          preserve_diagonals);
              }
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2                   = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_vector_drop_functor<SoC>(*A, *mergedA, threshold, dist2, results, rowTranslation, colTranslation);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(block_diagonalize,
                                            dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(block_diagonalize,
                                            dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2                   = DistanceLaplacian::TensorMaterialDistanceFunctor(*mergedA, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_vector_drop_functor<SoC>(*A, *mergedA, threshold, dist2, results, rowTranslation, colTranslation);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(block_diagonalize,
                                            dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(block_diagonalize,
                                            dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              }
            }
          } else if (distanceLaplacianAlgoStr == "unscaled cut") {
            TEUCHOS_ASSERT(false);
          } else if (distanceLaplacianAlgoStr == "scaled cut") {
            TEUCHOS_ASSERT(false);
          } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
            TEUCHOS_ASSERT(false);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << distanceLaplacianAlgoStr << "\"");
          }
        } else {
          if (distanceLaplacianAlgoStr == "default") {
            const auto SoC = Misc::SmoothedAggregationMeasure;

            if (distanceLaplacianMetric == "unweighted") {
              auto dist2                   = DistanceLaplacian::UnweightedDistanceFunctor(*mergedA, coords);
              auto dist_laplacian_dropping = DistanceLaplacian::make_vector_drop_functor<SoC>(*A, *mergedA, threshold, dist2, results, rowTranslation, colTranslation);

              if (aggregationMayCreateDirichlet) {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          // drop_boundaries,
                                          preserve_diagonals,
                                          mark_singletons_as_boundary);
              } else {
                MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                          // drop_boundaries,
                                          preserve_diagonals);
              }
            } else if (distanceLaplacianMetric == "material") {
              auto material = Get<RCP<MultiVector>>(currentLevel, "Material");
              if (material->getNumVectors() == 1) {
                GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;

                auto dist2                   = DistanceLaplacian::ScalarMaterialDistanceFunctor(*A, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_vector_drop_functor<SoC>(*A, *mergedA, threshold, dist2, results, rowTranslation, colTranslation);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              } else {
                TEUCHOS_TEST_FOR_EXCEPTION(coords->getNumVectors() * coords->getNumVectors() != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");

                {
                  std::stringstream ss;
                  ss << "material tensor mean =" << std::endl;
                  size_t k = 0;
                  for (size_t i = 0; i < coords->getNumVectors(); ++i) {
                    ss << "   ";
                    for (size_t j = 0; j < coords->getNumVectors(); ++j) {
                      ss << material->getVector(k)->meanValue() << " ";
                      ++k;
                    }
                    ss << std::endl;
                  }
                  GetOStream(Runtime0) << ss.str();
                }

                auto dist2                   = DistanceLaplacian::TensorMaterialDistanceFunctor(*mergedA, coords, material);
                auto dist_laplacian_dropping = DistanceLaplacian::make_vector_drop_functor<SoC>(*A, *mergedA, threshold, dist2, results, rowTranslation, colTranslation);

                if (aggregationMayCreateDirichlet) {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals,
                                            mark_singletons_as_boundary);
                } else {
                  MueLu_runDroppingFunctors(dist_laplacian_dropping,
                                            drop_boundaries,
                                            preserve_diagonals);
                }
              }
            }
          } else if (distanceLaplacianAlgoStr == "unscaled cut") {
            TEUCHOS_ASSERT(false);
          } else if (distanceLaplacianAlgoStr == "scaled cut") {
            TEUCHOS_ASSERT(false);
          } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
            TEUCHOS_ASSERT(false);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut|scaled cut symmetric), not \"" << distanceLaplacianAlgoStr << "\"");
          }
        }

      } else if (algo == "block diagonal") {
        RCP<LocalOrdinalVector> BlockNumber = Get<RCP<LocalOrdinalVector>>(currentLevel, "BlockNumber");
        auto block_diagonalize              = Misc::BlockDiagonalizeVectorFunctor(*A, *BlockNumber, nonUniqueMap, results, rowTranslation, colTranslation);

        std::cout << "fact line1870" << std::endl;

        MueLu_runDroppingFunctors(block_diagonalize);
      } else {
        TEUCHOS_ASSERT(false);
      }
    } else {
      Kokkos::deep_copy(results, KEEP);
      // MueLu_runDroppingFunctors(drop_boundaries);
      auto no_op = Misc::NoOpFunctor<LocalOrdinal>();
      MueLu_runDroppingFunctors(no_op);
    }
#undef MueLu_runDroppingFunctors
  }
  LocalOrdinal nnz_filtered = nnz.first;
  LocalOrdinal nnz_graph    = nnz.second;
  GO numTotal               = lclA.nnz();
  GO numDropped             = numTotal - nnz_filtered;
  // We now know the number of entries of filtered A and have the final rowptr.

  //////////////////////////////////////////////////////////////////////
  // Pass 4: Create local matrix for filtered A
  //
  // Dropped entries are optionally lumped to the diagonal.

  RCP<Matrix> filteredA;
  RCP<LWGraph_kokkos> graph;
  {
    SubFactoryMonitor mFill(*this, "Filtered matrix fill", currentLevel);

    local_matrix_type lclFilteredA;
    if (reuseGraph) {
      lclFilteredA = local_matrix_type("filteredA", lclA.graph, lclA.numCols());
    } else {
      auto colidx  = entries_type("entries", nnz_filtered);
      auto values  = values_type("values", nnz_filtered);
      lclFilteredA = local_matrix_type("filteredA",
                                       lclA.numRows(), lclA.numCols(),
                                       nnz_filtered,
                                       values, filtered_rowptr, colidx);
    }

    local_graph_type lclGraph;
    {
      auto colidx = entries_type("entries", nnz_graph);
      lclGraph    = local_graph_type(colidx, graph_rowptr);
    }

    if (lumping) {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, true, true>(lclA, blkPartSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_reuse", range, fillFunctor);
      } else {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, true, false>(lclA, blkPartSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_noreuse", range, fillFunctor);
      }
    } else {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, false, true>(lclA, blkSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_reuse", range, fillFunctor);
      } else {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, false, false>(lclA, blkSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_noreuse", range, fillFunctor);
      }
    }

    filteredA = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclFilteredA, A->getRowMap(), A->getColMap(), A->getDomainMap(), A->getRangeMap());
    filteredA->SetFixedBlockSize(blkSize);

    if (reuseEigenvalue) {
      // Reuse max eigenvalue from A
      // It is unclear what eigenvalue is the best for the smoothing, but we already may have
      // the D^{-1}A estimate in A, may as well use it.
      // NOTE: ML does that too
      filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
    } else {
      filteredA->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
    }

    graph = rcp(new LWGraph_kokkos(lclGraph, uniqueMap, nonUniqueMap, "amalgamated graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
  }

  LO dofsPerNode = blkSize;

  Set(currentLevel, "DofsPerNode", dofsPerNode);
  Set(currentLevel, "Graph", graph);
  Set(currentLevel, "A", filteredA);

  return std::make_tuple(numDropped, boundaryNodes);
}

}  // namespace MueLu
#endif  // MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
