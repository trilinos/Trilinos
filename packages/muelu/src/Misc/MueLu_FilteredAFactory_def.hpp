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
#ifndef MUELU_FILTEREDAFACTORY_DEF_HPP
#define MUELU_FILTEREDAFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_FilteredAFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"

// Variable to enable lots of debug output
#define MUELU_FILTEREDAFACTORY_LOTS_OF_PRINTING 0

namespace MueLu {

template <class T>
void sort_and_unique(T& array) {
  std::sort(array.begin(), array.end());
  std::unique(array.begin(), array.end());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("filtered matrix: use lumping");
  SET_VALID_ENTRY("filtered matrix: reuse graph");
  SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");
  SET_VALID_ENTRY("filtered matrix: use root stencil");
  SET_VALID_ENTRY("filtered matrix: use spread lumping");
  SET_VALID_ENTRY("filtered matrix: spread lumping diag dom growth factor");
  SET_VALID_ENTRY("filtered matrix: spread lumping diag dom cap");
  SET_VALID_ENTRY("filtered matrix: Dirichlet threshold");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A used for filtering");
  validParamList->set<RCP<const FactoryBase> >("Graph", Teuchos::null, "Generating factory for coalesced filtered graph");
  validParamList->set<RCP<const FactoryBase> >("Filtering", Teuchos::null, "Generating factory for filtering boolean");

  // Only need these for the "use root stencil" option
  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Generating factory of the aggregates");
  validParamList->set<RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "Filtering");
  Input(currentLevel, "Graph");
  const ParameterList& pL = GetParameterList();
  if (pL.isParameter("filtered matrix: use root stencil") && pL.get<bool>("filtered matrix: use root stencil") == true) {
    Input(currentLevel, "Aggregates");
    Input(currentLevel, "UnAmalgamationInfo");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Matrix filtering", currentLevel);

  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");
  if (Get<bool>(currentLevel, "Filtering") == false) {
    GetOStream(Runtime0) << "Filtered matrix is not being constructed as no filtering is being done" << std::endl;
    Set(currentLevel, "A", A);
    return;
  }

  const ParameterList& pL = GetParameterList();
  bool lumping            = pL.get<bool>("filtered matrix: use lumping");
  if (lumping)
    GetOStream(Runtime0) << "Lumping dropped entries" << std::endl;

  bool use_spread_lumping = pL.get<bool>("filtered matrix: use spread lumping");
  if (use_spread_lumping && (!lumping))
    throw std::runtime_error("Must also request 'filtered matrix: use lumping' in order to use spread lumping");

  if (use_spread_lumping) {
    GetOStream(Runtime0) << "using spread lumping " << std::endl;
  }

  double DdomAllowGrowthRate = 1.1;
  double DdomCap             = 2.0;
  if (use_spread_lumping) {
    DdomAllowGrowthRate = pL.get<double>("filtered matrix: spread lumping diag dom growth factor");
    DdomCap             = pL.get<double>("filtered matrix: spread lumping diag dom cap");
  }
  bool use_root_stencil = lumping && pL.get<bool>("filtered matrix: use root stencil");
  if (use_root_stencil)
    GetOStream(Runtime0) << "Using root stencil for dropping" << std::endl;
  double dirichlet_threshold = pL.get<double>("filtered matrix: Dirichlet threshold");
  if (dirichlet_threshold >= 0.0)
    GetOStream(Runtime0) << "Filtering Dirichlet threshold of " << dirichlet_threshold << std::endl;

  if (use_root_stencil || pL.get<bool>("filtered matrix: reuse graph"))
    GetOStream(Runtime0) << "Reusing graph" << std::endl;
  else
    GetOStream(Runtime0) << "Generating new graph" << std::endl;

  RCP<GraphBase> G = Get<RCP<GraphBase> >(currentLevel, "Graph");
  if (MUELU_FILTEREDAFACTORY_LOTS_OF_PRINTING) {
    FILE* f         = fopen("graph.dat", "w");
    size_t numGRows = G->GetNodeNumVertices();
    for (size_t i = 0; i < numGRows; i++) {
      // Set up filtering array
      ArrayView<const LO> indsG = G->getNeighborVertices(i);
      for (size_t j = 0; j < (size_t)indsG.size(); j++) {
        fprintf(f, "%d %d 1.0\n", (int)i, (int)indsG[j]);
      }
    }
    fclose(f);
  }

  RCP<ParameterList> fillCompleteParams(new ParameterList);
  fillCompleteParams->set("No Nonlocal Changes", true);

  RCP<Matrix> filteredA;
  if (use_root_stencil) {
    filteredA = MatrixFactory::Build(A->getCrsGraph());
    filteredA->fillComplete(fillCompleteParams);
    filteredA->resumeFill();
    BuildNewUsingRootStencil(*A, *G, dirichlet_threshold, currentLevel, *filteredA, use_spread_lumping, DdomAllowGrowthRate, DdomCap);
    filteredA->fillComplete(fillCompleteParams);

  } else if (pL.get<bool>("filtered matrix: reuse graph")) {
    filteredA = MatrixFactory::Build(A->getCrsGraph());
    filteredA->resumeFill();
    BuildReuse(*A, *G, (lumping != use_spread_lumping), dirichlet_threshold, *filteredA);
    // only lump inside BuildReuse if lumping is true and use_spread_lumping is false
    // note: they use_spread_lumping cannot be true if lumping is false

    if (use_spread_lumping) ExperimentalLumping(*A, *filteredA, DdomAllowGrowthRate, DdomCap);
    filteredA->fillComplete(fillCompleteParams);

  } else {
    filteredA = MatrixFactory::Build(A->getRowMap(), A->getColMap(), A->getLocalMaxNumRowEntries());
    BuildNew(*A, *G, (lumping != use_spread_lumping), dirichlet_threshold, *filteredA);
    // only lump inside BuildNew if lumping is true and use_spread_lumping is false
    // note: they use_spread_lumping cannot be true if lumping is false
    if (use_spread_lumping) ExperimentalLumping(*A, *filteredA, DdomAllowGrowthRate, DdomCap);
    filteredA->fillComplete(A->getDomainMap(), A->getRangeMap(), fillCompleteParams);
  }

  if (MUELU_FILTEREDAFACTORY_LOTS_OF_PRINTING) {
    Xpetra::IO<SC, LO, GO, NO>::Write("filteredA.dat", *filteredA);

    // original filtered A  and actual A
    Xpetra::IO<SC, LO, GO, NO>::Write("A.dat", *A);
    RCP<Matrix> origFilteredA = MatrixFactory::Build(A->getRowMap(), A->getColMap(), A->getLocalMaxNumRowEntries());
    BuildNew(*A, *G, lumping, dirichlet_threshold, *origFilteredA);
    if (use_spread_lumping) ExperimentalLumping(*A, *origFilteredA, DdomAllowGrowthRate, DdomCap);
    origFilteredA->fillComplete(A->getDomainMap(), A->getRangeMap(), fillCompleteParams);
    Xpetra::IO<SC, LO, GO, NO>::Write("origFilteredA.dat", *origFilteredA);
  }

  filteredA->SetFixedBlockSize(A->GetFixedBlockSize());

  if (pL.get<bool>("filtered matrix: reuse eigenvalue")) {
    // Reuse max eigenvalue from A
    // It is unclear what eigenvalue is the best for the smoothing, but we already may have
    // the D^{-1}A estimate in A, may as well use it.
    // NOTE: ML does that too
    filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
  }

  Set(currentLevel, "A", filteredA);
}

// Epetra's API allows direct access to row array.
// Tpetra's API does not, providing only ArrayView<const .>
// But in most situations we are currently interested in, it is safe to assume
// that the view is to the actual data. So this macro directs the code to do
// const_cast, and modify the entries directly. This allows us to avoid
// replaceLocalValues() call which is quite expensive due to all the searches.
//#define ASSUME_DIRECT_ACCESS_TO_ROW // See github issue 10883#issuecomment-1256676340

// Both Epetra and Tpetra matrix-matrix multiply use the following trick:
// if an entry of the left matrix is zero, it does not compute or store the
// zero value.
//
// This trick allows us to bypass constructing a new matrix. Instead, we
// make a deep copy of the original one, and fill it in with zeros, which
// are ignored during the prolongator smoothing.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildReuse(const Matrix& A, const GraphBase& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const {
  using TST = typename Teuchos::ScalarTraits<SC>;
  SC zero   = TST::zero();

  size_t blkSize = A.GetFixedBlockSize();

  ArrayView<const LO> inds;
  ArrayView<const SC> valsA;
#ifdef ASSUME_DIRECT_ACCESS_TO_ROW
  ArrayView<SC> vals;
#else
  Array<SC> vals;
#endif

  Array<char> filter(std::max(blkSize * G.GetImportMap()->getLocalNumElements(),
                              A.getColMap()->getLocalNumElements()),
                     0);

  size_t numGRows = G.GetNodeNumVertices();
  for (size_t i = 0; i < numGRows; i++) {
    // Set up filtering array
    ArrayView<const LO> indsG = G.getNeighborVertices(i);
    for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
      for (size_t k = 0; k < blkSize; k++)
        filter[indsG[j] * blkSize + k] = 1;

    for (size_t k = 0; k < blkSize; k++) {
      LO row = i * blkSize + k;

      A.getLocalRowView(row, inds, valsA);

      size_t nnz = inds.size();
      if (nnz == 0)
        continue;

#ifdef ASSUME_DIRECT_ACCESS_TO_ROW
      // Transform ArrayView<const SC> into ArrayView<SC>
      ArrayView<const SC> vals1;
      filteredA.getLocalRowView(row, inds, vals1);
      vals = ArrayView<SC>(const_cast<SC*>(vals1.getRawPtr()), nnz);

      memcpy(vals.getRawPtr(), valsA.getRawPtr(), nnz * sizeof(SC));
#else
      vals = Array<SC>(valsA);
#endif

      SC ZERO = Teuchos::ScalarTraits<SC>::zero();
      //      SC ONE = Teuchos::ScalarTraits<SC>::one();
      SC A_rowsum = ZERO, F_rowsum = ZERO;
      for (LO l = 0; l < (LO)inds.size(); l++)
        A_rowsum += valsA[l];

      if (lumping == false) {
        for (size_t j = 0; j < nnz; j++)
          if (!filter[inds[j]])
            vals[j] = zero;

      } else {
        LO diagIndex = -1;
        SC diagExtra = zero;

        for (size_t j = 0; j < nnz; j++) {
          if (filter[inds[j]]) {
            if (inds[j] == row) {
              // Remember diagonal position
              diagIndex = j;
            }
            continue;
          }

          diagExtra += vals[j];

          vals[j] = zero;
        }

        // Lump dropped entries
        // NOTE
        //  * Does it make sense to lump for elasticity?
        //  * Is it different for diffusion and elasticity?
        // SC diagA = ZERO;
        if (diagIndex != -1) {
          // diagA = vals[diagIndex];
          vals[diagIndex] += diagExtra;
          if (dirichletThresh >= 0.0 && TST::real(vals[diagIndex]) <= dirichletThresh) {
            //            printf("WARNING: row %d diag(Afiltered) = %8.2e diag(A)=%8.2e\n",row,vals[diagIndex],diagA);
            for (LO l = 0; l < (LO)nnz; l++)
              F_rowsum += vals[l];
            //            printf("       : A rowsum = %8.2e F rowsum = %8.2e\n",A_rowsum,F_rowsum);
            vals[diagIndex] = TST::one();
          }
        }
      }

#ifndef ASSUME_DIRECT_ACCESS_TO_ROW
      // Because we used a column map in the construction of the matrix
      // we can just use insertLocalValues here instead of insertGlobalValues
      filteredA.replaceLocalValues(row, inds, vals);
#endif
    }

    // Reset filtering array
    for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
      for (size_t k = 0; k < blkSize; k++)
        filter[indsG[j] * blkSize + k] = 0;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildNew(const Matrix& A, const GraphBase& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const {
  using TST = typename Teuchos::ScalarTraits<SC>;
  SC zero   = Teuchos::ScalarTraits<SC>::zero();

  size_t blkSize = A.GetFixedBlockSize();

  ArrayView<const LO> indsA;
  ArrayView<const SC> valsA;
  Array<LO> inds;
  Array<SC> vals;

  Array<char> filter(blkSize * G.GetImportMap()->getLocalNumElements(), 0);

  size_t numGRows = G.GetNodeNumVertices();
  for (size_t i = 0; i < numGRows; i++) {
    // Set up filtering array
    ArrayView<const LO> indsG = G.getNeighborVertices(i);
    for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
      for (size_t k = 0; k < blkSize; k++)
        filter[indsG[j] * blkSize + k] = 1;

    for (size_t k = 0; k < blkSize; k++) {
      LO row = i * blkSize + k;

      A.getLocalRowView(row, indsA, valsA);

      size_t nnz = indsA.size();
      if (nnz == 0)
        continue;

      inds.resize(indsA.size());
      vals.resize(valsA.size());

      size_t numInds = 0;
      if (lumping == false) {
        for (size_t j = 0; j < nnz; j++)
          if (filter[indsA[j]]) {
            inds[numInds] = indsA[j];
            vals[numInds] = valsA[j];
            numInds++;
          }

      } else {
        LO diagIndex = -1;
        SC diagExtra = zero;

        for (size_t j = 0; j < nnz; j++) {
          if (filter[indsA[j]]) {
            inds[numInds] = indsA[j];
            vals[numInds] = valsA[j];

            // Remember diagonal position
            if (inds[numInds] == row)
              diagIndex = numInds;

            numInds++;

          } else {
            diagExtra += valsA[j];
          }
        }

        // Lump dropped entries
        // NOTE
        //  * Does it make sense to lump for elasticity?
        //  * Is it different for diffusion and elasticity?
        if (diagIndex != -1) {
          vals[diagIndex] += diagExtra;
          if (dirichletThresh >= 0.0 && TST::real(vals[diagIndex]) <= dirichletThresh) {
            //              SC A_rowsum = ZERO, F_rowsum = ZERO;
            //              printf("WARNING: row %d diag(Afiltered) = %8.2e diag(A)=%8.2e\n",row,vals[diagIndex],diagA);
            //              for(LO l = 0; l < (LO)nnz; l++)
            //                F_rowsum += vals[l];
            //              printf("       : A rowsum = %8.2e F rowsum = %8.2e\n",A_rowsum,F_rowsum);
            vals[diagIndex] = TST::one();
          }
        }
      }
      inds.resize(numInds);
      vals.resize(numInds);

      // Because we used a column map in the construction of the matrix
      // we can just use insertLocalValues here instead of insertGlobalValues
      filteredA.insertLocalValues(row, inds, vals);
    }

    // Reset filtering array
    for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
      for (size_t k = 0; k < blkSize; k++)
        filter[indsG[j] * blkSize + k] = 0;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildNewUsingRootStencil(const Matrix& A, const GraphBase& G, double dirichletThresh, Level& currentLevel, Matrix& filteredA, bool use_spread_lumping, double DdomAllowGrowthRate, double DdomCap) const {
  using TST = typename Teuchos::ScalarTraits<SC>;
  using Teuchos::arcp_const_cast;
  SC ZERO    = Teuchos::ScalarTraits<SC>::zero();
  SC ONE     = Teuchos::ScalarTraits<SC>::one();
  LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  size_t numNodes = G.GetNodeNumVertices();
  size_t blkSize  = A.GetFixedBlockSize();
  size_t numRows  = A.getMap()->getLocalNumElements();
  ArrayView<const LO> indsA;
  ArrayView<const SC> valsA;
  ArrayRCP<const size_t> rowptr;
  ArrayRCP<const LO> inds;
  ArrayRCP<const SC> vals_const;
  ArrayRCP<SC> vals;

  // We're going to grab the vals array from filteredA and then blitz it with NAN as a placeholder for "entries that have
  // not yey been touched."  If I see an entry in the primary loop that has a zero, then I assume it has been nuked by
  // it's symmetric pair, so I add it to the diagonal.  If it has a NAN, process as normal.
  RCP<CrsMatrix> filteredAcrs = dynamic_cast<const CrsMatrixWrap*>(&filteredA)->getCrsMatrix();
  filteredAcrs->getAllValues(rowptr, inds, vals_const);
  vals = arcp_const_cast<SC>(vals_const);
  Array<bool> vals_dropped_indicator(vals.size(), false);

  // In the badAggNeighbors loop, if the entry has any number besides NAN, I add it to the diagExtra and then zero the guy.
  RCP<Aggregates> aggregates      = Get<RCP<Aggregates> >(currentLevel, "Aggregates");
  RCP<AmalgamationInfo> amalgInfo = Get<RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");
  LO numAggs                      = aggregates->GetNumAggregates();

  // Check map nesting
  RCP<const Map> rowMap = A.getRowMap();
  RCP<const Map> colMap = A.getColMap();
  bool goodMap          = MueLu::Utilities<SC, LO, GO, NO>::MapsAreNested(*rowMap, *colMap);
  TEUCHOS_TEST_FOR_EXCEPTION(!goodMap, Exceptions::RuntimeError, "FilteredAFactory: Maps are not nested");

  // Since we're going to symmetrize this
  Array<LO> diagIndex(numRows, INVALID);
  Array<SC> diagExtra(numRows, ZERO);

  // Lists of nodes in each aggregate
  struct {
    // GH: For now, copy everything to host until we properly set this factory to run device code
    // Instead, we'll copy data into HostMirrors and run the algorithms on host, saving optimization for later.
    typename Aggregates::LO_view ptr, nodes, unaggregated;
    typename Aggregates::LO_view::HostMirror ptr_h, nodes_h, unaggregated_h;
  } nodesInAgg;
  aggregates->ComputeNodesInAggregate(nodesInAgg.ptr, nodesInAgg.nodes, nodesInAgg.unaggregated);
  nodesInAgg.ptr_h          = Kokkos::create_mirror_view(nodesInAgg.ptr);
  nodesInAgg.nodes_h        = Kokkos::create_mirror_view(nodesInAgg.nodes);
  nodesInAgg.unaggregated_h = Kokkos::create_mirror_view(nodesInAgg.unaggregated);
  Kokkos::deep_copy(nodesInAgg.ptr_h, nodesInAgg.ptr);
  Kokkos::deep_copy(nodesInAgg.nodes_h, nodesInAgg.nodes);
  Kokkos::deep_copy(nodesInAgg.unaggregated_h, nodesInAgg.unaggregated);
  Teuchos::ArrayRCP<const LO> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);  // GH: this is needed on device, grab the pointer after we call ComputeNodesInAggregate

  LO graphNumCols = G.GetImportMap()->getLocalNumElements();
  Array<bool> filter(graphNumCols, false);

  // Loop over the unaggregated nodes. Blitz those rows. We don't want to smooth singletons.
  for (LO i = 0; i < (LO)nodesInAgg.unaggregated_h.extent(0); i++) {
    for (LO m = 0; m < (LO)blkSize; m++) {
      LO row = amalgInfo->ComputeLocalDOF(nodesInAgg.unaggregated_h(i), m);
      if (row >= (LO)numRows) continue;
      size_t index_start = rowptr[row];
      A.getLocalRowView(row, indsA, valsA);
      for (LO k = 0; k < (LO)indsA.size(); k++) {
        if (row == indsA[k]) {
          vals[index_start + k] = ONE;
          diagIndex[row]        = k;
        } else
          vals[index_start + k] = ZERO;
      }
    }
  }  // end nodesInAgg.unaggregated.extent(0);

  std::vector<LO> badCount(numAggs, 0);

  // Find the biggest aggregate size in *nodes*
  LO maxAggSize = 0;
  for (LO i = 0; i < numAggs; i++)
    maxAggSize = std::max(maxAggSize, nodesInAgg.ptr_h(i + 1) - nodesInAgg.ptr_h(i));

  // Loop over all the aggregates
  std::vector<LO> goodAggNeighbors(G.getLocalMaxNumRowEntries());
  std::vector<LO> badAggNeighbors(std::min(G.getLocalMaxNumRowEntries() * maxAggSize, numNodes));

  size_t numNewDrops   = 0;
  size_t numOldDrops   = 0;
  size_t numFixedDiags = 0;
  size_t numSymDrops   = 0;

  for (LO i = 0; i < numAggs; i++) {
    LO numNodesInAggregate = nodesInAgg.ptr_h(i + 1) - nodesInAgg.ptr_h(i);
    if (numNodesInAggregate == 0) continue;

    // Find the root *node*
    LO root_node = INVALID;
    for (LO k = nodesInAgg.ptr_h(i); k < nodesInAgg.ptr_h(i + 1); k++) {
      if (aggregates->IsRoot(nodesInAgg.nodes_h(k))) {
        root_node = nodesInAgg.nodes_h(k);
        break;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(root_node == INVALID,
                               Exceptions::RuntimeError, "MueLu::FilteredAFactory::BuildNewUsingRootStencil: Cannot find root node");

    // Find the list of "good" node neighbors (aka nodes which border the root node in the Graph G)
    ArrayView<const LO> goodNodeNeighbors = G.getNeighborVertices(root_node);

    // Now find the list of "good" aggregate neighbors (aka the aggregates neighbor the root node in the Graph G)
    goodAggNeighbors.resize(0);
    for (LO k = 0; k < (LO)goodNodeNeighbors.size(); k++) {
      goodAggNeighbors.push_back(vertex2AggId[goodNodeNeighbors[k]]);
    }
    sort_and_unique(goodAggNeighbors);

    // Now we get the list of "bad" aggregate neighbors (aka aggregates which border the
    // root node in the original matrix A, which are not goodNodeNeighbors).  Since we
    // don't have an amalgamated version of the original matrix, we use the matrix directly
    badAggNeighbors.resize(0);
    for (LO j = 0; j < (LO)blkSize; j++) {
      LO row = amalgInfo->ComputeLocalDOF(root_node, j);
      if (row >= (LO)numRows) continue;
      A.getLocalRowView(row, indsA, valsA);
      for (LO k = 0; k < (LO)indsA.size(); k++) {
        if ((indsA[k] < (LO)numRows) && (TST::magnitude(valsA[k]) != TST::magnitude(ZERO))) {
          LO node = amalgInfo->ComputeLocalNode(indsA[k]);
          LO agg  = vertex2AggId[node];
          if (!std::binary_search(goodAggNeighbors.begin(), goodAggNeighbors.end(), agg))
            badAggNeighbors.push_back(agg);
        }
      }
    }
    sort_and_unique(badAggNeighbors);

    // Go through the filtered graph and count the number of connections to the badAggNeighbors
    // if there are 2 or more of these connections, remove them from the bad list.

    for (LO k = nodesInAgg.ptr_h(i); k < nodesInAgg.ptr_h(i + 1); k++) {
      ArrayView<const LO> nodeNeighbors = G.getNeighborVertices(k);
      for (LO kk = 0; kk < nodeNeighbors.size(); kk++) {
        if ((vertex2AggId[nodeNeighbors[kk]] >= 0) && (vertex2AggId[nodeNeighbors[kk]] < numAggs))
          (badCount[vertex2AggId[nodeNeighbors[kk]]])++;
      }
    }
    std::vector<LO> reallyBadAggNeighbors(std::min(G.getLocalMaxNumRowEntries() * maxAggSize, numNodes));
    reallyBadAggNeighbors.resize(0);
    for (LO k = 0; k < (LO)badAggNeighbors.size(); k++) {
      if (badCount[badAggNeighbors[k]] <= 1) reallyBadAggNeighbors.push_back(badAggNeighbors[k]);
    }
    for (LO k = nodesInAgg.ptr_h(i); k < nodesInAgg.ptr_h(i + 1); k++) {
      ArrayView<const LO> nodeNeighbors = G.getNeighborVertices(k);
      for (LO kk = 0; kk < nodeNeighbors.size(); kk++) {
        if ((vertex2AggId[nodeNeighbors[kk]] >= 0) && (vertex2AggId[nodeNeighbors[kk]] < numAggs))
          badCount[vertex2AggId[nodeNeighbors[kk]]] = 0;
      }
    }

    // For each of the reallyBadAggNeighbors, we go and blitz their connections to dofs in this aggregate.
    // We remove the INVALID marker when we do this so we don't wind up doubling this up later
    for (LO b = 0; b < (LO)reallyBadAggNeighbors.size(); b++) {
      LO bad_agg = reallyBadAggNeighbors[b];
      for (LO k = nodesInAgg.ptr_h(bad_agg); k < nodesInAgg.ptr_h(bad_agg + 1); k++) {
        LO bad_node = nodesInAgg.nodes_h(k);
        for (LO j = 0; j < (LO)blkSize; j++) {
          LO bad_row = amalgInfo->ComputeLocalDOF(bad_node, j);
          if (bad_row >= (LO)numRows) continue;
          size_t index_start = rowptr[bad_row];
          A.getLocalRowView(bad_row, indsA, valsA);
          for (LO l = 0; l < (LO)indsA.size(); l++) {
            if (indsA[l] < (LO)numRows && vertex2AggId[amalgInfo->ComputeLocalNode(indsA[l])] == i && vals_dropped_indicator[index_start + l] == false) {
              vals_dropped_indicator[index_start + l] = true;
              vals[index_start + l]                   = ZERO;
              diagExtra[bad_row] += valsA[l];
              numSymDrops++;
            }
          }
        }
      }
    }

    // Now lets fill the rows in this aggregate and figure out the diagonal lumping
    // We loop over each node in the aggregate and then over the neighbors of that node

    for (LO k = nodesInAgg.ptr_h(i); k < nodesInAgg.ptr_h(i + 1); k++) {
      LO row_node = nodesInAgg.nodes_h(k);

      // Set up filtering array
      ArrayView<const LO> indsG = G.getNeighborVertices(row_node);
      for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
        filter[indsG[j]] = true;

      for (LO m = 0; m < (LO)blkSize; m++) {
        LO row = amalgInfo->ComputeLocalDOF(row_node, m);
        if (row >= (LO)numRows) continue;
        size_t index_start = rowptr[row];
        A.getLocalRowView(row, indsA, valsA);

        for (LO l = 0; l < (LO)indsA.size(); l++) {
          int col_node = amalgInfo->ComputeLocalNode(indsA[l]);
          bool is_good = filter[col_node];
          if (indsA[l] == row) {
            diagIndex[row]        = l;
            vals[index_start + l] = valsA[l];
            continue;
          }

          // If we've already dropped this guy (from symmetry above), then continue onward
          if (vals_dropped_indicator[index_start + l] == true) {
            if (is_good)
              numOldDrops++;
            else
              numNewDrops++;
            continue;
          }

          // FIXME: I'm assuming vertex2AggId is only length of the rowmap, so
          // we won'd do secondary dropping on off-processor neighbors
          if (is_good && indsA[l] < (LO)numRows) {
            int agg = vertex2AggId[col_node];
            if (std::binary_search(reallyBadAggNeighbors.begin(), reallyBadAggNeighbors.end(), agg))
              is_good = false;
          }

          if (is_good) {
            vals[index_start + l] = valsA[l];
          } else {
            if (!filter[col_node])
              numOldDrops++;
            else
              numNewDrops++;
            diagExtra[row] += valsA[l];
            vals[index_start + l]                   = ZERO;
            vals_dropped_indicator[index_start + l] = true;
          }
        }  // end for l "indsA.size()" loop

      }  // end m "blkSize" loop

      // Clear filtering array
      for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
        filter[indsG[j]] = false;

    }  // end k loop over number of nodes in this agg
  }    // end i loop over numAggs

  if (!use_spread_lumping) {
    // Now do the diagonal modifications in one, final pass
    for (LO row = 0; row < (LO)numRows; row++) {
      if (diagIndex[row] != INVALID) {
        size_t index_start       = rowptr[row];
        size_t diagIndexInMatrix = index_start + diagIndex[row];
        //        printf("diag_vals pre update =  %8.2e\n", vals[diagIndex] );
        vals[diagIndexInMatrix] += diagExtra[row];
        SC A_rowsum = ZERO, A_absrowsum = ZERO, F_rowsum = ZERO;

        if ((dirichletThresh >= 0.0 && TST::real(vals[diagIndexInMatrix]) <= dirichletThresh) || TST::real(vals[diagIndexInMatrix]) == ZERO) {
          if (MUELU_FILTEREDAFACTORY_LOTS_OF_PRINTING > 0) {
            A.getLocalRowView(row, indsA, valsA);
            //            SC diagA = valsA[diagIndex[row]];
            //            printf("WARNING: row %d (diagIndex=%d) diag(Afiltered) = %8.2e diag(A)=%8.2e numInds = %d\n",row,diagIndex[row],vals[diagIndexInMatrix],diagA,(LO)indsA.size());

            for (LO l = 0; l < (LO)indsA.size(); l++) {
              A_rowsum += valsA[l];
              A_absrowsum += std::abs(valsA[l]);
            }
            for (LO l = 0; l < (LO)indsA.size(); l++)
              F_rowsum += vals[index_start + l];
            //                printf("       : A rowsum = %8.2e |A| rowsum = %8.2e rowsum = %8.2e\n",A_rowsum,A_absrowsum,F_rowsum);
            if (MUELU_FILTEREDAFACTORY_LOTS_OF_PRINTING > 1) {
              //            printf("        Avals =");
              //            for(LO l = 0; l < (LO)indsA.size(); l++)
              //              printf("%d(%8.2e)[%d] ",(LO)indsA[l],valsA[l],(LO)l);
              //            printf("\n");
              //            printf("        Fvals =");
              //            for(LO l = 0; l < (LO)indsA.size(); l++)
              //              if(vals[index_start+l] != ZERO)
              //                printf("%d(%8.2e)[%d] ",(LO)indsA[l],vals[index_start+l],(LO)l);
            }
          }
          // Don't know what to do, so blitz the row and dump a one on the diagonal
          for (size_t l = rowptr[row]; l < rowptr[row + 1]; l++) {
            vals[l] = ZERO;
          }
          vals[diagIndexInMatrix] = TST::one();
          numFixedDiags++;
        }
      } else {
        GetOStream(Runtime0) << "WARNING: Row " << row << " has no diagonal " << std::endl;
      }
    } /*end row "numRows" loop"*/
  }

  // Copy all the goop out
  for (LO row = 0; row < (LO)numRows; row++) {
    filteredA.replaceLocalValues(row, inds(rowptr[row], rowptr[row + 1] - rowptr[row]), vals(rowptr[row], rowptr[row + 1] - rowptr[row]));
  }
  if (use_spread_lumping) ExperimentalLumping(A, filteredA, DdomAllowGrowthRate, DdomCap);

  size_t g_newDrops = 0, g_oldDrops = 0, g_fixedDiags = 0;

  MueLu_sumAll(A.getRowMap()->getComm(), numNewDrops, g_newDrops);
  MueLu_sumAll(A.getRowMap()->getComm(), numOldDrops, g_oldDrops);
  MueLu_sumAll(A.getRowMap()->getComm(), numFixedDiags, g_fixedDiags);
  GetOStream(Runtime0) << "Filtering out " << g_newDrops << " edges, in addition to the " << g_oldDrops << " edges dropped earlier" << std::endl;
  GetOStream(Runtime0) << "Fixing " << g_fixedDiags << " zero diagonal values" << std::endl;
}

// fancy lumping trying to not just move everything to the diagonal but to also consider moving
// some lumping to the kept off-diagonals. We basically aim to not increase the diagonal
// dominance in a row. In particular, the goal is that row i satisfies
//
//        lumpedDiagDomMeasure_i   <= rho2
// or
//        lumpedDiagDomMeasure <= rho*unlumpedDiagDomMeasure
//
// NOTE: THIS CODE assumes direct access to a row. See comments above concerning
//       ASSUME_DIRECT_ACCESS_TO_ROW
//
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExperimentalLumping(const Matrix& A, Matrix& filteredA, double irho, double irho2) const {
  using TST = typename Teuchos::ScalarTraits<SC>;
  SC zero   = TST::zero();
  SC one    = TST::one();

  ArrayView<const LO> inds;
  ArrayView<const SC> vals;
  ArrayView<const LO> finds;
  ArrayView<SC> fvals;

  SC PosOffSum, NegOffSum, PosOffDropSum, NegOffDropSum;
  SC diag, gamma, alpha;
  LO NumPosKept, NumNegKept;

  SC noLumpDdom;
  SC numer, denom;
  SC PosFilteredSum, NegFilteredSum;
  SC Target;

  SC rho  = as<Scalar>(irho);
  SC rho2 = as<Scalar>(irho2);

  for (LO row = 0; row < (LO)A.getRowMap()->getLocalNumElements(); row++) {
    noLumpDdom = as<Scalar>(10000.0);  // only used if diagonal is zero
                                       // the whole idea sort of breaks down
                                       // when the diagonal is zero. In particular,
                                       // the old diag dominance ratio is infinity
                                       // ... so what do we want for the new ddom
                                       // ratio. Do we want to allow the diagonal
                                       // to go negative, just to have a better ddom
                                       // ratio? This current choice essentially
                                       // changes 'Target' to a large number
                                       // meaning that we will allow the new
                                       // ddom number to be fairly large (because
                                       // the old one was infinity)

    ArrayView<const SC> tvals;
    A.getLocalRowView(row, inds, vals);
    size_t nnz = inds.size();
    if (nnz == 0) continue;
    filteredA.getLocalRowView(row, finds, tvals);  // assume 2 getLocalRowView()s
                                                   //  have things in same order
    fvals = ArrayView<SC>(const_cast<SC*>(tvals.getRawPtr()), nnz);

    LO diagIndex = -1, fdiagIndex = -1;

    PosOffSum     = zero;
    NegOffSum     = zero;
    PosOffDropSum = zero;
    NegOffDropSum = zero;
    diag          = zero;
    NumPosKept    = 0;
    NumNegKept    = 0;

    // first record diagonal, offdiagonal sums and off diag dropped sums
    for (size_t j = 0; j < nnz; j++) {
      if (inds[j] == row) {
        diagIndex = j;
        diag      = vals[j];
      } else {  // offdiagonal
        if (TST::real(vals[j]) > TST::real(zero))
          PosOffSum += vals[j];
        else
          NegOffSum += vals[j];
      }
    }
    PosOffDropSum = PosOffSum;
    NegOffDropSum = NegOffSum;
    NumPosKept    = 0;
    NumNegKept    = 0;
    LO j          = 0;
    for (size_t jj = 0; jj < (size_t)finds.size(); jj++) {
      while (inds[j] != finds[jj]) j++;  // assumes that finds is in the same order as
                                         // inds ... but perhaps has some entries missing
      if (finds[jj] == row)
        fdiagIndex = jj;
      else {
        if (TST::real(vals[j]) > TST::real(zero)) {
          PosOffDropSum -= fvals[jj];
          if (TST::real(fvals[jj]) != TST::real(zero)) NumPosKept++;
        } else {
          NegOffDropSum -= fvals[jj];
          if (TST::real(fvals[jj]) != TST::real(zero)) NumNegKept++;
        }
      }
    }

    // measure of diagonal dominance if no lumping is done.
    if (TST::magnitude(diag) != TST::magnitude(zero))
      noLumpDdom = (PosOffSum - NegOffSum) / diag;

    // Target is an acceptable diagonal dominance ratio
    // which should really be larger than 1

    Target = rho * noLumpDdom;
    if (TST::magnitude(Target) <= TST::magnitude(rho)) Target = rho2;

    PosFilteredSum = PosOffSum - PosOffDropSum;
    NegFilteredSum = NegOffSum - NegOffDropSum;
    // Note: PosNotFilterdSum is not equal to the sum of the
    // positive entries after lumping. It just reflects the
    // pos offdiag sum of the filtered matrix before lumping
    // and does not account for negative dropped terms lumped
    // to the positive kept terms.

    // dropped positive offdiags always go to the diagonal as these
    // always improve diagonal dominance.

    diag += PosOffDropSum;

    // now lets work on lumping dropped negative offdiags
    gamma = -NegOffDropSum - PosFilteredSum;

    if (TST::real(gamma) < TST::real(zero)) {
      // the total amount of negative dropping is less than PosFilteredSum,
      // so we can distribute this dropping to pos offdiags. After lumping
      // the sum of the pos offdiags is just -gamma so we just assign pos
      // offdiags proportional to vals[j]/PosFilteredSum
      // Note: in this case the diagonal is not changed as all lumping
      //       occurs to the pos offdiags

      if (fdiagIndex != -1) fvals[fdiagIndex] = diag;
      j = 0;
      for (LO jj = 0; jj < (LO)finds.size(); jj++) {
        while (inds[j] != finds[jj]) j++;  // assumes that finds is in the same order as
                                           // inds ... but perhaps has some entries missing
        if ((j != diagIndex) && (TST::real(vals[j]) > TST::real(zero)) && (TST::magnitude(fvals[jj]) != TST::magnitude(zero)))
          fvals[jj] = -gamma * (vals[j] / PosFilteredSum);
      }
    } else {
      // So there are more negative values that need lumping than kept
      // positive offdiags. Meaning there is enough negative lumping to
      // completely clear out all pos offdiags. If we lump all negs
      // to pos off diags, we'd actually change them to negative. We
      // only do this if we are desperate. Otherwise, we'll clear out
      // all the positive kept offdiags and try to lump the rest
      // somewhere else.  We defer the clearing of pos off diags
      // to see first if we are going to be desperate.

      bool flipPosOffDiagsToNeg = false;

      // Even if we lumped by zeroing positive offdiags, we are still
      // going to have more lumping to distribute to either
      //     1) the diagonal
      //     2) the kept negative offdiags
      //     3) the kept positive offdiags (desperate)

      // Let's first considering lumping the remaining neg offdiag stuff
      // to the diagonal ... if this does not increase the diagonal
      // dominance ratio too much (given  by rho).

      if ((TST::real(diag) > TST::real(gamma)) &&
          (TST::real((-NegFilteredSum) / (diag - gamma)) <= TST::real(Target))) {
        // 1st if term above insures that resulting diagonal (=diag-gamma)
        // is positive. . The left side of 2nd term is the diagonal dominance
        // if we lump the remaining stuff (gamma) to the diagonal. Recall,
        // that now there are no positive off-diags so the sum(abs(offdiags))
        // is just the negative of NegFilteredSum

        if (fdiagIndex != -1) fvals[fdiagIndex] = diag - gamma;
      } else if (NumNegKept > 0) {
        // need to do some lumping to neg offdiags to avoid a large
        // increase in diagonal dominance. We first compute alpha
        // which measures how much gamma should go to the
        // negative offdiags. The rest will go to the diagonal

        numer = -NegFilteredSum - Target * (diag - gamma);
        denom = gamma * (Target - TST::one());

        // make sure that alpha is between 0 and 1 ... and that it doesn't
        // result in a sign flip
        //    Note: when alpha is set to 1, then the diagonal is not modified
        //          and the negative offdiags just get shifted from those
        //          removed and those kept, meaning that the digaonal dominance
        //          should be the same as before
        //
        //          can alpha be negative? It looks like denom should always
        //          be positive. The 'if' statement above
        //          Normally, diag-gamma should also be positive (but if it
        //          is negative then numer is guaranteed to be positve).
        //          look at the 'if' above,
        //                if (( TST::real(diag) > TST::real(gamma)) &&
        //                ( TST::real((-NegFilteredSum)/(diag - gamma)) <= TST::real(Target))) {
        //
        //          Should guarantee that numer is positive. This is obvious when
        //          the second condition is false. When it is the first condition that
        //          is false, it follows that the two indiviudal terms in the numer
        //          formula must be positive.

        if (TST::magnitude(denom) < TST::magnitude(numer))
          alpha = TST::one();
        else
          alpha = numer / denom;
        if (TST::real(alpha) < TST::real(zero)) alpha = zero;
        if (TST::real(diag) < TST::real((one - alpha) * gamma)) alpha = TST::one();

        // first change the diagonal

        if (fdiagIndex != -1) fvals[fdiagIndex] = diag - (one - alpha) * gamma;

        // after lumping the sum of neg offdiags will be NegFilteredSum
        // + alpha*gamma. That is the remaining negative entries altered
        // by the percent (=alpha) of stuff (=gamma) that needs to be
        // lumped after taking into account lumping to pos offdiags

        // Do this by assigning a fraction of NegFilteredSum+alpha*gamma
        //  proportional to vals[j]/NegFilteredSum

        SC temp = (NegFilteredSum + alpha * gamma) / NegFilteredSum;
        j       = 0;
        for (LO jj = 0; jj < (LO)finds.size(); jj++) {
          while (inds[j] != finds[jj]) j++;  // assumes that finds is in the same order as
                                             // inds ... but perhaps has some entries missing
          if ((jj != fdiagIndex) && (TST::magnitude(fvals[jj]) != TST::magnitude(zero)) &&
              (TST::real(vals[j]) < TST::real(zero)))
            fvals[jj] = temp * vals[j];
        }
      } else {  // desperate case
        // So we don't have any kept negative offdiags  ...

        if (NumPosKept > 0) {
          // luckily we can push this stuff to the pos offdiags
          // which now makes them negative
          flipPosOffDiagsToNeg = true;

          j = 0;
          for (LO jj = 0; jj < (LO)finds.size(); jj++) {
            while (inds[j] != finds[jj]) j++;  // assumes that finds is in the same order as
                                               // inds ... but perhaps has some entries missing
            if ((j != diagIndex) && (TST::magnitude(fvals[jj]) != TST::magnitude(zero)) &&
                (TST::real(vals[j]) > TST::real(zero)))
              fvals[jj] = -gamma / ((SC)NumPosKept);
          }
        }
        // else abandon rowsum preservation and do nothing
      }
      if (!flipPosOffDiagsToNeg) {  // not desperate so we now zero out
                                    // all pos terms including some
                                    // not originally filtered
                                    // but zeroed due to lumping
        j = 0;
        for (LO jj = 0; jj < (LO)finds.size(); jj++) {
          while (inds[j] != finds[jj]) j++;  // assumes that finds is in the same order as
                                             // inds ... but perhaps has some entries missing
          if ((jj != fdiagIndex) && (TST::real(vals[j]) > TST::real(zero))) fvals[jj] = zero;
        }
      }
    }  // positive gamma else

  }  // loop over all rows
}

}  // namespace MueLu

#endif  // MUELU_FILTEREDAFACTORY_DEF_HPP
