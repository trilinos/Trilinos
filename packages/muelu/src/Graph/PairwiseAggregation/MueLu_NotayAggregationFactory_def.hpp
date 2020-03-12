
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
#ifndef MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_spgemm.hpp"

#include "MueLu_NotayAggregationFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: pairwise: size");
    SET_VALID_ENTRY("aggregation: compute aggregate qualities");
    SET_VALID_ENTRY("aggregation: Dirichlet threshold");
#undef SET_VALID_ENTRY



    /*
    validParamList->getEntry("aggregation: ordering").setValidator(
	  rcp(new validatorType(Teuchos::tuple<std::string>("natural", "graph", "random", "cuthill-mckee"), "aggregation: ordering")));
    */

    // general variables needed in AggregationFactory
    validParamList->set< RCP<const FactoryBase> >("A",           null, "Generating factory of the matrix");
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
    validParamList->set< RCP<const FactoryBase> >("AggregateQualities", null, "Generating factory for variable \'AggregateQualities\'");


    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");
    Input(currentLevel, "A");
    Input(currentLevel, "DofsPerNode");

    const ParameterList& pL = GetParameterList();

    if (pL.get<bool>("aggregation: compute aggregate qualities")) {
        Input(currentLevel, "AggregateQualities");
    }


  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    const ParameterList& pL = GetParameterList();

    // Parameters
    int aggregateTargetSize = 2;
    if (pL.isParameter("aggregation: pairwise: size"))
      aggregateTargetSize = pL.get<int>("aggregation: pairwise: size");
    TEUCHOS_TEST_FOR_EXCEPTION(aggregateTargetSize != 2 && aggregateTargetSize != 4 && aggregateTargetSize != 8,
			       Exceptions::RuntimeError,
			       "NotayAggregationFactory::Build(): aggregateTargetSize needs to be a power of two");


    RCP<const GraphBase> graph = Get< RCP<GraphBase> >(currentLevel, "Graph");
    RCP<const Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    // Setup aggregates & aggStat objects
    RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("PW");

    // Create intermediate prolongator
    typename Matrix::local_matrix_type intermediateP;
    typename Matrix::local_matrix_type coarseLocalA;

    const LO numRows = graph->GetNodeNumVertices();

    // construct aggStat information
    std::vector<unsigned> aggStat(numRows, READY);


    const int DofsPerNode =  Get<int>(currentLevel,"DofsPerNode");
    TEUCHOS_TEST_FOR_EXCEPTION(DofsPerNode != 1, Exceptions::RuntimeError,
			       "Pairwise only supports one dof per node");

    // This follows the paper:
    // Notay, "Aggregation-based algebraic multigrid for convection-diffusion equations", SISC 34(3), pp. A2288-2316.


    // FIXME: Do the ordering


    // Get the party stated
    LO numNonAggregatedNodes = numRows, numDirichletNodes = 0;
    BuildInitialAggregates(pL, A, *aggregates, aggStat, numNonAggregatedNodes, numDirichletNodes);
    TEUCHOS_TEST_FOR_EXCEPTION(0 < numNonAggregatedNodes, Exceptions::RuntimeError,
                               "Initial pairwise aggregation failed to aggregate all nodes");

    // Compute the intermediate prolongator
    BuildIntermediateProlongator(numRows, numDirichletNodes, aggregates, intermediateP);

    // Compute the coarse local matrix
    BuildCoarseLocalMatrix(A->getLocalMatrix(), intermediateP, coarseLocalA);


    // DO stuff
    Set(currentLevel, "Aggregates", aggregates);
    GetOStream(Statistics0) << aggregates->description() << std::endl;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildInitialAggregates(const Teuchos::ParameterList& params,
                         const RCP<const Matrix>& A,
                         Aggregates& aggregates,
                         std::vector<unsigned>& aggStat,
                         LO& numNonAggregatedNodes,
                         LO& numDirichletNodes) const {

    Monitor m(*this, "BuildInitialAggregates");
    using STS = Teuchos::ScalarTraits<Scalar>;
    using MT  = typename STS::magnitudeType;
    using RealValuedVector = Xpetra::Vector<MT,LocalOrdinal,GlobalOrdinal,Node>;

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    const SC SC_ZERO    = Teuchos::ScalarTraits<SC>::zero();
    const MT MT_ZERO    = Teuchos::ScalarTraits<MT>::zero();
    const MT MT_ONE     = Teuchos::ScalarTraits<MT>::one();
    const MT MT_TWO     = MT_ONE + MT_ONE;
    const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
    const LO LO_ZERO    = Teuchos::OrdinalTraits<LO>::zero();

    const MT kappa       = STS::magnitude(static_cast<SC>(params.get<double>("aggregation: Dirichlet threshold")));
    TEUCHOS_TEST_FOR_EXCEPTION(kappa <= 2, Exceptions::RuntimeError,
			       "Pairwise requires kappa > 2 otherwise all rows are considered as Dirichlet rows.");
    const MT kappa_init  = kappa / (kappa - MT_TWO);

    const LO numRows = aggStat.size();
    const int myRank  = A->getMap()->getComm()->getRank();

    // NOTE: Assumes 1 dof per node.  This constraint is enforced in Build(), and so we're not doing again here.
    // This should probably be fixed at some point.

    // Extract diagonal, rowsums, etc
    RCP<Vector> ghostedDiag = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
    RCP<Vector> ghostedRowSum = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDeletedRowsum(*A);
    RCP<RealValuedVector> ghostedAbsRowSum = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedAbsDeletedRowsum(*A);
    const ArrayRCP<const SC> D     = ghostedDiag->getData(0);
    const ArrayRCP<const SC> S     = ghostedRowSum->getData(0);
    const ArrayRCP<const MT> AbsRs = ghostedAbsRowSum->getData(0);

    // Aggregates stuff
    ArrayRCP<LO>  vertex2AggId_rcp = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO>  procWinner_rcp   = aggregates.GetProcWinner()  ->getDataNonConst(0);
    ArrayView<LO> vertex2AggId     = vertex2AggId_rcp();
    ArrayView<LO> procWinner       = procWinner_rcp();

    // Algorithm 4.2

    // 0,1 : Initialize: Flag boundary conditions
    // Modification: We assume symmetry here aij = aji

    //    printf("numRows = %d, A->getRowMap()->getNodeNumElements() = %d\n",(int)numRows,(int) A->getRowMap()->getNodeNumElements());

    for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getNodeNumElements()); ++row) {
      MT aii = STS::magnitude(D[row]);
      MT rowsum = AbsRs[row];

      if(aii >= kappa_init * rowsum) {
        *out << "Flagging index " << row << "as dirichlet "
          "aii >= kappa*rowsum = " << aii << " >= " << kappa_init << " " << rowsum << std::endl;
	aggStat[row] = IGNORED;
	--numNonAggregatedNodes;
        ++numDirichletNodes;
      }
    }

    // FIXME: Add ordering here or pass it in from Build()

    // 2 : Iteration
    LO aggIndex = LO_ZERO;
    for(LO current_idx = 0; current_idx < numRows; current_idx++) {
      // If we're aggregated already, skip this guy
      if(aggStat[current_idx] != READY)
        continue;

      MT best_mu = MT_ZERO;
      LO best_idx = LO_INVALID;

      size_t nnz = A->getNumEntriesInLocalRow(current_idx);
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      A->getLocalRowView(current_idx, indices, vals);

      MT aii = STS::real(D[current_idx]);
      MT si  = STS::real(S[current_idx]);
      for (LO colidx = 0; colidx < static_cast<LO>(nnz); colidx++) {
        // Skip aggregated neighbors, off-rank neighbors, hard zeros and self
        LO col = indices[colidx];
        SC val = vals[colidx];
        if(current_idx == col || aggStat[col] != READY || col > numRows || val == SC_ZERO)
          continue;


	MT aij = STS::real(val);
	MT ajj = STS::real(D[col]);
	MT sj  = STS::real(S[col]);
	if(aii - si + ajj - sj >= MT_ZERO) {
          // Modification: We assume symmetry here aij = aji
	  MT mu_top    = MT_TWO / ( MT_ONE / aii + MT_ONE / ajj);
	  MT mu_bottom =  - aij + MT_ONE / ( MT_ONE / (aii - si) + MT_ONE / (ajj - sj) );
	  MT mu = mu_top / mu_bottom;
	  if (best_idx == LO_INVALID ||  mu < best_mu) {
	    best_mu  = mu;
	    best_idx = col;
	  }
          *out << "[" << current_idx << "] Column SUCCESS " << col << ": "
              << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
              << " = " << aii - si + ajj - sj << ", mu = " << mu << std::endl;
	}
        else {
          *out << "[" << current_idx << "] Column FAILED " << col << ": "
              << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
              << " = " << aii - si + ajj - sj << std::endl;
        }
      }// end column for loop

      if(best_idx == LO_INVALID) {
        *out << "No node buddy found for index " << current_idx
            << " [agg " << aggIndex << "]\n" << std::endl;
        // We found no potential node-buddy, so let's just make this a singleton
        // NOTE: The behavior of what to do if you have no un-aggregated neighbors is not specified in
        // the paper

        aggStat[current_idx] = ONEPT;
        vertex2AggId[current_idx] = aggIndex;
        procWinner[current_idx]   = myRank;
	numNonAggregatedNodes--;
        aggIndex++;

      } else {
        // We have a buddy, so aggregate, either as a singleton or as a pair, depending on mu
        if(best_mu <= kappa) {
          *out << "Node buddies (" << current_idx << "," << best_idx << ") [agg " << aggIndex << "]" << std::endl;

          aggStat[current_idx] = AGGREGATED;
          aggStat[best_idx]    = AGGREGATED;
          vertex2AggId[current_idx] = aggIndex;
          vertex2AggId[best_idx]    = aggIndex;
          procWinner[current_idx]   = myRank;
          procWinner[best_idx]      = myRank;
          numNonAggregatedNodes-=2;
          aggIndex++;

        } else {
          *out << "No buddy found for index " << current_idx << ","
            " but aggregating as singleton [agg " << aggIndex << "]" << std::endl;

          aggStat[current_idx] = ONEPT;
          vertex2AggId[current_idx] = aggIndex;
          procWinner[current_idx]   = myRank;
          numNonAggregatedNodes--;
          aggIndex++;
        } // best_mu
      } // best_idx
    }// end Algorithm 4.2

    *out << "vertex2aggid :";
    for(int i = 0; i < static_cast<int>(vertex2AggId.size()); ++i) {
      *out << i << "(" << vertex2AggId[i] << ")";
    }
    *out << std::endl;

    // update aggregate object
    aggregates.SetNumAggregates(aggIndex);
    aggregates.AggregatesCrossProcessors(false);
    aggregates.ComputeAggregateSizes(true/*forceRecompute*/);
  } // BuildInitialAggregates

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildIntermediateProlongator(const LocalOrdinal numRows,
                               const LocalOrdinal numDirichletNodes,
                               const RCP<Aggregates>& aggregates,
                               typename Matrix::local_matrix_type& intermediateP) const {
    Monitor m(*this, "BuildIntermediateProlongator");

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    using local_matrix_type = typename Matrix::local_matrix_type;
    using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using values_type       = typename local_matrix_type::values_type;
    using size_type         = typename local_graph_type::size_type;
    using col_index_type    = typename local_graph_type::data_type;
    using array_layout      = typename local_graph_type::array_layout;
    using device_type       = typename local_graph_type::device_type;
    using memory_traits     = typename local_graph_type::memory_traits;
    using row_pointer_type  = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
    using col_indices_type  = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;

    const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    const int intermediatePnnz = numRows - numDirichletNodes;
    row_pointer_type rowPtr("intermediateP row pointer", numRows + 1);
    col_indices_type colInd("intermediateP column indices", intermediatePnnz);
    values_type      values("intermediateP values", intermediatePnnz);
    typename row_pointer_type::HostMirror rowPtr_h = Kokkos::create_mirror_view(rowPtr);
    typename col_indices_type::HostMirror colInd_h = Kokkos::create_mirror_view(colInd);

    ArrayView<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0).view(0, numRows);

    rowPtr_h(0) = 0;
    for(int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      // Skip Dirichlet nodes
      if(vertex2AggId[rowIdx] == LO_INVALID) {
        rowPtr_h(rowIdx + 1) = rowPtr_h(rowIdx);
      } else {
        rowPtr_h(rowIdx + 1) = rowPtr_h(rowIdx) + 1;
        colInd_h(rowPtr_h(rowIdx)) = vertex2AggId[rowIdx];
      }
    }

    Kokkos::deep_copy(rowPtr, rowPtr_h);
    Kokkos::deep_copy(colInd, colInd_h);
    Kokkos::deep_copy(values, Kokkos::ArithTraits<typename values_type::value_type>::one());

    // *out << "numDirichletNodes=" << numDirichletNodes
    //      << ", numAggregates=" << aggregates->GetNumAggregates()
    //      << ", numNnz=" << numRows - numDirichletNodes << std::endl;

    intermediateP = local_matrix_type("intermediateP",
                                      numRows, aggregates->GetNumAggregates(), intermediatePnnz,
                                      values, rowPtr, colInd);

    *out << printLocalMatrix("intermediateP", intermediateP);
  } // BuildIntermediateProlongator

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildCoarseLocalMatrix(const typename Matrix::local_matrix_type& A,
                         const typename Matrix::local_matrix_type& intermediateP,
                         typename Matrix::local_matrix_type& coarseA) const {
    Monitor m(*this, "BuildCoarseLocalMatrix");

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    using local_matrix_type = typename Matrix::local_matrix_type;
    using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using values_type       = typename local_matrix_type::values_type;
    using size_type         = typename local_graph_type::size_type;
    using col_index_type    = typename local_graph_type::data_type;
    using array_layout      = typename local_graph_type::array_layout;
    using device_type       = typename local_graph_type::device_type;
    using memory_traits     = typename local_graph_type::memory_traits;
    using row_pointer_type  = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
    using col_indices_type  = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;

    const int numRows = A.numRows();

    // Extract on rank part of A
    // Simply check that the column index is less than the number of local rows
    // otherwise remove it.
    local_matrix_type onrankA;

    row_pointer_type rowPtr("onrankA row pointer", numRows + 1);
    typename row_pointer_type::HostMirror rowPtr_h = Kokkos::create_mirror_view(rowPtr);
    typename local_matrix_type::row_map_type::HostMirror origRowPtr_h
      = Kokkos::create_mirror_view(A.graph.row_map);
    typename local_graph_type::entries_type::HostMirror origColind_h
      = Kokkos::create_mirror_view(A.graph.entries);
    typename values_type::HostMirror origValues_h
      = Kokkos::create_mirror_view(A.values);
    Kokkos::deep_copy(origRowPtr_h, A.graph.row_map);
    Kokkos::deep_copy(origColind_h, A.graph.entries);
    Kokkos::deep_copy(origValues_h, A.values);

    // Compute the number of nnz entries per row
    rowPtr_h(0) = 0;
    for(int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      for(size_type entryIdx = origRowPtr_h(rowIdx); entryIdx < origRowPtr_h(rowIdx + 1); ++entryIdx) {
        if(origColind_h(entryIdx) < numRows) {rowPtr_h(rowIdx + 1) += 1;}
      }
      rowPtr_h(rowIdx + 1) = rowPtr_h(rowIdx + 1) + rowPtr_h(rowIdx);
    }
    Kokkos::deep_copy(rowPtr, rowPtr_h);

    const LO nnzOnrankA = rowPtr_h(numRows);

    // Now use nnz per row to allocate matrix views and store column indices and values
    col_indices_type colInd("onrankA column indices", rowPtr_h(numRows));
    values_type      values("onrankA values", rowPtr_h(numRows));
    typename col_indices_type::HostMirror colInd_h = Kokkos::create_mirror_view(colInd);
    typename values_type::HostMirror      values_h = Kokkos::create_mirror_view(values);
    int entriesInRow;
    for(int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      entriesInRow = 0;
      for(size_type entryIdx = origRowPtr_h(rowIdx); entryIdx < origRowPtr_h(rowIdx + 1); ++entryIdx) {
        if(origColind_h(entryIdx) < numRows) {
          colInd_h(rowPtr_h(rowIdx) + entriesInRow) = origColind_h(entryIdx);
          values_h(rowPtr_h(rowIdx) + entriesInRow) = origValues_h(entryIdx);
          ++entriesInRow;
        }
      }
    }
    Kokkos::deep_copy(colInd, colInd_h);
    Kokkos::deep_copy(values, values_h);

    onrankA = local_matrix_type("on-rank A", numRows, numRows,
                                nnzOnrankA, values, rowPtr, colInd);

    *out << printLocalMatrix("onrankA", onrankA);

    local_matrix_type AP;

    localSpGEMM(onrankA, intermediateP, "AP", AP);

    // Note 03/11/20, lbv: does kh need to destroy and recreate the spgemm handle
    // I am not sure but doing it for safety in case it stashes data from the previous
    // spgemm computation...

    // Compute Ac = Pt * AP
    // Two steps needed:
    //   1. compute Pt
    //   2. perform multiplication

    // Step 1 compute Pt
    // Obviously this requires the same amount of storage as P except for the rowPtr
    row_pointer_type rowPtrPt(Kokkos::ViewAllocateWithoutInitializing("Pt row pointer"),
                              intermediateP.numCols() + 1);
    col_indices_type colIndPt(Kokkos::ViewAllocateWithoutInitializing("Pt column indices"),
                              intermediateP.nnz());
    values_type      valuesPt(Kokkos::ViewAllocateWithoutInitializing("Pt values"),
                              intermediateP.nnz());

    typename row_pointer_type::HostMirror rowPtrPt_h = Kokkos::create_mirror_view(rowPtrPt);
    Kokkos::deep_copy(rowPtrPt_h, 0);
    for(size_type entryIdx = 0; entryIdx < intermediateP.nnz(); ++entryIdx) {
      rowPtrPt_h(intermediateP.graph.entries(entryIdx) + 1) += 1;
    }
    for(LO rowIdx = 0; rowIdx < intermediateP.numCols(); ++rowIdx) {
      rowPtrPt_h(rowIdx + 1) += rowPtrPt_h(rowIdx);
    }
    Kokkos::deep_copy(rowPtrPt, rowPtrPt_h);

    typename row_pointer_type::HostMirror rowPtrP_h = Kokkos::create_mirror_view(intermediateP.graph.row_map);
    Kokkos::deep_copy(rowPtrP_h, intermediateP.graph.row_map);
    typename col_indices_type::HostMirror colIndP_h = Kokkos::create_mirror_view(intermediateP.graph.entries);
    Kokkos::deep_copy(colIndP_h, intermediateP.graph.entries);
    typename values_type::HostMirror valuesP_h  = Kokkos::create_mirror_view(intermediateP.values);
    Kokkos::deep_copy(valuesP_h, intermediateP.values);
    typename col_indices_type::HostMirror colIndPt_h = Kokkos::create_mirror_view(colIndPt);
    typename values_type::HostMirror valuesPt_h = Kokkos::create_mirror_view(valuesPt);
    const col_index_type invalidColumnIndex = KokkosSparse::OrdinalTraits<col_index_type>::invalid();
    Kokkos::deep_copy(colIndPt_h, invalidColumnIndex);

    *out << "allocated Pt's column indices and values" << std::endl;

    col_index_type colIdx = 0;
    for(LO rowIdx = 0; rowIdx < intermediateP.numRows(); ++rowIdx) {
      for(size_type entryIdxP = rowPtrP_h(rowIdx); entryIdxP < rowPtrP_h(rowIdx + 1); ++entryIdxP) {
        colIdx = intermediateP.graph.entries(entryIdxP);
        for(size_type entryIdxPt = rowPtrPt_h(colIdx); entryIdxPt < rowPtrPt_h(colIdx + 1); ++entryIdxPt) {
          if(colIndPt_h(entryIdxPt) == invalidColumnIndex) {
            colIndPt_h(entryIdxPt) = rowIdx;
            valuesPt_h(entryIdxPt) = valuesP_h(entryIdxP);
            break;
          }
        } // Loop over entries in row of Pt
      } // Loop over entries in row of P
    } // Loop over rows of P

    Kokkos::deep_copy(colIndPt, colIndPt_h);
    Kokkos::deep_copy(valuesPt, valuesPt_h);

    local_matrix_type intermediatePt("intermediatePt",
                                     intermediateP.numCols(),
                                     intermediateP.numRows(),
                                     intermediateP.nnz(),
                                     valuesPt, rowPtrPt, colIndPt);

    *out << printLocalMatrix("intermediatePt", intermediatePt);
    *out << "computed Pt" << std::endl;


    // Create views for coarseA matrix
    localSpGEMM(intermediatePt, AP, "coarseA", coarseA);

    *out << printLocalMatrix("coarseA", coarseA);
    *out << "created Ac" << std::endl;

  } // BuildCoarseLocalMatrix

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  localSpGEMM(const typename Matrix::local_matrix_type& A,
              const typename Matrix::local_matrix_type& B,
              const std::string matrixLabel,
              typename Matrix::local_matrix_type& C) const {
    Monitor m(*this, "localSpGEMM");

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    using local_matrix_type = typename Matrix::local_matrix_type;
    using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using values_type       = typename local_matrix_type::values_type;
    using size_type         = typename local_graph_type::size_type;
    using col_index_type    = typename local_graph_type::data_type;
    using array_layout      = typename local_graph_type::array_layout;
    using device_type       = typename local_graph_type::device_type;
    using memory_traits     = typename local_graph_type::memory_traits;
    using row_pointer_type  = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
    using col_indices_type  = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;

    // Options
    int team_work_size = 16;
    std::string myalg("SPGEMM_KK_MEMORY");
    KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);
    KokkosKernels::Experimental::KokkosKernelsHandle<typename row_pointer_type::const_value_type,
                                                     typename col_indices_type::const_value_type,
                                                     typename values_type::const_value_type,
                                                     typename device_type::execution_space,
                                                     typename device_type::memory_space,
                                                     typename device_type::memory_space> kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);

    // Create views for AP matrix
    row_pointer_type rowPtrC(Kokkos::ViewAllocateWithoutInitializing("C row pointer"),
                             A.numRows() + 1);
    col_indices_type colIndC;
    values_type      valuesC;

    // Symbolic multiplication
    KokkosSparse::Experimental::spgemm_symbolic(&kh, A.numRows(),
                                                B.numRows(), B.numCols(),
                                                A.graph.row_map, A.graph.entries, false,
                                                B.graph.row_map, B.graph.entries, false,
                                                rowPtrC);

    // allocate column indices and values of AP
    size_t nnzC = kh.get_spgemm_handle()->get_c_nnz();
    if (nnzC) {
      colIndC = col_indices_type(Kokkos::ViewAllocateWithoutInitializing("C column inds"), nnzC);
      valuesC = values_type(Kokkos::ViewAllocateWithoutInitializing("C values"), nnzC);
    }

    *out << "computed symbolic C" << std::endl;

    // Numeric multiplication
    KokkosSparse::Experimental::spgemm_numeric(&kh, A.numRows(),
                                               B.numRows(), B.numCols(),
                                               A.graph.row_map, A.graph.entries, A.values, false,
                                               B.graph.row_map, B.graph.entries, B.values, false,
                                               rowPtrC, colIndC, valuesC);
    kh.destroy_spgemm_handle();

    *out << "computed numeric C" << std::endl;

    C = local_matrix_type(matrixLabel, A.numRows(), B.numCols(), nnzC, valuesC, rowPtrC, colIndC);

  } // localSpGEMM

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  printLocalMatrix(const std::string& matrixLabel,
                   const typename Matrix::local_matrix_type& A) const {
    // This function is used for debugging purposes only and does not run on Cuda
    // It could implement the construction of a hostA the uses the copy constructor
    // of CrsMatrix with HostMirror views of the original A.
    // At the moment let us simply check the memory space of views and throw if it is not
    // host compatible.

    using size_type      = typename Matrix::local_matrix_type::staticcrsgraph_type::size_type;
    using col_index_type = typename Matrix::local_matrix_type::staticcrsgraph_type::data_type;

    std::ostringstream matrix_description;

    if(std::is_same<typename Matrix::local_matrix_type::values_type::memory_space,
       typename Matrix::local_matrix_type::values_type::HostMirror::memory_space>::value &&
       std::is_same<typename Matrix::local_matrix_type::values_type::data_type,
       typename Matrix::local_matrix_type::values_type::HostMirror::data_type>::value) {

      matrix_description << matrixLabel << ":" << std::endl
                         << "  - numRows=" << A.numRows() << std::endl
                         << "  - numCols=" << A.numCols() << std::endl
                         << "  - nnz=    " << A.nnz() << std::endl
                         << "  - rows:" << std::endl;
      for(col_index_type rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        matrix_description << "      row " << rowIdx
                           << ", entries " << A.graph.row_map(rowIdx)
                           << " to " << A.graph.row_map(rowIdx + 1) << " { ";
        for(size_type entryIdx = A.graph.row_map(rowIdx); entryIdx < A.graph.row_map(rowIdx + 1); ++entryIdx) {
          matrix_description << "(" << A.graph.entries(entryIdx) << ", " << A.values(entryIdx) << ") ";
        }
        matrix_description << "}" << std::endl;
      }
    } else {
      throw std::runtime_error("You can only print matrices that are on host");
    }

    return matrix_description.str();
  }


} //namespace MueLu

#endif /* MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_ */
