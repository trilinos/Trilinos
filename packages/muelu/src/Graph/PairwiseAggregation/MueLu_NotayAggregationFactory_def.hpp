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
#include "KokkosSparse_spmv.hpp"

#include "MueLu_NotayAggregationFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

namespace NotayUtils {
template <class LocalOrdinal>
LocalOrdinal RandomOrdinal(LocalOrdinal min, LocalOrdinal max) {
  return min + as<LocalOrdinal>((max - min + 1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
}

template <class LocalOrdinal>
void RandomReorder(Teuchos::Array<LocalOrdinal>& list) {
  typedef LocalOrdinal LO;
  LO n = Teuchos::as<LO>(list.size());
  for (LO i = 0; i < n - 1; i++)
    std::swap(list[i], list[RandomOrdinal(i, n - 1)]);
}
}  // namespace NotayUtils

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: pairwise: size");
  SET_VALID_ENTRY("aggregation: pairwise: tie threshold");
  SET_VALID_ENTRY("aggregation: compute aggregate qualities");
  SET_VALID_ENTRY("aggregation: Dirichlet threshold");
  SET_VALID_ENTRY("aggregation: ordering");
#undef SET_VALID_ENTRY

  // general variables needed in AggregationFactory
  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix");
  validParamList->set<RCP<const FactoryBase> >("Graph", null, "Generating factory of the graph");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
  validParamList->set<RCP<const FactoryBase> >("AggregateQualities", null, "Generating factory for variable \'AggregateQualities\'");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = GetParameterList();

  Input(currentLevel, "A");
  Input(currentLevel, "Graph");
  Input(currentLevel, "DofsPerNode");
  if (pL.get<bool>("aggregation: compute aggregate qualities")) {
    Input(currentLevel, "AggregateQualities");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);
  using STS = Teuchos::ScalarTraits<Scalar>;
  using MT  = typename STS::magnitudeType;

  const MT MT_TWO = Teuchos::ScalarTraits<MT>::one() + Teuchos::ScalarTraits<MT>::one();

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  const ParameterList& pL = GetParameterList();

  const MT kappa = static_cast<MT>(pL.get<double>("aggregation: Dirichlet threshold"));
  TEUCHOS_TEST_FOR_EXCEPTION(kappa <= MT_TWO,
                             Exceptions::RuntimeError,
                             "Pairwise requires kappa > 2"
                             " otherwise all rows are considered as Dirichlet rows.");

  // Parameters
  int maxNumIter = 3;
  if (pL.isParameter("aggregation: pairwise: size"))
    maxNumIter = pL.get<int>("aggregation: pairwise: size");
  TEUCHOS_TEST_FOR_EXCEPTION(maxNumIter < 1,
                             Exceptions::RuntimeError,
                             "NotayAggregationFactory::Build(): \"aggregation: pairwise: size\""
                             " must be a strictly positive integer");

  RCP<const GraphBase> graph = Get<RCP<GraphBase> >(currentLevel, "Graph");
  RCP<const Matrix> A        = Get<RCP<Matrix> >(currentLevel, "A");

  // Setup aggregates & aggStat objects
  RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
  aggregates->setObjectLabel("PW");

  const LO numRows = graph->GetNodeNumVertices();

  // construct aggStat information
  std::vector<unsigned> aggStat(numRows, READY);

  const int DofsPerNode = Get<int>(currentLevel, "DofsPerNode");
  TEUCHOS_TEST_FOR_EXCEPTION(DofsPerNode != 1, Exceptions::RuntimeError,
                             "Pairwise only supports one dof per node");

  // This follows the paper:
  // Notay, "Aggregation-based algebraic multigrid for convection-diffusion equations",
  // SISC 34(3), pp. A2288-2316.

  // Handle Ordering
  std::string orderingStr = pL.get<std::string>("aggregation: ordering");
  enum {
    O_NATURAL,
    O_RANDOM,
    O_CUTHILL_MCKEE,
  } ordering;

  ordering = O_NATURAL;
  if (orderingStr == "random")
    ordering = O_RANDOM;
  else if (orderingStr == "natural") {
  } else if (orderingStr == "cuthill-mckee" || orderingStr == "cm")
    ordering = O_CUTHILL_MCKEE;
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError, "Invalid ordering type");
  }

  // Get an ordering vector
  // NOTE: The orderingVector only orders *rows* of the matrix.  Off-proc columns
  // will get ignored in the aggregation phases, so we don't need to worry about
  // running off the end.
  Array<LO> orderingVector(numRows);
  for (LO i = 0; i < numRows; i++)
    orderingVector[i] = i;
  if (ordering == O_RANDOM)
    MueLu::NotayUtils::RandomReorder(orderingVector);
  else if (ordering == O_CUTHILL_MCKEE) {
    RCP<Xpetra::Vector<LO, LO, GO, NO> > rcmVector = MueLu::Utilities<SC, LO, GO, NO>::CuthillMcKee(*A);
    auto localVector                               = rcmVector->getData(0);
    for (LO i = 0; i < numRows; i++)
      orderingVector[i] = localVector[i];
  }

  // Get the party stated
  LO numNonAggregatedNodes = numRows, numDirichletNodes = 0;
  BuildInitialAggregates(pL, A, orderingVector(), kappa,
                         *aggregates, aggStat, numNonAggregatedNodes, numDirichletNodes);
  TEUCHOS_TEST_FOR_EXCEPTION(0 < numNonAggregatedNodes, Exceptions::RuntimeError,
                             "Initial pairwise aggregation failed to aggregate all nodes");
  LO numLocalAggregates = aggregates->GetNumAggregates();
  GetOStream(Statistics0) << "Init   : " << numLocalAggregates << " - "
                          << A->getLocalNumRows() / numLocalAggregates << std::endl;

  // Temporary data storage for further aggregation steps
  local_matrix_type intermediateP;
  local_matrix_type coarseLocalA;

  // Compute the on rank part of the local matrix
  // that the square submatrix that only contains
  // columns corresponding to local rows.
  LO numLocalDirichletNodes = numDirichletNodes;
  Array<LO> localVertex2AggId(aggregates->GetVertex2AggId()->getData(0).view(0, numRows));
  BuildOnRankLocalMatrix(A->getLocalMatrixDevice(), coarseLocalA);
  for (LO aggregationIter = 1; aggregationIter < maxNumIter; ++aggregationIter) {
    // Compute the intermediate prolongator
    BuildIntermediateProlongator(coarseLocalA.numRows(), numLocalDirichletNodes, numLocalAggregates,
                                 localVertex2AggId(), intermediateP);

    // Compute the coarse local matrix and coarse row sum
    BuildCoarseLocalMatrix(intermediateP, coarseLocalA);

    // Directly compute rowsum from A, rather than coarseA
    row_sum_type rowSum("rowSum", numLocalAggregates);
    {
      std::vector<std::vector<LO> > agg2vertex(numLocalAggregates);
      auto vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
      for (LO i = 0; i < (LO)numRows; i++) {
        if (aggStat[i] != AGGREGATED)
          continue;
        LO agg = vertex2AggId[i];
        agg2vertex[agg].push_back(i);
      }

      typename row_sum_type::HostMirror rowSum_h = Kokkos::create_mirror_view(rowSum);
      for (LO i = 0; i < numRows; i++) {
        // If not aggregated already, skip this guy
        if (aggStat[i] != AGGREGATED)
          continue;
        int agg                = vertex2AggId[i];
        std::vector<LO>& myagg = agg2vertex[agg];

        size_t nnz = A->getNumEntriesInLocalRow(i);
        ArrayView<const LO> indices;
        ArrayView<const SC> vals;
        A->getLocalRowView(i, indices, vals);

        SC mysum = Teuchos::ScalarTraits<Scalar>::zero();
        for (LO colidx = 0; colidx < static_cast<LO>(nnz); colidx++) {
          bool found = false;
          if (indices[colidx] < numRows) {
            for (LO j = 0; j < (LO)myagg.size(); j++)
              if (vertex2AggId[indices[colidx]] == agg)
                found = true;
          }
          if (!found) {
            *out << "- ADDING col " << indices[colidx] << " = " << vals[colidx] << std::endl;
            mysum += vals[colidx];
          } else {
            *out << "- NOT ADDING col " << indices[colidx] << " = " << vals[colidx] << std::endl;
          }
        }

        rowSum_h[agg] = mysum;
      }
      Kokkos::deep_copy(rowSum, rowSum_h);
    }

    // Get local orderingVector
    Array<LO> localOrderingVector(numRows);
    for (LO i = 0; i < numRows; i++)
      localOrderingVector[i] = i;
    if (ordering == O_RANDOM)
      MueLu::NotayUtils::RandomReorder(localOrderingVector);
    else if (ordering == O_CUTHILL_MCKEE) {
      RCP<Xpetra::Vector<LO, LO, GO, NO> > rcmVector = MueLu::Utilities<SC, LO, GO, NO>::CuthillMcKee(*A);
      auto localVector                               = rcmVector->getData(0);
      for (LO i = 0; i < numRows; i++)
        localOrderingVector[i] = localVector[i];
    }

    // Compute new aggregates
    numLocalAggregates    = 0;
    numNonAggregatedNodes = static_cast<LO>(coarseLocalA.numRows());
    std::vector<LO> localAggStat(numNonAggregatedNodes, READY);
    localVertex2AggId.resize(numNonAggregatedNodes, -1);
    BuildFurtherAggregates(pL, A, localOrderingVector, coarseLocalA, kappa, rowSum,
                           localAggStat, localVertex2AggId,
                           numLocalAggregates, numNonAggregatedNodes);

    // After the first initial pairwise aggregation
    // the Dirichlet nodes have been removed.
    numLocalDirichletNodes = 0;

    // Update the aggregates
    RCP<LOMultiVector> vertex2AggIdMV = aggregates->GetVertex2AggId();
    ArrayRCP<LO> vertex2AggId         = vertex2AggIdMV->getDataNonConst(0);
    for (size_t vertexIdx = 0; vertexIdx < A->getLocalNumRows(); ++vertexIdx) {
      LO oldAggIdx = vertex2AggId[vertexIdx];
      if (oldAggIdx != Teuchos::OrdinalTraits<LO>::invalid()) {
        vertex2AggId[vertexIdx] = localVertex2AggId[oldAggIdx];
      }
    }

    // We could probably print some better statistics at some point
    GetOStream(Statistics0) << "Iter " << aggregationIter << ": " << numLocalAggregates << " - "
                            << A->getLocalNumRows() / numLocalAggregates << std::endl;
  }
  aggregates->SetNumAggregates(numLocalAggregates);
  aggregates->AggregatesCrossProcessors(false);
  aggregates->ComputeAggregateSizes(true /*forceRecompute*/);

  // DO stuff
  Set(currentLevel, "Aggregates", aggregates);
  GetOStream(Statistics0) << aggregates->description() << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildInitialAggregates(const Teuchos::ParameterList& params,
                           const RCP<const Matrix>& A,
                           const Teuchos::ArrayView<const LO>& orderingVector,
                           const typename Teuchos::ScalarTraits<Scalar>::magnitudeType kappa,
                           Aggregates& aggregates,
                           std::vector<unsigned>& aggStat,
                           LO& numNonAggregatedNodes,
                           LO& numDirichletNodes) const {
  Monitor m(*this, "BuildInitialAggregates");
  using STS              = Teuchos::ScalarTraits<Scalar>;
  using MT               = typename STS::magnitudeType;
  using RealValuedVector = Xpetra::Vector<MT, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
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

  const MT kappa_init = kappa / (kappa - MT_TWO);
  const LO numRows    = aggStat.size();
  const int myRank    = A->getMap()->getComm()->getRank();

  // For finding "ties" where we fall back to the ordering.  Making this larger than
  // hard zero substantially increases code robustness.
  double tie_criterion = params.get<double>("aggregation: pairwise: tie threshold");
  double tie_less      = 1.0 - tie_criterion;
  double tie_more      = 1.0 + tie_criterion;

  // NOTE: Assumes 1 dof per node.  This constraint is enforced in Build(),
  // and so we're not doing again here.
  // This should probably be fixed at some point.

  // Extract diagonal, rowsums, etc
  // NOTE: The ghostedRowSum vector here has has the sign flipped from Notay's S
  RCP<Vector> ghostedDiag                = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixOverlappedDiagonal(*A);
  RCP<Vector> ghostedRowSum              = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixOverlappedDeletedRowsum(*A);
  RCP<RealValuedVector> ghostedAbsRowSum = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixOverlappedAbsDeletedRowsum(*A);
  const ArrayRCP<const SC> D             = ghostedDiag->getData(0);
  const ArrayRCP<const SC> S             = ghostedRowSum->getData(0);
  const ArrayRCP<const MT> AbsRs         = ghostedAbsRowSum->getData(0);

  // Aggregates stuff
  ArrayRCP<LO> vertex2AggId_rcp = aggregates.GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner_rcp   = aggregates.GetProcWinner()->getDataNonConst(0);
  ArrayView<LO> vertex2AggId    = vertex2AggId_rcp();
  ArrayView<LO> procWinner      = procWinner_rcp();

  // Algorithm 4.2

  // 0,1 : Initialize: Flag boundary conditions
  // Modification: We assume symmetry here aij = aji
  for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getLocalNumElements()); ++row) {
    MT aii    = STS::magnitude(D[row]);
    MT rowsum = AbsRs[row];

    if (aii >= kappa_init * rowsum) {
      *out << "Flagging index " << row << " as dirichlet "
                                          "aii >= kappa*rowsum = "
           << aii << " >= " << kappa_init << " " << rowsum << std::endl;
      aggStat[row] = IGNORED;
      --numNonAggregatedNodes;
      ++numDirichletNodes;
    }
  }

  // 2 : Iteration
  LO aggIndex = LO_ZERO;
  for (LO i = 0; i < numRows; i++) {
    LO current_idx = orderingVector[i];
    // If we're aggregated already, skip this guy
    if (aggStat[current_idx] != READY)
      continue;

    MT best_mu  = MT_ZERO;
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
      if (current_idx == col || col >= numRows || aggStat[col] != READY || val == SC_ZERO)
        continue;

      MT aij = STS::real(val);
      MT ajj = STS::real(D[col]);
      MT sj  = -STS::real(S[col]);  // NOTE: The ghostedRowSum vector here has has the sign flipped from Notay's S
      if (aii - si + ajj - sj >= MT_ZERO) {
        // Modification: We assume symmetry here aij = aji
        MT mu_top    = MT_TWO / (MT_ONE / aii + MT_ONE / ajj);
        MT mu_bottom = -aij + MT_ONE / (MT_ONE / (aii - si) + MT_ONE / (ajj - sj));
        MT mu        = mu_top / mu_bottom;

        // Modification: Explicitly check the tie criterion here
        if (mu > MT_ZERO && (best_idx == LO_INVALID || mu < best_mu * tie_less ||
                             (mu < best_mu * tie_more && orderingVector[col] < orderingVector[best_idx]))) {
          best_mu  = mu;
          best_idx = col;
          *out << "[" << current_idx << "] Column     UPDATED " << col << ": "
               << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
               << " = " << aii - si + ajj - sj << ", aij = " << aij << ", mu = " << mu << std::endl;
        } else {
          *out << "[" << current_idx << "] Column NOT UPDATED " << col << ": "
               << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
               << " = " << aii - si + ajj - sj << ", aij = " << aij << ", mu = " << mu << std::endl;
        }
      } else {
        *out << "[" << current_idx << "] Column     FAILED " << col << ": "
             << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
             << " = " << aii - si + ajj - sj << std::endl;
      }
    }  // end column for loop

    if (best_idx == LO_INVALID) {
      *out << "No node buddy found for index " << current_idx
           << " [agg " << aggIndex << "]\n"
           << std::endl;
      // We found no potential node-buddy, so let's just make this a singleton
      // NOTE: The behavior of what to do if you have no un-aggregated neighbors is not specified in
      // the paper

      aggStat[current_idx]      = ONEPT;
      vertex2AggId[current_idx] = aggIndex;
      procWinner[current_idx]   = myRank;
      numNonAggregatedNodes--;
      aggIndex++;

    } else {
      // We have a buddy, so aggregate, either as a singleton or as a pair, depending on mu
      if (best_mu <= kappa) {
        *out << "Node buddies (" << current_idx << "," << best_idx << ") [agg " << aggIndex << "]" << std::endl;

        aggStat[current_idx]      = AGGREGATED;
        aggStat[best_idx]         = AGGREGATED;
        vertex2AggId[current_idx] = aggIndex;
        vertex2AggId[best_idx]    = aggIndex;
        procWinner[current_idx]   = myRank;
        procWinner[best_idx]      = myRank;
        numNonAggregatedNodes -= 2;
        aggIndex++;

      } else {
        *out << "No buddy found for index " << current_idx << ","
                                                              " but aggregating as singleton [agg "
             << aggIndex << "]" << std::endl;

        aggStat[current_idx]      = ONEPT;
        vertex2AggId[current_idx] = aggIndex;
        procWinner[current_idx]   = myRank;
        numNonAggregatedNodes--;
        aggIndex++;
      }  // best_mu
    }    // best_idx
  }      // end Algorithm 4.2

  *out << "vertex2aggid :";
  for (int i = 0; i < static_cast<int>(vertex2AggId.size()); ++i) {
    *out << i << "(" << vertex2AggId[i] << ")";
  }
  *out << std::endl;

  // update aggregate object
  aggregates.SetNumAggregates(aggIndex);
}  // BuildInitialAggregates

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildFurtherAggregates(const Teuchos::ParameterList& params,
                           const RCP<const Matrix>& A,
                           const Teuchos::ArrayView<const LO>& orderingVector,
                           const typename Matrix::local_matrix_type& coarseA,
                           const typename Teuchos::ScalarTraits<Scalar>::magnitudeType kappa,
                           const Kokkos::View<typename Kokkos::ArithTraits<Scalar>::val_type*,
                                              Kokkos::LayoutLeft,
                                              typename Matrix::local_matrix_type::device_type>& rowSum,
                           std::vector<LocalOrdinal>& localAggStat,
                           Teuchos::Array<LocalOrdinal>& localVertex2AggID,
                           LO& numLocalAggregates,
                           LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildFurtherAggregates");

  // Set debug outputs based on environment variable
  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  using value_type             = typename local_matrix_type::value_type;
  const value_type KAT_zero    = Kokkos::ArithTraits<value_type>::zero();
  const magnitude_type MT_zero = Teuchos::ScalarTraits<magnitude_type>::zero();
  const magnitude_type MT_one  = Teuchos::ScalarTraits<magnitude_type>::one();
  const magnitude_type MT_two  = MT_one + MT_one;
  const LO LO_INVALID          = Teuchos::OrdinalTraits<LO>::invalid();

  // For finding "ties" where we fall back to the ordering.  Making this larger than
  // hard zero substantially increases code robustness.
  double tie_criterion = params.get<double>("aggregation: pairwise: tie threshold");
  double tie_less      = 1.0 - tie_criterion;
  double tie_more      = 1.0 + tie_criterion;

  typename row_sum_type::HostMirror rowSum_h = Kokkos::create_mirror_view(rowSum);
  Kokkos::deep_copy(rowSum_h, rowSum);

  // Extracting the diagonal of a KokkosSparse::CrsMatrix
  // is not currently provided in kokkos-kernels so here
  // is an ugly way to get that done...
  const LO numRows = static_cast<LO>(coarseA.numRows());
  typename local_matrix_type::values_type::HostMirror diagA_h("diagA host", numRows);
  typename local_matrix_type::row_map_type::HostMirror row_map_h = Kokkos::create_mirror_view(coarseA.graph.row_map);
  Kokkos::deep_copy(row_map_h, coarseA.graph.row_map);
  typename local_matrix_type::index_type::HostMirror entries_h = Kokkos::create_mirror_view(coarseA.graph.entries);
  Kokkos::deep_copy(entries_h, coarseA.graph.entries);
  typename local_matrix_type::values_type::HostMirror values_h = Kokkos::create_mirror_view(coarseA.values);
  Kokkos::deep_copy(values_h, coarseA.values);
  for (LO rowIdx = 0; rowIdx < numRows; ++rowIdx) {
    for (LO entryIdx = static_cast<LO>(row_map_h(rowIdx));
         entryIdx < static_cast<LO>(row_map_h(rowIdx + 1));
         ++entryIdx) {
      if (rowIdx == static_cast<LO>(entries_h(entryIdx))) {
        diagA_h(rowIdx) = values_h(entryIdx);
      }
    }
  }

  for (LO currentIdx = 0; currentIdx < numRows; ++currentIdx) {
    if (localAggStat[currentIdx] != READY) {
      continue;
    }

    LO bestIdx               = Teuchos::OrdinalTraits<LO>::invalid();
    magnitude_type best_mu   = Teuchos::ScalarTraits<magnitude_type>::zero();
    const magnitude_type aii = Teuchos::ScalarTraits<value_type>::real(diagA_h(currentIdx));
    const magnitude_type si  = Teuchos::ScalarTraits<value_type>::real(rowSum_h(currentIdx));
    for (auto entryIdx = row_map_h(currentIdx); entryIdx < row_map_h(currentIdx + 1); ++entryIdx) {
      const LO colIdx = static_cast<LO>(entries_h(entryIdx));
      if (currentIdx == colIdx || colIdx >= numRows || localAggStat[colIdx] != READY || values_h(entryIdx) == KAT_zero) {
        continue;
      }

      const magnitude_type aij = Teuchos::ScalarTraits<value_type>::real(values_h(entryIdx));
      const magnitude_type ajj = Teuchos::ScalarTraits<value_type>::real(diagA_h(colIdx));
      const magnitude_type sj  = -Teuchos::ScalarTraits<value_type>::real(rowSum_h(colIdx));  // NOTE: The ghostedRowSum vector here has has the sign flipped from Notay's S
      if (aii - si + ajj - sj >= MT_zero) {
        const magnitude_type mu_top    = MT_two / (MT_one / aii + MT_one / ajj);
        const magnitude_type mu_bottom = -aij + MT_one / (MT_one / (aii - si) + MT_one / (ajj - sj));
        const magnitude_type mu        = mu_top / mu_bottom;

        // Modification: Explicitly check the tie criterion here
        if (mu > MT_zero && (bestIdx == LO_INVALID || mu < best_mu * tie_less ||
                             (mu < best_mu * tie_more && orderingVector[colIdx] < orderingVector[bestIdx]))) {
          best_mu = mu;
          bestIdx = colIdx;
          *out << "[" << currentIdx << "] Column     UPDATED " << colIdx << ": "
               << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
               << " = " << aii - si + ajj - sj << ", aij = " << aij << " mu = " << mu << std::endl;
        } else {
          *out << "[" << currentIdx << "] Column NOT UPDATED " << colIdx << ": "
               << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
               << " = " << aii - si + ajj - sj << ", aij = " << aij << ", mu = " << mu << std::endl;
        }
      } else {
        *out << "[" << currentIdx << "] Column      FAILED " << colIdx << ": "
             << "aii - si + ajj - sj = " << aii << " - " << si << " + " << ajj << " - " << sj
             << " = " << aii - si + ajj - sj << std::endl;
      }
    }  // end loop over row entries

    if (bestIdx == Teuchos::OrdinalTraits<LO>::invalid()) {
      localAggStat[currentIdx]      = ONEPT;
      localVertex2AggID[currentIdx] = numLocalAggregates;
      --numNonAggregatedNodes;
      ++numLocalAggregates;
    } else {
      if (best_mu <= kappa) {
        *out << "Node buddies (" << currentIdx << "," << bestIdx << ") [agg " << numLocalAggregates << "]" << std::endl;

        localAggStat[currentIdx]      = AGGREGATED;
        localVertex2AggID[currentIdx] = numLocalAggregates;
        --numNonAggregatedNodes;

        localAggStat[bestIdx]      = AGGREGATED;
        localVertex2AggID[bestIdx] = numLocalAggregates;
        --numNonAggregatedNodes;

        ++numLocalAggregates;
      } else {
        *out << "No buddy found for index " << currentIdx << ","
                                                             " but aggregating as singleton [agg "
             << numLocalAggregates << "]" << std::endl;

        localAggStat[currentIdx]      = ONEPT;
        localVertex2AggID[currentIdx] = numLocalAggregates;
        --numNonAggregatedNodes;
        ++numLocalAggregates;
      }
    }
  }  // end loop over matrix rows

}  // BuildFurtherAggregates

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildOnRankLocalMatrix(const typename Matrix::local_matrix_type& localA,
                           typename Matrix::local_matrix_type& onrankA) const {
  Monitor m(*this, "BuildOnRankLocalMatrix");

  // Set debug outputs based on environment variable
  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  using local_graph_type = typename local_matrix_type::staticcrsgraph_type;
  using values_type      = typename local_matrix_type::values_type;
  using size_type        = typename local_graph_type::size_type;
  using col_index_type   = typename local_graph_type::data_type;
  using array_layout     = typename local_graph_type::array_layout;
  using memory_traits    = typename local_graph_type::memory_traits;
  using row_pointer_type = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
  using col_indices_type = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;
  // Extract on rank part of A
  // Simply check that the column index is less than the number of local rows
  // otherwise remove it.

  const int numRows = static_cast<int>(localA.numRows());
  row_pointer_type rowPtr("onrankA row pointer", numRows + 1);
  typename row_pointer_type::HostMirror rowPtr_h                   = Kokkos::create_mirror_view(rowPtr);
  typename local_graph_type::row_map_type::HostMirror origRowPtr_h = Kokkos::create_mirror_view(localA.graph.row_map);
  typename local_graph_type::entries_type::HostMirror origColind_h = Kokkos::create_mirror_view(localA.graph.entries);
  typename values_type::HostMirror origValues_h                    = Kokkos::create_mirror_view(localA.values);
  Kokkos::deep_copy(origRowPtr_h, localA.graph.row_map);
  Kokkos::deep_copy(origColind_h, localA.graph.entries);
  Kokkos::deep_copy(origValues_h, localA.values);

  // Compute the number of nnz entries per row
  rowPtr_h(0) = 0;
  for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
    for (size_type entryIdx = origRowPtr_h(rowIdx); entryIdx < origRowPtr_h(rowIdx + 1); ++entryIdx) {
      if (origColind_h(entryIdx) < numRows) {
        rowPtr_h(rowIdx + 1) += 1;
      }
    }
    rowPtr_h(rowIdx + 1) = rowPtr_h(rowIdx + 1) + rowPtr_h(rowIdx);
  }
  Kokkos::deep_copy(rowPtr, rowPtr_h);

  const LO nnzOnrankA = rowPtr_h(numRows);

  // Now use nnz per row to allocate matrix views and store column indices and values
  col_indices_type colInd("onrankA column indices", rowPtr_h(numRows));
  values_type values("onrankA values", rowPtr_h(numRows));
  typename col_indices_type::HostMirror colInd_h = Kokkos::create_mirror_view(colInd);
  typename values_type::HostMirror values_h      = Kokkos::create_mirror_view(values);
  int entriesInRow;
  for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
    entriesInRow = 0;
    for (size_type entryIdx = origRowPtr_h(rowIdx); entryIdx < origRowPtr_h(rowIdx + 1); ++entryIdx) {
      if (origColind_h(entryIdx) < numRows) {
        colInd_h(rowPtr_h(rowIdx) + entriesInRow) = origColind_h(entryIdx);
        values_h(rowPtr_h(rowIdx) + entriesInRow) = origValues_h(entryIdx);
        ++entriesInRow;
      }
    }
  }
  Kokkos::deep_copy(colInd, colInd_h);
  Kokkos::deep_copy(values, values_h);

  onrankA = local_matrix_type("onrankA", numRows, numRows,
                              nnzOnrankA, values, rowPtr, colInd);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildIntermediateProlongator(const LocalOrdinal numRows,
                                 const LocalOrdinal numDirichletNodes,
                                 const LocalOrdinal numLocalAggregates,
                                 const Teuchos::ArrayView<const LocalOrdinal>& localVertex2AggID,
                                 typename Matrix::local_matrix_type& intermediateP) const {
  Monitor m(*this, "BuildIntermediateProlongator");

  // Set debug outputs based on environment variable
  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_PAIRWISEAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  using local_graph_type = typename local_matrix_type::staticcrsgraph_type;
  using values_type      = typename local_matrix_type::values_type;
  using size_type        = typename local_graph_type::size_type;
  using col_index_type   = typename local_graph_type::data_type;
  using array_layout     = typename local_graph_type::array_layout;
  using memory_traits    = typename local_graph_type::memory_traits;
  using row_pointer_type = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
  using col_indices_type = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;

  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  const int intermediatePnnz = numRows - numDirichletNodes;
  row_pointer_type rowPtr("intermediateP row pointer", numRows + 1);
  col_indices_type colInd("intermediateP column indices", intermediatePnnz);
  values_type values("intermediateP values", intermediatePnnz);
  typename row_pointer_type::HostMirror rowPtr_h = Kokkos::create_mirror_view(rowPtr);
  typename col_indices_type::HostMirror colInd_h = Kokkos::create_mirror_view(colInd);

  rowPtr_h(0) = 0;
  for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
    // Skip Dirichlet nodes
    if (localVertex2AggID[rowIdx] == LO_INVALID) {
      rowPtr_h(rowIdx + 1) = rowPtr_h(rowIdx);
    } else {
      rowPtr_h(rowIdx + 1)       = rowPtr_h(rowIdx) + 1;
      colInd_h(rowPtr_h(rowIdx)) = localVertex2AggID[rowIdx];
    }
  }

  Kokkos::deep_copy(rowPtr, rowPtr_h);
  Kokkos::deep_copy(colInd, colInd_h);
  Kokkos::deep_copy(values, Kokkos::ArithTraits<typename values_type::value_type>::one());

  intermediateP = local_matrix_type("intermediateP",
                                    numRows, numLocalAggregates, intermediatePnnz,
                                    values, rowPtr, colInd);
}  // BuildIntermediateProlongator

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildCoarseLocalMatrix(const typename Matrix::local_matrix_type& intermediateP,
                           typename Matrix::local_matrix_type& coarseA) const {
  Monitor m(*this, "BuildCoarseLocalMatrix");

  using local_graph_type = typename local_matrix_type::staticcrsgraph_type;
  using values_type      = typename local_matrix_type::values_type;
  using size_type        = typename local_graph_type::size_type;
  using col_index_type   = typename local_graph_type::data_type;
  using array_layout     = typename local_graph_type::array_layout;
  using memory_traits    = typename local_graph_type::memory_traits;
  using row_pointer_type = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
  using col_indices_type = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;

  local_matrix_type AP;
  localSpGEMM(coarseA, intermediateP, "AP", AP);

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
  values_type valuesPt(Kokkos::ViewAllocateWithoutInitializing("Pt values"),
                       intermediateP.nnz());

  typename row_pointer_type::HostMirror rowPtrPt_h = Kokkos::create_mirror_view(rowPtrPt);
  typename col_indices_type::HostMirror entries_h  = Kokkos::create_mirror_view(intermediateP.graph.entries);
  Kokkos::deep_copy(entries_h, intermediateP.graph.entries);
  Kokkos::deep_copy(rowPtrPt_h, 0);
  for (size_type entryIdx = 0; entryIdx < intermediateP.nnz(); ++entryIdx) {
    rowPtrPt_h(entries_h(entryIdx) + 1) += 1;
  }
  for (LO rowIdx = 0; rowIdx < intermediateP.numCols(); ++rowIdx) {
    rowPtrPt_h(rowIdx + 1) += rowPtrPt_h(rowIdx);
  }
  Kokkos::deep_copy(rowPtrPt, rowPtrPt_h);

  typename row_pointer_type::HostMirror rowPtrP_h = Kokkos::create_mirror_view(intermediateP.graph.row_map);
  Kokkos::deep_copy(rowPtrP_h, intermediateP.graph.row_map);
  typename col_indices_type::HostMirror colIndP_h = Kokkos::create_mirror_view(intermediateP.graph.entries);
  Kokkos::deep_copy(colIndP_h, intermediateP.graph.entries);
  typename values_type::HostMirror valuesP_h = Kokkos::create_mirror_view(intermediateP.values);
  Kokkos::deep_copy(valuesP_h, intermediateP.values);
  typename col_indices_type::HostMirror colIndPt_h = Kokkos::create_mirror_view(colIndPt);
  typename values_type::HostMirror valuesPt_h      = Kokkos::create_mirror_view(valuesPt);
  const col_index_type invalidColumnIndex          = KokkosSparse::OrdinalTraits<col_index_type>::invalid();
  Kokkos::deep_copy(colIndPt_h, invalidColumnIndex);

  col_index_type colIdx = 0;
  for (LO rowIdx = 0; rowIdx < intermediateP.numRows(); ++rowIdx) {
    for (size_type entryIdxP = rowPtrP_h(rowIdx); entryIdxP < rowPtrP_h(rowIdx + 1); ++entryIdxP) {
      colIdx = entries_h(entryIdxP);
      for (size_type entryIdxPt = rowPtrPt_h(colIdx); entryIdxPt < rowPtrPt_h(colIdx + 1); ++entryIdxPt) {
        if (colIndPt_h(entryIdxPt) == invalidColumnIndex) {
          colIndPt_h(entryIdxPt) = rowIdx;
          valuesPt_h(entryIdxPt) = valuesP_h(entryIdxP);
          break;
        }
      }  // Loop over entries in row of Pt
    }    // Loop over entries in row of P
  }      // Loop over rows of P

  Kokkos::deep_copy(colIndPt, colIndPt_h);
  Kokkos::deep_copy(valuesPt, valuesPt_h);

  local_matrix_type intermediatePt("intermediatePt",
                                   intermediateP.numCols(),
                                   intermediateP.numRows(),
                                   intermediateP.nnz(),
                                   valuesPt, rowPtrPt, colIndPt);

  // Create views for coarseA matrix
  localSpGEMM(intermediatePt, AP, "coarseA", coarseA);
}  // BuildCoarseLocalMatrix

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    localSpGEMM(const typename Matrix::local_matrix_type& A,
                const typename Matrix::local_matrix_type& B,
                const std::string matrixLabel,
                typename Matrix::local_matrix_type& C) const {
  using local_graph_type = typename local_matrix_type::staticcrsgraph_type;
  using values_type      = typename local_matrix_type::values_type;
  using size_type        = typename local_graph_type::size_type;
  using col_index_type   = typename local_graph_type::data_type;
  using array_layout     = typename local_graph_type::array_layout;
  using memory_space     = typename device_type::memory_space;
  using memory_traits    = typename local_graph_type::memory_traits;
  using row_pointer_type = Kokkos::View<size_type*, array_layout, device_type, memory_traits>;
  using col_indices_type = Kokkos::View<col_index_type*, array_layout, device_type, memory_traits>;

  // Options
  int team_work_size = 16;
  std::string myalg("SPGEMM_KK_MEMORY");
  KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);
  KokkosKernels::Experimental::KokkosKernelsHandle<typename row_pointer_type::const_value_type,
                                                   typename col_indices_type::const_value_type,
                                                   typename values_type::const_value_type,
                                                   execution_space,
                                                   memory_space,
                                                   memory_space>
      kh;
  kh.create_spgemm_handle(alg_enum);
  kh.set_team_work_size(team_work_size);

  // Create views for AP matrix
  row_pointer_type rowPtrC(Kokkos::ViewAllocateWithoutInitializing("C row pointer"),
                           A.numRows() + 1);
  col_indices_type colIndC;
  values_type valuesC;

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

  // Numeric multiplication
  KokkosSparse::Experimental::spgemm_numeric(&kh, A.numRows(),
                                             B.numRows(), B.numCols(),
                                             A.graph.row_map, A.graph.entries, A.values, false,
                                             B.graph.row_map, B.graph.entries, B.values, false,
                                             rowPtrC, colIndC, valuesC);
  kh.destroy_spgemm_handle();

  C = local_matrix_type(matrixLabel, A.numRows(), B.numCols(), nnzC, valuesC, rowPtrC, colIndC);

}  // localSpGEMM

}  // namespace MueLu

#endif /* MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_ */
