// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATES_DEF_HPP
#define MUELU_AGGREGATES_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_LWGraph_kokkos.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Utilities_decl.hpp"
#include "MueLu_Aggregates_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::Aggregates(const LWGraph& graph) {
  numAggregates_       = 0;
  numGlobalAggregates_ = 0;

  vertex2AggId_ = LOMultiVectorFactory::Build(graph.GetImportMap(), 1);
  vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

  procWinner_ = LOVectorFactory::Build(graph.GetImportMap());
  procWinner_->putScalar(MUELU_UNASSIGNED);

  isRoot_ = Teuchos::ArrayRCP<bool>(graph.GetImportMap()->getLocalNumElements(), false);

  // slow but safe, force TentativePFactory to build column map for P itself
  aggregatesIncludeGhosts_ = true;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::
    Aggregates(LWGraph_kokkos graph) {
  numAggregates_       = 0;
  numGlobalAggregates_ = 0;

  vertex2AggId_ = LOMultiVectorFactory::Build(graph.GetImportMap(), 1);
  vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

  procWinner_ = LOVectorFactory::Build(graph.GetImportMap());
  procWinner_->putScalar(MUELU_UNASSIGNED);

  isRoot_ = Teuchos::ArrayRCP<bool>(graph.GetImportMap()->getLocalNumElements(), false);

  // slow but safe, force TentativePFactory to build column map for P itself
  aggregatesIncludeGhosts_ = true;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::
    Aggregates(const RCP<const Map>& map) {
  numAggregates_       = 0;
  numGlobalAggregates_ = 0;

  vertex2AggId_ = LOMultiVectorFactory::Build(map, 1);
  vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

  procWinner_ = LOVectorFactory::Build(map);
  procWinner_->putScalar(MUELU_UNASSIGNED);

  isRoot_ = Teuchos::ArrayRCP<bool>(map->getLocalNumElements(), false);

  // slow but safe, force TentativePFactory to build column map for P itself
  aggregatesIncludeGhosts_ = true;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Aggregates<LocalOrdinal, GlobalOrdinal, Node>::aggregates_sizes_type::const_type
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::ComputeAggregateSizes(bool forceRecompute) const {
  if (aggregateSizes_.size() && !forceRecompute) {
    return aggregateSizes_;

  } else {
    // It is necessary to initialize this to 0
    aggregates_sizes_type aggregateSizes("aggregates", numAggregates_);

    int myPID = GetMap()->getComm()->getRank();

    auto vertex2AggId = vertex2AggId_->getDeviceLocalView(Xpetra::Access::ReadOnly);
    auto procWinner   = procWinner_->getDeviceLocalView(Xpetra::Access::ReadOnly);

    typename AppendTrait<decltype(aggregateSizes_), Kokkos::Atomic>::type aggregateSizesAtomic = aggregateSizes;
    Kokkos::parallel_for(
        "MueLu:Aggregates:ComputeAggregateSizes:for", range_type(0, procWinner.size()),
        KOKKOS_LAMBDA(const LO i) {
          if (procWinner(i, 0) == myPID)
            aggregateSizesAtomic(vertex2AggId(i, 0))++;
        });

    aggregateSizes_ = aggregateSizes;

    return aggregateSizes;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ArrayRCP<LocalOrdinal>
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::
    ComputeAggregateSizesArrayRCP(bool forceRecompute) const {
  auto aggregateSizes = this->ComputeAggregateSizes(forceRecompute);

  // if this is the first time this is called, setup the host mirror and fill it
  if (!aggregateSizesHost_.is_allocated()) {
    aggregateSizesHost_ = Kokkos::create_mirror_view(aggregateSizes);
    Kokkos::deep_copy(aggregateSizesHost_, aggregateSizes);
  } else {
    // otherwise, only update if we forced a recompute
    if (forceRecompute)
      Kokkos::deep_copy(aggregateSizesHost_, aggregateSizes);
  }

  // put the data in an ArrayRCP, but do not give it ownership of the data
  Teuchos::ArrayRCP<LocalOrdinal> aggregateSizesArrayRCP(aggregateSizesHost_.data(), 0, aggregateSizesHost_.extent(0), false);

  return aggregateSizesArrayRCP;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Aggregates<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::GetGraph() const {
  using row_map_type = typename local_graph_type::row_map_type;
  using entries_type = typename local_graph_type::entries_type;
  using size_type    = typename local_graph_type::size_type;

  auto numAggregates = numAggregates_;

  if (static_cast<LO>(graph_.numRows()) == numAggregates)
    return graph_;

  auto vertex2AggId = vertex2AggId_->getDeviceLocalView(Xpetra::Access::ReadOnly);
  auto procWinner   = procWinner_->getDeviceLocalView(Xpetra::Access::ReadOnly);
  auto sizes        = ComputeAggregateSizes();

  // FIXME_KOKKOS: replace by ViewAllocateWithoutInitializing + rows(0) = 0.
  typename row_map_type::non_const_type rows("Agg_rows", numAggregates + 1);  // rows(0) = 0 automatically

  // parallel_scan (exclusive)
  Kokkos::parallel_scan(
      "MueLu:Aggregates:GetGraph:compute_rows", range_type(0, numAggregates),
      KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
        update += sizes(i);
        if (final_pass)
          rows(i + 1) = update;
      });

  decltype(rows) offsets(Kokkos::ViewAllocateWithoutInitializing("Agg_offsets"), numAggregates + 1);  // +1 is just for ease
  Kokkos::deep_copy(offsets, rows);

  int myPID = GetMap()->getComm()->getRank();

  size_type numNNZ;
  {
    Kokkos::View<size_type, device_type> numNNZ_device                    = Kokkos::subview(rows, numAggregates);
    typename Kokkos::View<size_type, device_type>::HostMirror numNNZ_host = Kokkos::create_mirror_view(numNNZ_device);
    Kokkos::deep_copy(numNNZ_host, numNNZ_device);
    numNNZ = numNNZ_host();
  }
  typename entries_type::non_const_type cols(Kokkos::ViewAllocateWithoutInitializing("Agg_cols"), numNNZ);
  size_t realnnz = 0;
  Kokkos::parallel_reduce(
      "MueLu:Aggregates:GetGraph:compute_cols", range_type(0, procWinner.size()),
      KOKKOS_LAMBDA(const LO i, size_t& nnz) {
        if (procWinner(i, 0) == myPID) {
          typedef typename std::remove_reference<decltype(offsets(0))>::type atomic_incr_type;
          auto idx  = Kokkos::atomic_fetch_add(&offsets(vertex2AggId(i, 0)), atomic_incr_type(1));
          cols(idx) = i;
          nnz++;
        }
      },
      realnnz);
  TEUCHOS_TEST_FOR_EXCEPTION(realnnz != numNNZ, Exceptions::RuntimeError,
                             "MueLu: Internal error: Something is wrong with aggregates graph construction: numNNZ = " << numNNZ << " != " << realnnz << " = realnnz");

  graph_ = local_graph_type(cols, rows);

  return graph_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void Aggregates<LocalOrdinal, GlobalOrdinal, Node>::ComputeNodesInAggregate(LO_view& aggPtr, LO_view& aggNodes, LO_view& unaggregated) const {
  LO numAggs                                          = GetNumAggregates();
  LO numNodes                                         = vertex2AggId_->getLocalLength();
  auto vertex2AggId                                   = vertex2AggId_->getDeviceLocalView(Xpetra::Access::ReadOnly);
  typename aggregates_sizes_type::const_type aggSizes = ComputeAggregateSizes(true);
  LO INVALID                                          = Teuchos::OrdinalTraits<LO>::invalid();

  aggPtr   = LO_view("aggPtr", numAggs + 1);
  aggNodes = LO_view("aggNodes", numNodes);
  LO_view aggCurr("agg curr", numAggs + 1);

  // Construct the "rowptr" and the counter
  Kokkos::parallel_scan(
      "MueLu:Aggregates:ComputeNodesInAggregate:scan", range_type(0, numAggs + 1),
      KOKKOS_LAMBDA(const LO aggIdx, LO& aggOffset, bool final_pass) {
        LO count = 0;
        if (aggIdx < numAggs)
          count = aggSizes(aggIdx);
        if (final_pass) {
          aggPtr(aggIdx)  = aggOffset;
          aggCurr(aggIdx) = aggOffset;
          if (aggIdx == numAggs)
            aggCurr(numAggs) = 0;  // use this for counting unaggregated nodes
        }
        aggOffset += count;
      });

  // Preallocate unaggregated to the correct size
  LO numUnaggregated = 0;
  Kokkos::parallel_reduce(
      "MueLu:Aggregates:ComputeNodesInAggregate:unaggregatedSize", range_type(0, numNodes),
      KOKKOS_LAMBDA(const LO nodeIdx, LO& count) {
        if (vertex2AggId(nodeIdx, 0) == INVALID)
          count++;
      },
      numUnaggregated);
  unaggregated = LO_view("unaggregated", numUnaggregated);

  // Stick the nodes in each aggregate's spot
  Kokkos::parallel_for(
      "MueLu:Aggregates:ComputeNodesInAggregate:for", range_type(0, numNodes),
      KOKKOS_LAMBDA(const LO nodeIdx) {
        LO aggIdx = vertex2AggId(nodeIdx, 0);
        if (aggIdx != INVALID) {
          // atomic postincrement aggCurr(aggIdx) each time
          aggNodes(Kokkos::atomic_fetch_add(&aggCurr(aggIdx), 1)) = nodeIdx;
        } else {
          // same, but using last entry of aggCurr for unaggregated nodes
          unaggregated(Kokkos::atomic_fetch_add(&aggCurr(numAggs), 1)) = nodeIdx;
        }
      });
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string Aggregates<LocalOrdinal, GlobalOrdinal, Node>::description() const {
  if (numGlobalAggregates_ == -1)
    return BaseClass::description() + "{nGlobalAggregates = not computed}";
  else
    return BaseClass::description() + "{nGlobalAggregates = " + toString(numGlobalAggregates_) + "}";
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void Aggregates<LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Statistics1) {
    if (numGlobalAggregates_ == -1)
      out0 << "Global number of aggregates: not computed " << std::endl;
    else
      out0 << "Global number of aggregates: " << numGlobalAggregates_ << std::endl;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal Aggregates<LocalOrdinal, GlobalOrdinal, Node>::GetNumGlobalAggregatesComputeIfNeeded() {
  if (numGlobalAggregates_ != -1) {
    LO nAggregates = GetNumAggregates();
    GO nGlobalAggregates;
    MueLu_sumAll(vertex2AggId_->getMap()->getComm(), (GO)nAggregates, nGlobalAggregates);
    SetNumGlobalAggregates(nGlobalAggregates);
  }
  return numGlobalAggregates_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
Aggregates<LocalOrdinal, GlobalOrdinal, Node>::GetMap() const {
  return vertex2AggId_->getMap();
}

}  // namespace MueLu

#endif  // MUELU_AGGREGATES_DEF_HPP
