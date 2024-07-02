// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_ZOLTAN2_PARTITIONER_DEF_HPP
#define IFPACK2_ZOLTAN2_PARTITIONER_DEF_HPP

#if defined(HAVE_IFPACK2_ZOLTAN2)
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Zoltan2Partitioner_decl.hpp"

namespace Ifpack2 {

//==============================================================================
// Constructor
template<class GraphType>
Zoltan2Partitioner<GraphType>::
Zoltan2Partitioner (const Teuchos::RCP<const row_graph_type>& graph) :
  OverlappingPartitioner<GraphType> (graph), zoltan2AlgoName_ ("phg")
{}


template<class GraphType>
Zoltan2Partitioner<GraphType>::~Zoltan2Partitioner() {}


template<class GraphType>
void
Zoltan2Partitioner<GraphType>::
setPartitionParameters (Teuchos::ParameterList& List) {
  // Default is Parallel Hypergraph
  zoltan2AlgoName_ = List.get<std::string>("zoltan2: algorithm", zoltan2AlgoName_);
}

template<class GraphType>
void Zoltan2Partitioner<GraphType>::computePartitions()
{
  // Create an input adapter for the Tpetra matrix.
  Zoltan2GraphAdapterType zoltan2_graph(this->Graph_);

  // Specify partitioning parameters
  Teuchos::ParameterList zoltan2_params;
  zoltan2_params.set("partitioning_approach", "partition");
  zoltan2_params.set("num_local_parts", this->NumLocalParts_);
  if (zoltan2AlgoName_ == "parmetis") {
    zoltan2_params.set("algorithm", "parmetis");
    zoltan2_params.set("symmetrize_input", "transpose"); // not sure if this does anything, and may fail with non-symm graph
    zoltan2_params.set("partitioning_objective", "minimize_cut_edge_weight");
  } else {
    zoltan2_params.set("algorithm", zoltan2AlgoName_);
  }

  // Create and solve partitioning problem
  Zoltan2::PartitioningProblem<Zoltan2GraphAdapterType>
    problem(&zoltan2_graph, &zoltan2_params, this->Graph_->getComm());
  problem.solve();

  // Save partition
  auto parts = problem.getSolution().getPartListView();
  for (size_t i = 0; i < this->Graph_->getLocalNumRows (); ++i) {
    this->Partition_[i] = parts[i];
  }
}


}// namespace Ifpack2

#define IFPACK2_ZOLTAN2PARTITIONER_INSTANT(LO,GO,N) \
  template class Ifpack2::Zoltan2Partitioner<Tpetra::CrsGraph< LO, GO, N > >; \
  template class Ifpack2::Zoltan2Partitioner<Tpetra::RowGraph< LO, GO, N > >;

#endif // HAVE_IFPACK2_ZOLTAN2
#endif // IFPACK2_ZOLTAN2PARTITIONER_DEF_HPP
