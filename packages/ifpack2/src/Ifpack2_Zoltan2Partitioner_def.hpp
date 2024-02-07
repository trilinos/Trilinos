/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
  zoltan2AlgoName_ ("phg"), OverlappingPartitioner<GraphType> (graph)
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
