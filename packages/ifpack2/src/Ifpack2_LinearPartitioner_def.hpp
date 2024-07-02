// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_LINEAR_PARTITIONER_DEF_HPP
#define IFPACK2_LINEAR_PARTITIONER_DEF_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_LinearPartitioner_decl.hpp"

namespace Ifpack2 {

//==============================================================================
// Constructor
template<class GraphType>
LinearPartitioner<GraphType>::
LinearPartitioner (const Teuchos::RCP<const row_graph_type>& graph) :
  OverlappingPartitioner<GraphType> (graph)
{}


template<class GraphType>
LinearPartitioner<GraphType>::~LinearPartitioner() {}


template<class GraphType>
void
LinearPartitioner<GraphType>::
setPartitionParameters (Teuchos::ParameterList& /* List */) {}


template<class GraphType>
void LinearPartitioner<GraphType>::computePartitions()
{
  using Teuchos::as;
  // Partition_ is an array of local_ordinal_type.  local_ordinal_type
  // may be signed or unsigned.  NumLocalParts_ is int, and needs to
  // be signed, since negative values are significant.  Comparisons
  // between signed and unsigned integers often result in compiler
  // warnings, which is why we use as() for explicit conversions
  // below.  We also use as() because in a debug build, it checks for
  // overflow.
  const int mod = as<int> (this->Graph_->getLocalNumRows () / 
                           this->NumLocalParts_);
  for (size_t i = 0; i < this->Graph_->getLocalNumRows (); ++i) {
    this->Partition_[i] = as<local_ordinal_type> (i / mod);
    if (this->Partition_[i] >= as<local_ordinal_type> (this->NumLocalParts_)) {
      this->Partition_[i] = this->NumLocalParts_ - 1;
    }
  }
}


}// namespace Ifpack2

#define IFPACK2_LINEARPARTITIONER_INSTANT(LO,GO,N) \
  template class Ifpack2::LinearPartitioner<Tpetra::CrsGraph< LO, GO, N > >; \
  template class Ifpack2::LinearPartitioner<Tpetra::RowGraph< LO, GO, N > >;

#endif // IFPACK2_LINEARPARTITIONER_DEF_HPP
