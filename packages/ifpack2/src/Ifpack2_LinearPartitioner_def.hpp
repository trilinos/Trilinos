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
  const int mod = as<int> (this->Graph_->getNodeNumRows () / 
                           this->NumLocalParts_);
  for (size_t i = 0; i < this->Graph_->getNodeNumRows (); ++i) {
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
