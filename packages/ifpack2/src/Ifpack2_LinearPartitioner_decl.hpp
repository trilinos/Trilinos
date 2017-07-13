/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_LINEARPARTITIONER_DECL_HPP
#define IFPACK2_LINEARPARTITIONER_DECL_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"

namespace Ifpack2 {

/// \class LinearPartitioner
/*! \brief A class to define linear partitions
    \tparam GraphType Specialization of Tpetra::RowGraph or Tpetra::CrsGraph.

    This class places the rows of the local graph into nonintersecting groups.  The number of groups is
    given by \c NumLocalParts_, a member of the base class OverlappingPartitioner.  Local row i is placed
    in the group floor(i/NumLocalParts_), with the exception that any leftover rows are placed in the largest
    group NumLocalParts_-1.
*/
template<class GraphType>
class LinearPartitioner : public OverlappingPartitioner<GraphType> {
public:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> 
    row_graph_type;

  //! Constructor.
  LinearPartitioner (const Teuchos::RCP<const row_graph_type>& graph);

  //! Destructor.
  virtual ~LinearPartitioner ();

  //! Set the partitioner's parameters (none for linear partitioning).
  void setPartitionParameters (Teuchos::ParameterList& List);

  //! Compute the partitions.
  void computePartitions ();
};

}// namespace Ifpack2

#endif // IFPACK2_LINEARPARTITIONER_DECL_HPP
