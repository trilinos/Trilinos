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

#ifndef IFPACK2_ZOLTAN2PARTITIONER_DECL_HPP
#define IFPACK2_ZOLTAN2PARTITIONER_DECL_HPP

#if defined(HAVE_IFPACK2_ZOLTAN2)
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"

#include "Zoltan2_PartitioningProblem.hpp"
#include "Zoltan2_TpetraRowGraphAdapter.hpp"

namespace Ifpack2 {

/// \class Zoltan2Partitioner
/*! \brief A class to define Zoltan2-based partitions
    \tparam GraphType Specialization of Tpetra::RowGraph or Tpetra::CrsGraph.

    This class places the rows of the local graph into nonintersecting groups.  The number of groups is
    given by \c NumLocalParts_, a member of the base class OverlappingPartitioner.
    Partition is computed by calling Zoltan2
*/
template<class GraphType>
class Zoltan2Partitioner : public OverlappingPartitioner<GraphType> {
public:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> 
    row_graph_type;
  typedef Zoltan2::TpetraRowGraphAdapter<row_graph_type>
    Zoltan2GraphAdapterType;

  //! Constructor.
  Zoltan2Partitioner (const Teuchos::RCP<const row_graph_type>& graph);

  //! Destructor.
  virtual ~Zoltan2Partitioner ();

  //! Set the partitioner's parameters (none for linear partitioning).
  void setPartitionParameters (Teuchos::ParameterList& List);

  //! Compute the partitions.
  void computePartitions ();

private:
  std::string zoltan2AlgoName_;
};

}// namespace Ifpack2

#endif // HAVE_IFPACK2_XPETRA && HAVE_IFPACK2_ZOLTAN2
#endif // IFPACK2_ZOLTAN2PARTITIONER_DECL_HPP
