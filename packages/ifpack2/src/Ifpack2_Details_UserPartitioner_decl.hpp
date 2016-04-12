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

#ifndef IFPACK2_USER_PARTITIONER_DECL_HPP
#define IFPACK2_USER_PARTITIONER_DECL_HPP

/// \file Ifpack2_Details_UserPartitioner_decl.hpp
/// \brief Declaration of a user-defined partitioner in which the user
///   can define a partition of the graph.  The partition may have locally
///   overlapping parts.
/// \author Tom Benson
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It declares a user-defined partitioner to mirror the one in Ifpack.

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"

namespace Ifpack2 {
namespace Details {

/// \class UserPartitioner
/// \brief Partition in which the user can define a nonoverlapping
///   partition of the graph in any way they choose.
///
/// See Ifpack2::Details::UserPartitioner::setPartitionParameters for a list of supported parameters.
/// \tparam GraphType Specialization of Tpetra::CrsGraph or
///   Tpetra::RowGraph.
template<class GraphType>
class UserPartitioner : public OverlappingPartitioner<GraphType> {
public:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> 
    row_graph_type;

  //! Constructor.
  UserPartitioner (const Teuchos::RCP<const row_graph_type>& graph);

  //! Destructor.
  virtual ~UserPartitioner();

  //! @brief Sets all the parameters for the partitioner.
  /// The only valid parameters are:
  ///   <ul>
  ///     <li> "partitioner: map" (Teuchos::ArrayRCP<local ordinal>)
  /// Use this option if the parts (aka blocks) do not overlap.  The ith entry in the ArrayRCP is the part (block) number
  /// that row i belongs to.  In this case, you are specifying OverlappingPartitioner::Partition_.
  ///     <li> "partitioner: parts" (Teuchos::Array<Teuchos::ArrayRCP<local ordinal>>)
  /// Use this option if the parts (aka blocks) overlap.  The i'th entry in the Array is an ArrayRCP that contains all the
  /// rows in part (block) i. In this case, you are specifying OverlappingPartitioner::Parts_.
  ///   </ul>
  /// You may set only one of these parameters.  Setting both will results in a runtime exception.
  void setPartitionParameters (Teuchos::ParameterList& List);

  //! Compute the partitions.
  void computePartitions ();

private:
  Teuchos::ArrayRCP<local_ordinal_type> map_;
  //! @brief True if user has provided list of parts.  False otherwise.
  bool userProvidedParts_;
  //! @brief True if user has provided map.  False otherwise.
  bool userProvidedMap_;
};

}// namespace Details
}// namespace Ifpack2

#endif // IFPACK2_USER_PARTITIONER_DECL_HPP
