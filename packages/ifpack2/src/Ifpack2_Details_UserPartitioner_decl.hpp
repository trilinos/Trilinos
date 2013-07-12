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

/// \file Ifpack2_Details_Tpetra_RowGraph_decl.hpp
/// \brief Declaration of a user-defined partitioner in which the user can define a nonoverlapping partition of the graph.
/// \author Tom Benson
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It declares a user-defined partitioner to mirror the one in Ifpack.

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"

#include "Teuchos_ParameterList.hpp"

namespace Ifpack2{
/// \namespace Details
/// \brief Ifpack2 implementation details
///
/// This namespace contains implementation details of Ifpack2.
/// It is <i>not</i> meant for users.  Users should not rely on
/// anything in this namespace.

namespace Details{
/// \class UserPartitioner
/// \brief Allow the user to specify any nonoverlapping partition of the graph that they may choose.
///

template<class GraphType>
class UserPartitioner : public OverlappingPartitioner<GraphType> {

public:
  typedef typename GraphType::local_ordinal_type LocalOrdinal;
  typedef typename GraphType::global_ordinal_type GlobalOrdinal;
  typedef typename GraphType::node_type Node;

  //! Constructor.
  UserPartitioner(const Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph);

  //! Destructor.
  virtual ~UserPartitioner();

  //! Sets all the parameters for the partitioner.
  void setPartitionParameters(Teuchos::ParameterList& List);

  //! Computes the partitions. Returns 0 if successful.
  void computePartitions();

protected:
  Teuchos::ArrayRCP<LocalOrdinal> Map_;

};
}// namespace Details
}// namespace Ifpack2

#endif // IFPACK2_USER_PARTITIONER_DECL_HPP
