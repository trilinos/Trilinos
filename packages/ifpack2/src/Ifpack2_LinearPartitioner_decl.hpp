/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_LINEARPARTITIONER_DECL_HPP
#define IFPACK2_LINEARPARTITIONER_DECL_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"

#include "Teuchos_ParameterList.hpp"

namespace Ifpack2{
//! Ifpack2::LinearPartitioner: A class to define linear partitions.

template<class GraphType>
class LinearPartitioner : public OverlappingPartitioner<GraphType> {

public:
  typedef typename GraphType::local_ordinal_type LocalOrdinal;
  typedef typename GraphType::global_ordinal_type GlobalOrdinal;
  typedef typename GraphType::node_type Node;

  //! Constructor.
  LinearPartitioner(const Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph);

  //! Destructor.
  virtual ~LinearPartitioner();

  //! Sets all the parameters for the partitioner (none for linear partioning).
  void setPartitionParameters(Teuchos::ParameterList& List);

  //! Computes the partitions. Returns 0 if successful.
  void computePartitions();

};
}// namespace Ifpack2

#endif // IFPACK2_LINEARPARTITIONER_DECL_HPP
