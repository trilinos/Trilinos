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

#ifndef IFPACK2_USER_PARTITIONER_DEF_HPP
#define IFPACK2_USER_PARTITIONER_DEF_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_UserPartitioner_decl.hpp"
#include "Ifpack2_OverlappingPartitioner_def.hpp"

namespace Ifpack2 {

//==============================================================================
// Constructor
template<class GraphType>
UserPartitioner<GraphType>::UserPartitioner(const Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph) :
  OverlappingPartitioner<GraphType>(Graph)
{
}
//==============================================================================
// Destructor.
template<class GraphType>
UserPartitioner<GraphType>::~UserPartitioner()
{
}

//==============================================================================
  //! Sets all the parameters for the partitioner (none for linear partioning).
template<class GraphType>
void UserPartitioner<GraphType>::setPartitionParameters(Teuchos::ParameterList& List)
{
  Map_ = List.get("partitioner: map",Map_);
  if (Map_ == Teuchos::null)
    {
      throw std::runtime_error("AAAAARRRGGGGGHHHh!");

    }

}

//==============================================================================
template<class GraphType>
void UserPartitioner<GraphType>::UserPartitioner::computePartitions()
{
  
  if (Map_ == Teuchos::null)
    {
      // Throw some exception
      throw std::runtime_error("AAAAARRRGGGGGHHHh!");
    }

  for (size_t ii=0 ; ii < this->Graph_->getNodeNumRows() ; ++ii)
    {
      this->Partition_[ii] = Map_[ii];
    }

  //IGNORING ALL THE SINGLETON B.S. IN IFPACK_USERPARTITIONER.CPP BY
  //ORDERS OF CHRIS SIEFERT

}


}// namespace Ifpack2

#endif // IFPACK2_USERPARTITIONER_DEF_HPP

/* Have to assign a partition number to everyone, but don't want to assign partition number to velocity - we really just want things that are not partitioned in the initial partition. We can assign a dummy value for an initial partition (something like -1 or something greater than num partitions or something Teuchos::OrdinalTraits<LocalOrdinal>::invalid();*/
