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

#ifndef IFPACK2_USER_PARTITIONER_DEF_HPP
#define IFPACK2_USER_PARTITIONER_DEF_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Details_UserPartitioner_decl.hpp"
#include "Ifpack2_OverlappingPartitioner_def.hpp"

namespace Ifpack2 {
  namespace Details {

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

}// namespace Details
}// namespace Ifpack2

#endif // IFPACK2_USERPARTITIONER_DEF_HPP


