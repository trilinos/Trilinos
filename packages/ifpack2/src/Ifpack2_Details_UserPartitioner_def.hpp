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

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Details_UserPartitioner_decl.hpp>
#include <Ifpack2_OverlappingPartitioner.hpp>

namespace Ifpack2 {
namespace Details {

template<class GraphType>
UserPartitioner<GraphType>::
UserPartitioner (const Teuchos::RCP<const row_graph_type>& graph) :
  OverlappingPartitioner<GraphType> (graph),
  userProvidedParts_(false),
  userProvidedMap_(false)
{}

template<class GraphType>
UserPartitioner<GraphType>::~UserPartitioner() {}

template<class GraphType>
void 
UserPartitioner<GraphType>::
setPartitionParameters (Teuchos::ParameterList& List)
{
  userProvidedParts_ = List.isParameter("partitioner: parts");
  userProvidedMap_   = List.isParameter("partitioner: map");
  if (userProvidedParts_ && userProvidedMap_) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Ifpack2::UserPartitioner::setPartitionParameters: "
                               "you may specify only one of \"partitioner: parts\" and \"partitioner: map\","
                               " not both");
  }
  if (userProvidedMap_) {
    map_ = List.get ("partitioner: map", map_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      map_.is_null (), std::runtime_error, "Ifpack2::UserPartitioner::"
      "setPartitionParameters: map_ is null.");
  }
  if (userProvidedParts_) {
    typedef Teuchos::Array<typename Teuchos::ArrayRCP<typename GraphType::local_ordinal_type>> parts_type;
    //FIXME JJH 9-Dec-2015 - This is a potentially expensive deep copy.  Use an ArrayRCP<ArrayRCP<int>> instead.
    this->Parts_ = List.get<parts_type>("partitioner: parts");
    List.remove("partitioner: parts"); //JJH I observed that iterating through the ParameterList is much more
                                       //expensive when "partitioner: collection" is a Parameter object.
                                       //Thus, I remove it once it's fetched from the list.
  }
}

template<class GraphType>
void UserPartitioner<GraphType>::computePartitions ()
{
  if (userProvidedParts_) {
    this->NumLocalParts_ = this->Parts_.size();
    //The user has explicitly defined Parts_, so there is no need to for Partition_.
    this->Partition_.resize(0);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      map_.is_null (), std::logic_error, "Ifpack2::UserPartitioner::"
      "computePartitions: map_ is null.");
    const size_t localNumRows = this->Graph_->getNodeNumRows ();
    for (size_t ii = 0; ii < localNumRows; ++ii) {
      this->Partition_[ii] = map_[ii];
    }
  }
  // IGNORING ALL THE SINGLETON STUFF IN IFPACK_USERPARTITIONER.CPP BY
  // ORDERS OF CHRIS SIEFERT
}

}// namespace Details
}// namespace Ifpack2

#endif // IFPACK2_USERPARTITIONER_DEF_HPP


