// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  int partCount = 0;
  bool  globalProvidedParts;
  userProvidedParts_ = List.isParameter("partitioner: parts");
  globalProvidedParts = List.isParameter("partitioner: global ID parts");
  if (  userProvidedParts_)   partCount++;
  if (globalProvidedParts ) { partCount++;  userProvidedParts_ = true; }

  userProvidedMap_   = List.isParameter("partitioner: map");
  if (  userProvidedMap_)   partCount++;
  if (partCount > 1) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Ifpack2::UserPartitioner::setPartitionParameters: "
                               "you may specify only one of \"partitioner: parts\", \"partitioner: global ID parts\" and \"partitioner: map\","
                               " not two of them");
  }
  if (userProvidedMap_) {
    map_ = List.get ("partitioner: map", map_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      map_.is_null (), std::runtime_error, "Ifpack2::UserPartitioner::"
      "setPartitionParameters: map_ is null.");
  }
  if (userProvidedParts_) {
    typedef Teuchos::Array<typename Teuchos::ArrayRCP<typename GraphType::local_ordinal_type>> parts_type;
    typedef Teuchos::Array<typename Teuchos::ArrayRCP<typename GraphType::global_ordinal_type>> gparts_type;
    //FIXME JJH 9-Dec-2015 - This is a potentially expensive deep copy.  Use an ArrayRCP<ArrayRCP<int>> instead.

    // convert glboal ids to local ids before storing in Parts_
    if (globalProvidedParts ) { 
      gparts_type gparts;
      gparts = List.get<gparts_type>("partitioner: global ID parts");
      typedef Teuchos::RCP<   Tpetra::Map<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> const >  map_type;
      map_type OverlapMap  = List.get<map_type>("OverlapRowMap");
      this->Parts_.resize(gparts.size());
      for (int ii = 0; ii < (int)  gparts.size(); ii++) {
        this->Parts_[ii].resize(gparts[ii].size());
        for (int jj = 0; jj < (int)  gparts[ii].size(); jj++) {
          local_ordinal_type itmp =  (int) OverlapMap->getLocalElement( gparts[ii][jj] );
          if (itmp == Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) printf("global id %d (%d,%d) not found\n",(int) gparts[ii][jj],(int) ii, (int) jj);
          TEUCHOS_TEST_FOR_EXCEPTION(itmp == Teuchos::OrdinalTraits<local_ordinal_type>::invalid(), std::runtime_error,  " \"partitioner: global ID parts\" requires that all global IDs within a block reside on the MPI rank owning the block. This can be done using overlapping Schwarz with a BLOCK_RELAXATION subdomain solver making sure that \"schwarz: overlap level\" is sufficient. Note: zeroed out Dirichlet columns will never be included in overlap parts of domains.");
          this->Parts_[ii][jj] = itmp;
        }
      }
      List.remove("partitioner: global ID parts"); //JJH I observed that iterating through the ParameterList is much more
                                                   //expensive when "partitioner: collection" is a Parameter object.
                                                  //Thus, I remove it once it's fetched from the list.
    }
    else  {
      this->Parts_ = List.get<parts_type>("partitioner: parts");
      List.remove("partitioner: parts"); //JJH I observed that iterating through the ParameterList is much more
                                         //expensive when "partitioner: collection" is a Parameter object.
                                         //Thus, I remove it once it's fetched from the list.
    }
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
    const size_t localNumRows = this->Graph_->getLocalNumRows ();
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


