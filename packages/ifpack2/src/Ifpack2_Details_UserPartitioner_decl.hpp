// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  /// The ith entry in the ArrayRCP is the part (block) number
  /// that row i belongs to.  In this case, you are specifying OverlappingPartitioner::Partition_.
  ///     <li> "partitioner: parts" (Teuchos::Array<Teuchos::ArrayRCP<local ordinal>>)
  /// The i'th entry in the Array is an ArrayRCP that contains all the
  /// rows (local IDs) in part (block) i. In this case, you are specifying OverlappingPartitioner::Parts_.
  ///     <li> "partitioner: global ID parts" (Teuchos::Array<Teuchos::ArrayRCP<global ordinal>>)
  /// The i'th entry in the Array is an ArrayRCP that contains all the
  /// rows (global IDs) in part (block) i. In this case, you are specifying a translated version of 
  //  OverlappingPartitioner::Parts_.
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
