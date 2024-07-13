// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
