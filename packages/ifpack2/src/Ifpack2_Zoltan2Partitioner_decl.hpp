// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#endif // HAVE_IFPACK2_ZOLTAN2
#endif // IFPACK2_ZOLTAN2PARTITIONER_DECL_HPP
