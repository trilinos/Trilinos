// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_LINEPARTITIONER_DECL_HPP
#define IFPACK2_LINEPARTITIONER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Ifpack2 {

/// \class LinePartitioner
/*! \brief Ifpack2::LinePartitioner: A class to define partitions into a set of lines.

These "line" partitions could then be used in to do block Gauss-Seidel smoothing, for instance.

The current implementation uses a local line detection inspired by the work of Mavriplis
for convection-diffusion (AIAA Journal, Vol 37(10), 1999).

Here we use coordinate information to pick "close" points if they are sufficiently far away
from the "far" points.  We also make sure the line can never double back on itself, so that
the associated sub-matrix could (in theory) be handed off to a fast triangular solver.  This
implementation doesn't actual do that, however.

This implementation is derived from the related routine in ML.

Supported parameters:
  <ul>
    <li> \c "partitioner: line detection threshold": if \f$||x_j - x_i||^2 < thresh * \max_k||x_k - x_i||^2\f$, then the points are close enough to line smooth <Scalar>
    <li> \c "partitioner: coordinates"  : coordinates of local nodes  < Teuchos::MultiVector<double> >
    <li> \c "partitioner: PDE equations": number of equations per node <int>
  </ul>
*/
/// \tparam GraphType Specialization of Tpetra::RowGraph or Tpetra::CrsGraph.

  template<class GraphType,class Scalar>
class LinePartitioner : public OverlappingPartitioner<GraphType> {
public:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type>  row_graph_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
  typedef Tpetra::MultiVector<magnitude_type,local_ordinal_type, global_ordinal_type, node_type>  multivector_type;

  typedef typename row_graph_type::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename row_graph_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;

  //! Constructor.
  LinePartitioner(const Teuchos::RCP<const row_graph_type>& graph);

  //! Destructor.
  virtual ~LinePartitioner();

  //! Set the partitioner's parameters (none for linear partitioning).
  void setPartitionParameters(Teuchos::ParameterList& List);

  //! Compute the partitions.
  void computePartitions();

private:
    // Useful functions
  int Compute_Blocks_AutoLine(Teuchos::ArrayView<local_ordinal_type> blockIndices) const;
  void local_automatic_line_search(int NumEqns, Teuchos::ArrayView <local_ordinal_type> blockIndices, local_ordinal_type last, local_ordinal_type next, local_ordinal_type LineID, double tol, Teuchos::Array<local_ordinal_type> itemp, Teuchos::Array<magnitude_type> dtemp) const;



  // User data
  int NumEqns_;
  Teuchos::RCP<multivector_type> coord_;
  double threshold_;

};

}// namespace Ifpack2

#endif // IFPACK2_LINEPARTITIONER_DECL_HPP
