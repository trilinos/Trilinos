/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_LINEPARTITIONER_DECL_HPP
#define IFPACK2_LINEPARTITIONER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner.hpp"
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
  typedef Tpetra::MultiVector<double,local_ordinal_type, global_ordinal_type, node_type>  multivector_type;


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
  void local_automatic_line_search(int NumEqns, Teuchos::ArrayView <local_ordinal_type> blockIndices, local_ordinal_type last, local_ordinal_type next, local_ordinal_type LineID, double tol, Teuchos::Array<local_ordinal_type> itemp, Teuchos::Array<double> dtemp) const;



  // User data
  int NumEqns_;
  Teuchos::RCP<multivector_type> coord_;
  double threshold_;

};

}// namespace Ifpack2

#endif // IFPACK2_LINEPARTITIONER_DECL_HPP
