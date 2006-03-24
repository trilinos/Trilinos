//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_Rebalance_hpp_
#define _Isorropia_Rebalance_hpp_

#include <Isorropia_configdefs.hpp>

#ifdef HAVE_EPETRA
class Epetra_Vector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA
/** Create a balanced copy of an input Epetra_CrsMatrix.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.

  The caller is responsible for deleting the returned matrix.
*/
Epetra_CrsMatrix* create_balanced_copy(const Epetra_CrsMatrix& input_matrix);
#endif

#ifdef HAVE_EPETRA
/** Create a balanced copy of an input Epetra_CrsMatrix, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_matrix.

  The caller is responsible for deleting the returned matrix.
*/
Epetra_CrsMatrix* create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                                       const Epetra_Vector& row_weights);
#endif //HAVE_EPETRA

#ifdef HAVE_EPETRA
/** Create a balanced copy of an input Epetra_CrsGraph.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.

  The caller is responsible for deleting the returned graph.
*/
Epetra_CrsGraph* create_balanced_copy(const Epetra_CrsGraph& input_graph);
#endif //HAVE_EPETRA

#ifdef HAVE_EPETRA
/** Create a balanced copy of an input Epetra_CrsGraph, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_graph.

  The caller is responsible for deleting the returned graph.
*/
Epetra_CrsGraph* create_balanced_copy(const Epetra_CrsGraph& input_graph,
                                       const Epetra_Vector& row_weights);
#endif //HAVE_EPETRA

}//namespace Isorropia

#endif

