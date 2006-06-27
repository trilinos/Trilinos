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

#ifndef _ispatest_epetra_test_utils_hpp_
#define _ispatest_epetra_test_utils_hpp_

#include <Isorropia_configdefs.hpp>

#ifdef HAVE_EPETRA

class Epetra_CrsGraph;
class Epetra_CrsMatrix;

/** ispatest is the namespace that contains isorropia's test-utilities.
  These test-utilities are for internal testing, and are not generally
  expected to be seen or used by external users.
*/
namespace ispatest {

/** Fill a matrix with the specified number of nonzeros per row, using
  matrix.InsertGlobalValues.
  Call FillComplete on the matrix before returning. If any negative error
  code is returned by an Epetra method, that will be the return value
  of this function.
*/
int fill_matrix(Epetra_CrsMatrix& matrix,
                int numNonzerosPerRow,
                bool verbose);

/** Fill a graph with the specified number of nonzeros per row.
  Call FillComplete on the graph before returning. If any non-zero error
  code is returned by an Epetra method, that will be the return value
  of this function.
*/
int fill_graph(Epetra_CrsGraph& graph,
                int numNonzerosPerRow,
                bool verbose);
}//namespace ispatest

#endif //HAVE_EPTERA

#endif

