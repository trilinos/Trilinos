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

#include <ispatest_epetra_utils.hpp>
#include <Isorropia_Exception.hpp>

#ifdef HAVE_EPETRA

#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

/** ispatest is the namespace that contains isorropia's test-utilities
*/
namespace ispatest {

int fill_matrix(Epetra_CrsMatrix& matrix,
                int numNonzerosPerRow,
                bool verbose)
{
  int err = 0;
  const Epetra_Map& rowmap = matrix.RowMap();
  int num_my_rows = rowmap.NumMyElements();
  int num_global_rows = rowmap.NumGlobalElements();

  std::vector<int> indices(numNonzerosPerRow);
  std::vector<double> coefs(numNonzerosPerRow);

  for(int i=0; i<num_my_rows; ++i) {
    int global_row = rowmap.GID(i);
    int first_col = global_row - numNonzerosPerRow/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (num_global_rows - numNonzerosPerRow)) {
      first_col = num_global_rows - numNonzerosPerRow;
    }

    for(int j=0; j<numNonzerosPerRow; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    err = matrix.InsertGlobalValues(global_row, numNonzerosPerRow,
                                    &coefs[0], &indices[0]);
    if (err < 0) {
      err = matrix.ReplaceGlobalValues(global_row, numNonzerosPerRow,
                                      &coefs[0], &indices[0]);
      if (err < 0) {
        return(err);
      }
    }
  }

  err = matrix.FillComplete();
  return(err);
}

int fill_graph(Epetra_CrsGraph& graph,
                int numNonzerosPerRow,
                bool verbose)
{
  int err = 0;
  const Epetra_BlockMap& rowmap = graph.RowMap();
  int num_my_rows = rowmap.NumMyElements();
  int num_global_rows = rowmap.NumGlobalElements();

  std::vector<int> indices(numNonzerosPerRow);
  std::vector<double> coefs(numNonzerosPerRow);

  for(int i=0; i<num_my_rows; ++i) {
    int global_row = rowmap.GID(i);
    int first_col = global_row - numNonzerosPerRow/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (num_global_rows - numNonzerosPerRow)) {
      first_col = num_global_rows - numNonzerosPerRow;
    }

    for(int j=0; j<numNonzerosPerRow; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    err = graph.InsertGlobalIndices(global_row, numNonzerosPerRow,
                                     &indices[0]);
    if (err < 0) {
      return(err);
    }
  }

  err = graph.FillComplete();
  return(err);
}

}//namespace ispatest

#endif //HAVE_EPETRA

