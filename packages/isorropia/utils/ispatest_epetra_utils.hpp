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

************************************************************************
*/
//@HEADER

#ifndef _ispatest_epetra_test_utils_hpp_
#define _ispatest_epetra_test_utils_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <string>
#include <vector>
#include <ostream>

#ifdef HAVE_EPETRA

class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Comm;

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

/** Verify that a matrix is a valid Epetra_CrsMatrix by attempting
  to multiply with it.  Return true if successful, false otherwise.

  A is the input matrix, and "x" and "y" will be created and then
  y = Ax will be computed.

  "x" will be created with the domain map of the "A",
  "y" will be created with the range map.

  Return true if successful, false otherwise.
*/
bool test_matrix_vector_multiply(Epetra_CrsMatrix &A);

/** Verify that a matrix is a valid Epetra_RowMatrix by attempting
  to multiply with it.  Return true if successful, false otherwise.

  y = Rx will be calculated, where "y" will be created with the
  OperatorRangeMap() of "R" and "x" will be created with 
  the OperatorDomainMap() of "R".

  Return true if successful, false otherwise.
*/
bool test_row_matrix_vector_multiply(Epetra_RowMatrix &R);

/** Verify that a matrix is a valid Epetra_CrsGraph by attempting
  to multiply with it.  Return true if successful, false otherwise.

  An Epetra_CrsMatrix (A) will be created from the input graph G, and
  then y = Ax will be computed.

  "x" will be created with the domain map of "G", and "y" will be 
  created with the range map of "G".

  Return true if successful, false otherwise.
*/
bool test_matrix_vector_multiply(Epetra_CrsGraph &G);

/** Verify that the Epetra_RowMatrix in a Epetra_LinearProblem is
    valid by attempting to multiply with it.

  An Epetra_RowMatrix (R) will be extracted from the linear
  problem, and then y = Rx will be computed.

  "y" will be created with the OperatorRangeMap() of "R" and "x" 
   will be created with the OperatorDomainMap() of "R".

  This test does not use the vectors in the linear problem, it uses
  a made up vector.  (Maybe it should use rhs and lhs?)

  Return true if successful, false otherwise.
*/
bool test_matrix_vector_multiply(Epetra_LinearProblem &LP);

/** Method to create an Epetra_Map from an Epetra_BlockMap, when using
    methods that require an Epetra_Map but you only have an Epetra_BlockMap.
 
    Caller should delete the returned map when done with it.
*/
Epetra_Map *map_from_blockmap(const Epetra_BlockMap &b);

/** Read in the file "fname".  It should be a text file with a list of
    1, 2 or 3-dimensional floating point coordinates.  Each point is
    on one line of the file, and coordinates are separated by white space.
    Blank lines are allowed.

    Fill the vectors x, y, and z with the x, y and z coordinates.  (z will
    have size 0 for 1 and 2 dimensional coordinates, y will have size 0
    for 1 dimensional coordinates.)

    Return the dimension of the coordinates, or 0 if none could be read.
 */

int readCoordFile(const std::string &fname,
         std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);

/** Create a multivector of 1 dimensional weights with the same
    distribution as the supplied map.  The weight for each object
    will be computed with the supplied function.  The arguments for the
    supplied function are the global ID of the object, the rank of the
    calling process, the global number of objects, and the number of
    processes.
*/
Epetra_MultiVector *makeWeights(const Epetra_BlockMap &map, double (*wFunc)(const int, const int, const int, const int));

/** A function for makeWeights, returns 1.0 for the object's weight.
  */
double unitWeights(const int id, const int me, const int nids, const int nprocs);

/** A function for makeWeights, weight is based on global ID, lower in the middle and
    higher at both ends.  Graph of weights for IDs would be "v" shaped.
  */
double veeWeights(const int id, const int me, const int nids, const int nprocs);

/** A function for makeWeights, weights alternate 1.0, 2.0, 1.0, etc. by process ID.
  */
double alternateWeights(const int id, const int me, const int nids, const int nprocs);

/** Read in a file of 1, 2 or 3 dimensional coordinates and create a
    multivector with a standard linear distribution.
*/
Epetra_MultiVector *file2multivector(const Epetra_Comm &comm, const std::string &fname);

/** Print out the contents of the multivector by process, with optional title.
    This is more compact the Epetra_MultiVector::Print().
  */
int printMultiVector(const Epetra_MultiVector mv, std::ostream &os, 
                     const char *s="MultiVector contents");

} //namespace ispatest

#endif //HAVE_EPETRA

#endif

