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

#ifndef _Isorropia_Epetra_hpp_
#define _Isorropia_Epetra_hpp_

#include <Isorropia_configdefs.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
class Epetra_Comm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_LinearProblem;
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {
  class Partitioner;

/** The Epetra namespace contains Isorropia's Epetra-specific
  functions.
*/
namespace Epetra {

#ifdef HAVE_EPETRA

/** Given an input matrix-graph that is to be repartitioned, and a parameter-
    list (possibly specifying the partitioning package/method etc.),
    create an instance of Isorropia::Partitioner. This is a factory
    function, the run-time type of the returned Partitioner is
    Isorropia::EpetraPartitioner.
*/
Teuchos::RefCountPtr<Isorropia::Partitioner>
create_partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
		   const Teuchos::ParameterList& paramlist);

/** Given an input row-matrix that is to be repartitioned, and a parameter-
    list (possibly specifying the partitioning package/method etc.),
    create an instance of Isorropia::Partitioner. This is a factory
    function, the run-time type of the returned Partitioner is
    Isorropia::EpetraPartitioner.
*/
Teuchos::RefCountPtr<Isorropia::Partitioner>
create_partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
		   const Teuchos::ParameterList& paramlist);

/** Given a Partitioner object, create a target map representing the
   new partitioning.
   This function calls partitioner.compute_partitioning() if it has not
   already been called.
*/
Teuchos::RefCountPtr<Epetra_Map>
create_target_map(const Epetra_Comm& comm, Partitioner& partitioner);

/** Create a balanced copy of an input Epetra_CrsMatrix.

  This function represents a basic case, not accepting weights.
  By default, it uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  If the ParameterList object contains a string parameter with
  name == "Balancing package" and value == "Zoltan", then the Zoltan
  library is used to perform the balancing.

  If the Zoltan package is specified, then the parameter-list will also
  be checked for Zoltan parameters such as "LB_METHOD". Valid values for
  the LB_METHOD parameter include "GRAPH" and "HYPERGRAPH", and others.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     Teuchos::ParameterList& paramlist);

/** Create a balanced copy of an input Epetra_RowMatrix.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix);

/** Create a balanced copy of an input Epetra_CrsMatrix, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_matrix.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                       const Epetra_Vector& row_weights);

/** Create a balanced copy of an input Epetra_RowMatrix, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_matrix.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                       const Epetra_Vector& row_weights);

/** Create a balanced copy of an input Epetra_CrsGraph.

  This function represents a basic case, not accepting weights.
  By default, it uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  If the ParameterList object contains a string parameter with
  name == "Balancing package" and value == "Zoltan", then the Zoltan
  library is used to perform the balancing.

  If the Zoltan package is specified, then the parameter-list will also
  be checked for Zoltan parameters such as "LB_METHOD". Valid values for
  the LB_METHOD parameter include "GRAPH" and "HYPERGRAPH", and others.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     Teuchos::ParameterList& paramlist);

/** Create a balanced copy of an input Epetra_CrsGraph, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_graph.
*/

Teuchos::RefCountPtr<Epetra_CrsGraph>
  create_balanced_copy(const Epetra_CrsGraph& input_graph,
                       const Epetra_Vector& row_weights);

/** Create a balanced copy of an input Epetra_LinearProblem.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros in the matrix equal in each partition. I.e., it is equivalent
  to specifying a weighted matrix rebalance where the weights are the number
  of nonzeros in each row. The vectors are then redistributed to match the
  row-distribution of the matrix.

  Important note: It is the caller's responsibility to destroy the
  matrix and vector attributes of the returned Epetra_LinearProblem.
*/
Teuchos::RefCountPtr<Epetra_LinearProblem>
  create_balanced_copy(const Epetra_LinearProblem& input_problem);

/** Return a new Epetra_CrsMatrix object constructed with target_rowmap,
  and with the contents of input_matrix imported into it.

  The caller is responsible for deleting the returned object.

  @param input_matrix Source/input object.

  @param target_rowmap Target rowmap, required to be compatible with
     input_matrix.RowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  redistribute_rows(const Epetra_CrsMatrix& input_matrix,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_CrsMatrix object constructed with target_rowmap,
  and with the contents of input_matrix imported into it.

  The caller is responsible for deleting the returned object.

  @param input_matrix Source/input object.

  @param target_rowmap Target rowmap, required to be compatible with
     input_matrix.RowMatrixRowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  redistribute_rows(const Epetra_RowMatrix& input_matrix,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_CrsGraph object constructed with target_rowmap,
  and with the contents of input_graph imported into it.

  @param input_graph Source/input object.

  @param target_rowmap Target rowmap, required to be compatible with
     input_graph.RowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_CrsGraph>
  redistribute_rows(const Epetra_CrsGraph& input_graph,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_MultiVector object constructed with target_map,
  and with the contents of 'input' imported into it.

  @param input Source/input object.

  @param target_map Target map, required to be compatible with
     input.Map() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_MultiVector>
  redistribute(const Epetra_MultiVector& input,
               const Epetra_BlockMap& target_map,
               Epetra_Import* importer=0);

/** Return a new Epetra_Vector object constructed with target_map,
  and with the contents of 'input' imported into it.

  @param input Source/input object.

  @param target_map Target map, required to be compatible with
     input.RowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_Vector>
  redistribute(const Epetra_Vector& input,
               const Epetra_Map& target_map,
               Epetra_Import* importer=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_matrix. The returned vector will have
  the same size and distribution as input_matrix's row-map.
*/
Epetra_Vector* create_row_weights_nnz(const Epetra_RowMatrix& input_matrix);

/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_graph. The returned vector will have
  the same size and distribution as input_graph's row-map.
*/
Epetra_Vector* create_row_weights_nnz(const Epetra_CrsGraph& input_graph);

/** Calculate a new partitioning, and fill output containers with new
    elements for the local partition, as well as export and import lists.

    \param input_map Input map describing the existing or 'old' partitioning.

    \param weights Input vector giving a weight for each element in input_map.
    weights.Map() is required to be the same size and layout as input_map.

    \param myNewElements Output vector containing all elements that will
    reside on the local partition in the new partitioning.

    \param exports Output map contains set of export elements, and maps them
    to the processors that they are to be exported to.

    \param imports Output map contains set of import elements, and maps them
    to the processors that they are to be imported from.

    \return Error-code, 0 if successful. This probably should be a void
    function, since a serious error will result in an exception-throw
    rather than an integer-code error-return.
*/
int
repartition(const Epetra_BlockMap& input_map,
	    const Epetra_Vector& weights,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports);

/** Given an Epetra_BlockMap object, fill a vector of length numprocs+1
  with each processor's starting offset into the Epetra_BlockMap's global
  set of elements (the last position will contain num-global-elements).
  Gather the vector of offsets onto all processors.
*/
void gather_all_proc_global_offsets(const Epetra_BlockMap& blkmap,
                                    std::vector<int>& all_proc_offsets);

#endif //DOXYGEN_SHOULD_SKIP_THIS
#endif //HAVE_EPETRA

}//namespace Epetra
}//namespace Isorropia

#endif

