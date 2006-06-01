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

#include <Isorropia_Rebalance.hpp>
#include <Isorropia_Zoltan_Rebalance.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra_utils.hpp>

#include <Teuchos_RefCountPtr.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA

Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
		     Teuchos::ParameterList& paramlist)
{
  std::string bal_package_str("Balancing package");
  std::string bal_package = paramlist.get(bal_package_str, "none_specified");
  if (bal_package == "Zoltan" || bal_package == "zoltan") {
#ifdef HAVE_EPETRAEXT_ZOLTAN
    return( Isorropia_Zoltan::create_balanced_copy(input_matrix, paramlist) );
#else
    throw Isorropia::Exception("Zoltan requested, but epetraext-zoltan not enabled.");
#endif
  }

  //first, create a weights vector which contains the number of nonzeros
  //per row in the input_matrix.
  Epetra_Vector* weights = 0;
  try {
    weights = Isorropia::Epetra_Utils::
       create_row_weights_nnz(input_matrix.Graph());
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //now, call the other overloading of 'create_balanced_copy', which
  //accepts a weights vector.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> balanced_matrix;
  try {
    balanced_matrix =
      Isorropia::create_balanced_copy(input_matrix, *weights);
    delete weights;
  }
  catch(std::exception& exc) {
    delete weights;
    throw exc;
  }

  return balanced_matrix;
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix)
{
  //first, create a weights vector which contains the number of nonzeros
  //per row in the input_matrix.
  Epetra_Vector* weights = 0;
  try {
    weights = Isorropia::Epetra_Utils::
       create_row_weights_nnz(input_matrix);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //now, call the other overloading of 'create_balanced_copy', which
  //accepts a weights vector.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> balanced_matrix;
  try {
    balanced_matrix =
      Isorropia::create_balanced_copy(input_matrix, *weights);
    delete weights;
  }
  catch(std::exception& exc) {
    delete weights;
    throw exc;
  }

  return balanced_matrix;
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     const Epetra_Vector& row_weights)
{
  //first, figure out what the balanced row-distribution should be, and
  //create a Epetra_Map that describes that row-distribution.
  const Epetra_Comm& comm = input_matrix.Comm();

  //bal_rowmap will be the result of the create_balanced_map function...
  Teuchos::RefCountPtr<Epetra_Map> bal_rowmap;
  try {
    bal_rowmap =
      Isorropia::Epetra_Utils::create_balanced_map(input_matrix.RowMap(),
                                                      row_weights);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //next, create a new Epetra_CrsMatrix (which will be the return-value of
  //this function) with the new row-distribution.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> balanced_matrix =
    Isorropia::redistribute_rows(input_matrix, *bal_rowmap);

  if (input_matrix.Filled()) {
    //If input_matrix.Filled(), the call FillComplete() on balanced_matrix.
    //Potential problem: what if matrix isn't square? Would it be
    //appropriate to use the domain-map and range-map from input_matrix??

    balanced_matrix->FillComplete();
  }

  return balanced_matrix;
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     const Epetra_Vector& row_weights)
{
  //first, figure out what the balanced row-distribution should be, and
  //create a Epetra_Map that describes that row-distribution.
  const Epetra_Comm& comm = input_matrix.Comm();

  //bal_rowmap that will be the result of the create_balanced_map function...
  Teuchos::RefCountPtr<Epetra_Map> bal_rowmap;
  try {
    bal_rowmap =
      Isorropia::Epetra_Utils::create_balanced_map(input_matrix.RowMatrixRowMap(),
						   row_weights);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //next, create a new Epetra_CrsMatrix (which will be the return-value of
  //this function) with the new row-distribution.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> balanced_matrix =
    Isorropia::redistribute_rows(input_matrix, *bal_rowmap);

  if (input_matrix.Filled()) {
    //If input_matrix.Filled(), the call FillComplete() on balanced_matrix.
    //Potential problem: what if matrix isn't square? Would it be
    //appropriate to use the domain-map and range-map from input_matrix??

    balanced_matrix->FillComplete();
  }

  return balanced_matrix;
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     Teuchos::ParameterList& paramlist)
{
  std::string bal_package_str("Balancing package");
  std::string bal_package = paramlist.get(bal_package_str, "none_specified");
  if (bal_package == "Zoltan" || bal_package == "zoltan") {
#ifdef HAVE_EPETRAEXT_ZOLTAN
    return( Isorropia_Zoltan::create_balanced_copy(input_graph, paramlist) );
#else
    throw Isorropia::Exception("Zoltan requested, but epetraext-zoltan not enabled.");
#endif
  }

  //first, create a weights vector which contains the number of nonzeros
  //per row in the input_graph.
  Epetra_Vector* weights = 0;
  try {
    weights = Isorropia::Epetra_Utils::create_row_weights_nnz(input_graph);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //now, call the other overloading of 'create_balanced_copy', which
  //accepts a weights vector.
  Teuchos::RefCountPtr<Epetra_CrsGraph> balanced_graph;
  try {
    balanced_graph =
      Isorropia::create_balanced_copy(input_graph, *weights);
    delete weights;
  }
  catch(std::exception& exc) {
    delete weights;
    throw exc;
  }

  return balanced_graph;
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     const Epetra_Vector& row_weights)
{
  //first, figure out what the balanced row-distribution should be, and
  //create a Epetra_Map that describes that row-distribution.
  const Epetra_Comm& comm = input_graph.Comm();
  //construct a dummy map that will be replaced by the result of the
  //create_balanced_map function...
  Teuchos::RefCountPtr<Epetra_Map> bal_rowmap;
  try {
    bal_rowmap =
      Isorropia::Epetra_Utils::create_balanced_map(input_graph.RowMap(),
                                                   row_weights);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //next, create a new Epetra_CrsGraph (which will be the return-value of
  //this function) with the new row-distribution.
  Teuchos::RefCountPtr<Epetra_CrsGraph> balanced_graph =
    Isorropia::redistribute_rows(input_graph, *bal_rowmap);

  if (input_graph.Filled()) {
    //If input_graph.Filled(), the call FillComplete() on balanced_graph.
    //Potential problem: what if graph isn't square? Would it be
    //appropriate to use the domain-map and range-map from input_graph??

    balanced_graph->FillComplete();
  }

  return balanced_graph;
}

Teuchos::RefCountPtr<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem)
{
  //first, create a weights vector which contains the number of nonzeros
  //per row in the input_graph.
  Epetra_Vector* weights = 0;
  try {
    weights =
      Isorropia::Epetra_Utils::create_row_weights_nnz(*input_problem.GetMatrix());
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //call another overloading of 'create_balanced_copy', which
  //accepts a weights vector and rebalances a matrix.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> balanced_matrix;
  try {
    balanced_matrix =
      Isorropia::create_balanced_copy(*input_problem.GetMatrix(), *weights);
  }
  catch(std::exception& exc) {
    delete weights;
    throw exc;
  }

  //redistribute input_problem's vector objects to match the row-distribution
  //of the new balanced_matrix.
  Teuchos::RefCountPtr<Epetra_MultiVector> x;
  Teuchos::RefCountPtr<Epetra_MultiVector> b;
  try {
    x = redistribute(*input_problem.GetLHS(), balanced_matrix->RowMap());
    b = redistribute(*input_problem.GetLHS(), balanced_matrix->RowMap());
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy(Epetra_LinearProblem): caught exception:");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  delete weights;

  Teuchos::RefCountPtr<Epetra_LinearProblem> linprob =
    Teuchos::rcp(new Epetra_LinearProblem(balanced_matrix.get(),
                                          x.get(), b.get()));

  //have the RefCountPtrs release ownership of the matrix and vectors,
  //otherwise those objects would be destroyed when this function returns.
  balanced_matrix.release();
  x.release();
  b.release();

  return( linprob );
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
redistribute_rows(const Epetra_CrsMatrix& input_matrix,
                  const Epetra_Map& target_rowmap,
                  Epetra_Import* importer)
{
  Epetra_Import* new_importer = 0;
  if (importer == 0) {
    new_importer = new Epetra_Import(target_rowmap, input_matrix.RowMap());
    importer = new_importer;
  }

  Teuchos::RefCountPtr<Epetra_CrsMatrix> target_matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, target_rowmap, 0));

  target_matrix->Import(input_matrix, *importer, Insert);

  //it is safe to delete new_importer even if it is NULL
  delete new_importer;

  return(target_matrix);
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
redistribute_rows(const Epetra_RowMatrix& input_matrix,
                  const Epetra_Map& target_rowmap,
                  Epetra_Import* importer)
{
  Epetra_Import* new_importer = 0;
  if (importer == 0) {
    new_importer = new Epetra_Import(target_rowmap, input_matrix.RowMatrixRowMap());
    importer = new_importer;
  }

  Teuchos::RefCountPtr<Epetra_CrsMatrix> target_matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, target_rowmap, 0));

  target_matrix->Import(input_matrix, *importer, Insert);

  //it is safe to delete new_importer even if it is NULL
  delete new_importer;

  return(target_matrix);
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
redistribute_rows(const Epetra_CrsGraph& input_graph,
                  const Epetra_Map& target_rowmap,
                  Epetra_Import* importer)
{
  Epetra_Import* new_importer = 0;
  if (importer == 0) {
    new_importer = new Epetra_Import(target_rowmap, input_graph.RowMap());
    importer = new_importer;
  }

  Teuchos::RefCountPtr<Epetra_CrsGraph> target_graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, target_rowmap, 0));

  target_graph->Import(input_graph, *importer, Insert);

  //it is safe to delete new_importer even if it is NULL
  delete new_importer;

  return(target_graph);
}

Teuchos::RefCountPtr<Epetra_MultiVector>
redistribute(const Epetra_MultiVector& input,
             const Epetra_BlockMap& target_map,
             Epetra_Import* importer)
{
  Epetra_Import* new_importer = 0;
  if (importer == 0) {
    new_importer = new Epetra_Import(target_map, input.Map());
    importer = new_importer;
  }

  Teuchos::RefCountPtr<Epetra_MultiVector> target_multivec =
    Teuchos::rcp(new Epetra_MultiVector(target_map, input.NumVectors(), false));

  target_multivec->Import(input, *importer, Insert);

  //it is safe to delete new_importer even if it is NULL
  delete new_importer;

  return(target_multivec);
}

Teuchos::RefCountPtr<Epetra_Vector>
redistribute(const Epetra_Vector& input,
             const Epetra_Map& target_map,
             Epetra_Import* importer)
{
  Epetra_Import* new_importer = 0;
  if (importer == 0) {
    new_importer = new Epetra_Import(target_map, input.Map());
    importer = new_importer;
  }

  Teuchos::RefCountPtr<Epetra_Vector> target_vec =
    Teuchos::rcp(new Epetra_Vector(target_map, false));

  target_vec->Import(input, *importer, Insert);

  //it is safe to delete new_importer even if it is NULL
  delete new_importer;

  return(target_vec);
}

#endif //HAVE_EPETRA

}//namespace Isorropia

