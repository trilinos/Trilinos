// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// Ifpack: Object-Oriented Algebraic Preconditioner Package
// Copyright (2002) NTESS
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Hypre.hpp>


#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include <string>
#include <stdio.h>
#include <map>
#include <HYPRE.h>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;
using Teuchos::RCP;
using Teuchos::rcp;


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, Construct, Scalar, LocalOrdinal, GlobalOrdinal, Node)
{
  if (!std::is_same<HYPRE_Real, Scalar>::value || !std::is_same<HYPRE_Int, GlobalOrdinal>::value)
    return;

  global_size_t num_rows_per_proc = 10;
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  
  RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  prec.initialize();
  TEST_EQUALITY(0,0);
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, BoomerAMG, Scalar, LocalOrdinal, GlobalOrdinal, Node)
{
  if (!std::is_same<HYPRE_Real, Scalar>::value || !std::is_same<HYPRE_Int, GlobalOrdinal>::value)
    return;

  const GlobalOrdinal num_rows_per_proc = 1000;
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  const double tol = 1e-7;

  auto rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  int NumProc = rowmap->getComm()->getSize();
  auto A = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap,-Teuchos::ScalarTraits<Scalar>::one());

  // Create the parameter list
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
  precList.set("HYPRE_BoomerAMGSetPrintLevel", 0);// print amg solution info
  precList.set("HYPRE_BoomerAMGSetCoarsenType", 6);
  precList.set("HYPRE_BoomerAMGSetRelaxType", 6);  //Sym G.S./Jacobi hybrid
  precList.set("HYPRE_BoomerAMGSetNumSweeps", 1);
  precList.set("HYPRE_BoomerAMGSetTol", 0.0);      // conv. tolerance zero
  precList.set("HYPRE_BoomerAMGSetMaxIter", 1);   //do only one iteration!
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "BoomerAMG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", true);
  
  // Create the preconditioner
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(A);
  preconditioner.setParameters(list);
  preconditioner.compute();
    
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.applyMat(KnownX, B);

  preconditioner.compute();
  preconditioner.apply(B, X);

  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, BoomerAMGNonContiguous, Scalar, LocalOrdinal, GlobalOrdinal, Node)
{
  if (!std::is_same<HYPRE_Real, Scalar>::value || !std::is_same<HYPRE_Int, GlobalOrdinal>::value)
    return;

  const GlobalOrdinal num_rows_per_proc = 1000;
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  const double tol = 1e-7;

  auto rowmap = tif_utest::create_odd_even_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  int NumProc = rowmap->getComm()->getSize();
  auto A = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap,-Teuchos::ScalarTraits<Scalar>::one());

  // Create the parameter list
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
  precList.set("HYPRE_BoomerAMGSetPrintLevel", 0);// print amg solution info
  precList.set("HYPRE_BoomerAMGSetCoarsenType", 6);
  precList.set("HYPRE_BoomerAMGSetRelaxType", 6);  //Sym G.S./Jacobi hybrid
  precList.set("HYPRE_BoomerAMGSetNumSweeps", 1);
  precList.set("HYPRE_BoomerAMGSetTol", 0.0);      // conv. tolerance zero
  precList.set("HYPRE_BoomerAMGSetMaxIter", 1);   //do only one iteration!
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "BoomerAMG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", true);
  

  // Create the preconditioner
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(A);
  preconditioner.setParameters(list);
  preconditioner.compute();
    
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.applyMat(KnownX, B);

  preconditioner.compute();
  preconditioner.apply(B, X);

  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}


// Tests the hypre interface's ability to work with both a preconditioner and linear
// solver via ApplyInverse
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  if (!std::is_same<HYPRE_Real, Scalar>::value || !std::is_same<HYPRE_Int, GlobalOrdinal>::value)
    return;

  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  const double tol = 1E-7;

  global_size_t num_rows_per_proc = 10;
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  
  RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
  int NumProc = rowmap->getComm()->getSize();

  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(crsmatrix);
  preconditioner.initialize();

  //
  // Create the solution vector and RHS
  //
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.compute();
  preconditioner.applyMat(KnownX, B);

  Teuchos::ParameterList list("New List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", 1e-9);                   // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 1);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 0);                  // needed to get run info later
  list.set("hypre: Solver functions",solverList);
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "ParaSails");
  list.set("hypre: SetPreconditioner", true);
  
  preconditioner.compute();
  preconditioner.apply(B, X);

  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}



TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, DiagonalMatrixInOrder, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  if (!std::is_same<HYPRE_Real, Scalar>::value || !std::is_same<HYPRE_Int, GlobalOrdinal>::value)
    return;

  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  const double tol = 1E-7;

  LocalOrdinal num_rows_per_proc = 10;
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsgraph);

  int NumProc = crsmatrix->getMap()->getComm()->getSize();

  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(crsmatrix);
  preconditioner.initialize();


  //
  // Create the initial guess and RHS
  //
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.compute();
  preconditioner.applyMat(KnownX, B);

  //
  // Create the parameter list
  //
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", tol);                    // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 0);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver", "PCG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", false);

  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  preconditioner.setParameters(list);
  preconditioner.compute();

  //
  // Solve the linear system
  //
  preconditioner.apply(B,X);
  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, DiagonalMatrixNonContiguous, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  if (!std::is_same<HYPRE_Real, Scalar>::value || !std::is_same<HYPRE_Int, GlobalOrdinal>::value)
    return;

  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  const double tol = 1E-7;

  LocalOrdinal num_rows_per_proc = 10;
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_odd_even_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsgraph);

  int NumProc = crsmatrix->getMap()->getComm()->getSize();

  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(crsmatrix);
  preconditioner.initialize();


  //
  // Create the initial guess and RHS
  //
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.compute();
  preconditioner.applyMat(KnownX, B);

  //
  // Create the parameter list
  //
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", tol);                    // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 0);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver", "PCG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", false);

  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  preconditioner.setParameters(list);
  preconditioner.compute();

  //
  // Solve the linear system
  //
  preconditioner.apply(B,X);
  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}


#define UNIT_TEST_GROUP_SC_LO_GO_NO(Scalar,LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, Construct, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, Apply, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, DiagonalMatrixInOrder, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, DiagonalMatrixNonContiguous, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, BoomerAMG, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, BoomerAMGNonContiguous, Scalar, LocalOrdinal,GlobalOrdinal,Node) 


#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLGN( UNIT_TEST_GROUP_SC_LO_GO_NO )

} // namespace (anonymous)
