/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
  global_size_t num_rows_per_proc = 10;
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  
  RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  prec.initialize();
  TEST_EQUALITY(0,0);
}

#ifdef LEAVE_THIS_STUFF_OFF
// Tests hypre's ability to work when A has a funky row map, but the vectors have
// a contiguous row map
TEUCHOS_UNIT_TEST( Ifpack_Hypre, ParameterList ){
  const double tol = 1e-7;

  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  Epetra_Map*            Map;
  // pointer to the matrix to be created
  Epetra_CrsMatrix*      Matrix;
  // container for parameters
  Teuchos::ParameterList GaleriList;
  // here we specify the global dimension of the problem
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  try
  {
    // Create the matrix
    Map = Galeri::CreateMap("Cartesian2D", Comm, GaleriList);
    Matrix   = Galeri::CreateCrsMatrix("Biharmonic2D", Map, GaleriList);

    // Create the parameter list
    Teuchos::ParameterList list("Preconditioner List");
    Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
    precList.set("HYPRE_BoomerAMGSetPrintLevel", 1);// print amg solution info
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
    Ifpack_Hypre preconditioner(Matrix);
    TEST_EQUALITY(preconditioner.SetParameters(list),0);
    TEST_EQUALITY(preconditioner.Compute(),0);

    // Create the RHS and solution vector
    int numVec = 2;
    Epetra_MultiVector X(preconditioner.OperatorDomainMap(), numVec);
    Epetra_MultiVector KnownX(preconditioner.OperatorDomainMap(), numVec);
    TEST_EQUALITY(KnownX.Random(),0);
    Epetra_MultiVector B(preconditioner.OperatorRangeMap(), numVec);
    TEST_EQUALITY(preconditioner.Apply(KnownX, B),0);

    TEST_EQUALITY(preconditioner.ApplyInverse(B,X),0);
    TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,Comm.NumProc())), true);

    //delete preconditioner;
    delete Map;
    delete Matrix;
  }
  catch (Galeri::Exception& rhs)
  {
    if (Comm.MyPID() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();
    }
  }
}


// Tests hypre's ability to work when A has a funky row map, but the vectors have
// a contiguous row map
TEUCHOS_UNIT_TEST( Ifpack_Hypre, ParameterList2 ){
  const double tol = 1e-7;

  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  Epetra_Map*            Map;
  // pointer to the matrix to be created
  Epetra_CrsMatrix*      Matrix;
  // container for parameters
  Teuchos::ParameterList GaleriList;
  // here we specify the global dimension of the problem
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  try
  {
    // Create the matrix
    Map = Galeri::CreateMap("Cartesian2D", Comm, GaleriList);
    Matrix   = Galeri::CreateCrsMatrix("Biharmonic2D", Map, GaleriList);

    // Create the parameter list
    Teuchos::ParameterList list("Preconditioner List");

    list.set("hypre: Solver", "PCG");
    Teuchos::ParameterList solverList = list.sublist("hypre: Solver functions");
    solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
    solverList.set("HYPRE_PCGSetTol", tol);                   // conv. tolerance
    solverList.set("HYPRE_PCGSetTwoNorm", 1);                  // use the two norm as the stopping criteria
    solverList.set("HYPRE_PCGSetPrintLevel", 0);               // print solve info
    solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
    list.set("hypre: Solver functions",solverList);

    list.set("hypre: Preconditioner", "BoomerAMG");
    Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
    precList.set("HYPRE_BoomerAMGSetPrintLevel", 1);// print amg solution info
    precList.set("HYPRE_BoomerAMGSetCoarsenType", 6);
    precList.set("HYPRE_BoomerAMGSetRelaxType", 6);  //Sym G.S./Jacobi hybrid
    precList.set("HYPRE_BoomerAMGSetNumSweeps", 1);
    precList.set("HYPRE_BoomerAMGSetTol", 0.0);      // conv. tolerance zero
    precList.set("HYPRE_BoomerAMGSetMaxIter", 1);   //do only one iteration!
    list.set("hypre: Preconditioner functions",precList);

    list.set("hypre: SolveOrPrecondition", "Solver");
    list.set("hypre: SetPreconditioner", true);

    // Create the preconditioner
    Ifpack_Hypre preconditioner(Matrix);
    TEST_EQUALITY(preconditioner.SetParameters(list),0);
    TEST_EQUALITY(preconditioner.Compute(),0);

    // Create the RHS and solution vector
    int numVec = 2;
    Epetra_MultiVector X(preconditioner.OperatorDomainMap(), numVec);
    Epetra_MultiVector KnownX(preconditioner.OperatorDomainMap(), numVec);
    TEST_EQUALITY(KnownX.Random(),0);
    Epetra_MultiVector B(preconditioner.OperatorRangeMap(), numVec);
    TEST_EQUALITY(preconditioner.Apply(KnownX, B),0);

    TEST_EQUALITY(preconditioner.ApplyInverse(B,X),0);
    TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,Comm.NumProc())), true);

    //delete preconditioner;
    delete Map;
    delete Matrix;
  }
  catch (Galeri::Exception& rhs)
  {
    if (Comm.MyPID() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();
    }
  }
}
#endif


// Tests the hypre interface's ability to work with both a preconditioner and linear
// solver via ApplyInverse
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  solverList.set("HYPRE_PCGSetPrintLevel", 0);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver functions",solverList);
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "ParaSails");
  list.set("hypre: SetPreconditioner", true);
  /*preconditioner.setParameters(list);
  using namespace Ifpack2;
    preconditioner.SetParameter(Hypre_Is_Solver, PCG); // Use a PCG Solver
  preconditioner.SetParameter(Hypre_Is_Preconditioner, ParaSails); // Use a ParaSails Preconditioner
  preconditioner.SetParameter(Hypre_Is_Solver, &HYPRE_ParCSRPCGSetMaxIter, 1000); // Set maximum iterations
  preconditioner.SetParameter(Hypre_Is_Solver, &HYPRE_ParCSRPCGSetTol, 1e-9); // Set a tolerance
  preconditioner.SetParameter(Hypre_Is_Solver); // Solve the problem
 
  preconditioner.SetParameter(true); // Use the preconditioner
  */
  preconditioner.compute();
  preconditioner.apply(B, X);


  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));

  //  int numIters;
  //  double residual;
  //  preconditioner.SetParameter(Solver, &HYPRE_ParCSRPCGGetNumIterations, &numIters);
  //  preconditioner.SetParameter(Solver, &HYPRE_ParCSRPCGGetFinalRelativeResidualNorm, &residual);
  preconditioner.CallFunctions();
}


#ifdef LEAVE_THIS_OFF
// This example uses contiguous maps, so hypre should not have problems
TEUCHOS_UNIT_TEST( Ifpack_Hypre, DiagonalMatrixInOrder ) {
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int numProcs = comm.NumProc();
  int myRank = comm.MyPID();

  //
  // Construct the contiguous map
  //
  Epetra_Map accMap(2*numProcs,0,comm);

  //
  // Construct the diagonal matrix
  //
  Epetra_CrsMatrix accMat(Copy,accMap,2,true);
  int col;
  double val;
  col = 2*myRank;
  val = 2*myRank+1;
  accMat.InsertGlobalValues(col,1,&val,&col);
  col = 2*myRank+1;
  val = 2*myRank+2;
  accMat.InsertGlobalValues(col,1,&val,&col);
  accMat.FillComplete();

  //
  // Create the initial guess and RHS
  //
  int numVec = 2;
  Epetra_MultiVector X(accMap, numVec);
  Epetra_MultiVector KnownX(accMap, numVec);
  KnownX.Random();
  Epetra_MultiVector B(accMap, numVec);
  accMat.Apply(KnownX, B);

  //
  // Create the parameter list
  //
  const double tol = 1e-7;
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
  Ifpack_Hypre prec(&accMat);
  TEST_EQUALITY(prec.SetParameters(list),0);
  TEST_EQUALITY(prec.Compute(),0);

  //
  // Solve the linear system
  //
  prec.apply(B,X);

  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,numProcs)), true);
}


// hypre does not like the distribution of the vectors in this example.
// Our interface should detect that and remap as necessary.
TEUCHOS_UNIT_TEST( Ifpack_Hypre, DiagonalMatrixOutOfOrder ) {
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int numProcs = comm.NumProc();
  int myRank = comm.MyPID();

  //
  // Construct the map [0 0 1 1]
  // Note that the rows will be out of order
  //
  int elementList[2];
  elementList[0] = 2*myRank+1;
  elementList[1] = 2*myRank;
  Epetra_Map ncMap(2*numProcs,2,elementList,0,comm);

  //
  // Construct the diagonal matrix
  //
  Epetra_CrsMatrix accMat(Copy,ncMap,2,true);
  int col;
  double val;
  col = 2*myRank+1;
  val = 2*myRank+2;
  accMat.InsertGlobalValues(col,1,&val,&col);
  col = 2*myRank;
  val = 2*myRank+1;
  accMat.InsertGlobalValues(col,1,&val,&col);
  accMat.FillComplete();

  //
  // Create the parameter list
  //
  const double tol = 1e-7;
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", tol);                   // conv. tolerance
  solverList.set("HYPRE_PCGSetTwoNorm", 1);                   // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 2);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver", "PCG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", false);


  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  Ifpack_Hypre prec(&accMat);
  TEST_EQUALITY(prec.SetParameters(list),0);
  TEST_EQUALITY(prec.Compute(),0);

  //
  // Create the initial guess and RHS
  // Note that these use non-contiguous maps
  //
  int numVec = 2;
  Epetra_MultiVector X(ncMap, numVec);
  Epetra_MultiVector KnownX(ncMap, numVec);
  KnownX.Random();
  Epetra_MultiVector B(ncMap, numVec);
  accMat.Apply(KnownX, B);

  //
  // Solve the linear system
  //
  TEST_EQUALITY(prec.ApplyInverse(B,X),0);

  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,numProcs)), true);
}



// Creates the identity matrix with a non-contiguous row map
// Even though the Epetra identity matrix has a map that hypre should not be happy with,
// hypre should be able to redistribute it.  It should also be able to accept the
// vectors we give it, since they're using the same distribution as the hypre matrix.
// This tests hypre's ability to perform as a linear solver via ApplyInverse.
TEUCHOS_UNIT_TEST( Ifpack_Hypre, NonContiguousRowMap ) {

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int numProcs = comm.NumProc();
  int myRank = comm.MyPID();

  //
  // Construct the map [0 1 0 1]
  //
  int elementList[2];
  elementList[0] = myRank;
  elementList[1] = myRank+numProcs;
  Epetra_Map badMap(2*numProcs,2,elementList,0,comm);
  Epetra_Map goodMap(2*numProcs,0,comm);

  //
  // Construct the identity matrix
  //
  Epetra_CrsMatrix badMat(Copy,badMap,2,true);
  int col;
  double val = 1.0;
  col = myRank;
  badMat.InsertGlobalValues(col,1,&val,&col);
  col = myRank+numProcs;
  badMat.InsertGlobalValues(col,1,&val,&col);
  badMat.FillComplete();

  //
  // Create the parameter list
  //
  const double tol = 1e-7;
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", tol);                   // conv. tolerance
  solverList.set("HYPRE_PCGSetTwoNorm", 1);                   // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 2);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver", "PCG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", false);

  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  Ifpack_Hypre prec(&badMat);
  TEST_EQUALITY(prec.SetParameters(list),0);
  TEST_EQUALITY(prec.Compute(),0);

  //
  // Create the initial guess and RHS
  //
  int numVec = 2;
  Epetra_MultiVector X(goodMap, numVec);
  Epetra_MultiVector RHS(goodMap, numVec);
  RHS.Random();

  //
  // Solve the linear system
  //
  TEST_EQUALITY(prec.ApplyInverse(RHS,X),0);

  Teuchos::ArrayRCP<const Scalar> v1v = X.get1dView();
  Teuchos::ArrayRCP<const Scalar> v2v = RHS.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,numProcs));

}
#endif

#define UNIT_TEST_GROUP_SC_LO_GO_NO(Scalar,LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, Construct, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, Apply, Scalar, LocalOrdinal,GlobalOrdinal,Node)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLGN_REAL( UNIT_TEST_GROUP_SC_LO_GO_NO )

} // namespace (anonymous)
