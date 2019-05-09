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

#include "Ifpack.h"
#include "Ifpack_Hypre.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "hypre_Helpers.hpp"

#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;


// Tests hypre interface's ability to initialize correctly
TEUCHOS_UNIT_TEST( Ifpack_Hypre, Construct ) {
  const int N = 10;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map RowMap(N, 0, comm);
  Epetra_CrsMatrix Matrix(Copy, RowMap, 1);
  for(int i = 0; i < N; i++){
    int indices[1];
    double values[1];
    indices[0] = i;
    values[0] = 1.0;
    Matrix.InsertGlobalValues(i, 1, values, indices);
  }
  Matrix.FillComplete();
  Ifpack_Hypre preconditioner(&Matrix);
  TEST_EQUALITY(preconditioner.Initialize(),0);
}


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
    RCP<FunctionParameter> functs[11];
    functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000));               // max iterations
    functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, tol));                   // conv. tolerance
    functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTwoNorm, 1));                  // use the two norm as the stopping criteria
    functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 0));               // print solve info
    functs[4] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1));                  // needed to get run info later
    functs[5] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 1)); // print amg solution info
    functs[6] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 6));
    functs[7] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6));  //Sym G.S./Jacobi hybrid
    functs[8] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, 1));
    functs[9] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, 0.0));      // conv. tolerance zero
    functs[10] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, 1));   //do only one iteration!

    list.set("Solver", PCG);
    list.set("Preconditioner", BoomerAMG);
    list.set("SolveOrPrecondition", Solver);
    list.set("SetPreconditioner", true);
    list.set("NumFunctions", 11);
    list.set<RCP<FunctionParameter>*>("Functions", functs);

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


// Tests the hypre interface's ability to work with both a preconditioner and linear
// solver via ApplyInverse
TEUCHOS_UNIT_TEST( Ifpack_Hypre, Ifpack ){
  const double tol = 1E-7;

  //
  // Create Laplace2D with a contiguous row distribution
  //
  Epetra_CrsMatrix* Crs_Matrix = newCrsMatrix(4);

  //
  // Create the hypre preconditioner
  //
  Ifpack Factory;
  RCP<Ifpack_Preconditioner> preconditioner= rcp(Factory.Create("Hypre", Crs_Matrix));
  TEST_EQUALITY(preconditioner->Initialize(), 0);
  int NumProc = Crs_Matrix->Comm().NumProc();
  int MyPID = Crs_Matrix->Comm().MyPID();

  //
  // Create the solution vector and RHS
  //
  int numVec = 2;
  Epetra_MultiVector X(preconditioner->OperatorRangeMap(), numVec);
  Epetra_MultiVector KnownX(preconditioner->OperatorRangeMap(), numVec);
  Epetra_MultiVector B(preconditioner->OperatorDomainMap(), numVec);
  TEST_EQUALITY(KnownX.Random(), 0);
  TEST_EQUALITY(preconditioner->Apply(KnownX, B), 0);

  Teuchos::ParameterList list("New List");
  RCP<FunctionParameter> functs[5];
  functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000));
  functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-9));
  functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1));
  functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 0));
  functs[4] = rcp(new FunctionParameter(Preconditioner, &HYPRE_ParaSailsSetLogging, 0));
  list.set("NumFunctions", 5);
  list.set<RCP<FunctionParameter>*>("Functions", functs);
  list.set("SolveOrPrecondition", Solver);
  list.set("Solver", PCG);
  list.set("Preconditioner", ParaSails);
  list.set("SetPreconditioner", true);
  TEST_EQUALITY(preconditioner->SetParameters(list), 0);
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, PCG); // Use a PCG Solver
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Preconditioner, ParaSails); // Use a ParaSails Preconditioner
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGSetMaxIter, 1000); // Set maximum iterations
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGSetTol, 1e-9); // Set a tolerance
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver); // Solve the problem
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(true); // Use the preconditioner
  TEST_EQUALITY(preconditioner->Compute(),0);
  preconditioner->Condest(Ifpack_Cheap, 1000, 1e-9, Crs_Matrix);
  TEST_EQUALITY(preconditioner->ApplyInverse(B, X),0);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,NumProc)), true);
  int numIters;
  double residual;
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGGetNumIterations, &numIters);
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGGetFinalRelativeResidualNorm, &residual);
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->CallFunctions();
  if(MyPID == 0) printf("It took %d iterations, and achieved %e residual.\n", numIters, residual);
  delete Crs_Matrix;
}


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
  RCP<FunctionParameter> functs[5];
  functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 100));               // max iterations
  functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, tol));                   // conv. tolerance
  functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTwoNorm, 1));                  // use the two norm as the stopping criteria
  functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 2));               // print solve info
  functs[4] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1));
  list.set("Solver", PCG);
  list.set("SolveOrPrecondition", Solver);
  list.set("SetPreconditioner", false);
  list.set("NumFunctions", 5);
  list.set<RCP<FunctionParameter>*>("Functions", functs);

  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  Ifpack_Hypre prec(&accMat);
  TEST_EQUALITY(prec.SetParameters(list),0);
  TEST_EQUALITY(prec.Compute(),0);

  //
  // Solve the linear system
  //
  TEST_EQUALITY(prec.ApplyInverse(B,X),0);

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
  RCP<FunctionParameter> functs[5];
  functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 100));               // max iterations
  functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, tol));                   // conv. tolerance
  functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTwoNorm, 1));                  // use the two norm as the stopping criteria
  functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 2));               // print solve info
  functs[4] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1));
  list.set("Solver", PCG);
  list.set("SolveOrPrecondition", Solver);
  list.set("SetPreconditioner", false);
  list.set("NumFunctions", 5);
  list.set<RCP<FunctionParameter>*>("Functions", functs);

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
  TEST_EQUALITY(prec.ApplyInverse(RHS,X),0);

  TEST_EQUALITY(EquivalentVectors(X, RHS, tol*10*pow(10.0,numProcs)), true);
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
  RCP<FunctionParameter> functs[11];
  functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1));               // max iterations
  functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, tol));                   // conv. tolerance
  functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTwoNorm, 1));                  // use the two norm as the stopping criteria
  functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 2));               // print solve info
  functs[4] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1));
  list.set("Solver", PCG);
  list.set("SolveOrPrecondition", Solver);
  list.set("SetPreconditioner", false);
  list.set("NumFunctions", 5);
  list.set<RCP<FunctionParameter>*>("Functions", functs);

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

  TEST_EQUALITY(EquivalentVectors(X, RHS, tol*10*pow(10.0,numProcs)), true);
}
