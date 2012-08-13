/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#include "AztecOO.h"
#include "AztecOO_Version.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

int main(int argc, char *argv[]) {
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  bool verbose = false;  // used to set verbose false on non-root processors
  bool verbose1 = false; // user's command-line argument
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') {
    verbose = true;
    if (Comm.MyPID()==0) verbose1 = true;
  }
  
  if (verbose1)
    cout << AztecOO_Version() << endl << endl;
  if (verbose)
    cout << Comm <<endl;
  
  long long NumGlobalElements = 5;
  int IndexBase = 0;

  Epetra_Map Map(NumGlobalElements, IndexBase, Comm);
  int NumMyElements = Map.NumMyElements();
  long long * MyGlobalElements = Map.MyGlobalElements64();

  Epetra_CrsMatrix A(Copy,Map,0);
  Epetra_Vector b(Map);
  Epetra_Vector x(Map);
  Epetra_Vector xx(Map);

  if (verbose)
    cout << "Proc = " << Comm.MyPID() << "  NumMyElements=" << NumMyElements << endl;
    
  for (int i=0; i<NumMyElements; i++) {
    long long index = MyGlobalElements[i]; // Global Diagonal location
    double value =  -pow(((double)10),((double) index));
    if (index==0) 
      value = 1.0; // First value will be positive 1, reminder will be negative powers of 10
    b[i] = value;  // RHS has same value as diagonal
    xx[i] = 1.0;   // Makes solution all ones
    x[i] = 0.0;     // Start with zero guess.
    A.InsertGlobalValues(index, 1, &value, &index);
  }

  A.FillComplete();	// Signal that data entries are complete.
   
  if (verbose1) cout << "A = " << endl;
  if (verbose) cout << A << endl;
  if (verbose1) cout << "xx = " << endl;
  if (verbose) cout << xx << endl;
  if (verbose1) cout << "x = " << endl;
  if (verbose) cout << x << endl;
  if (verbose1) cout << "b = " << endl;
  if (verbose) cout << b << endl;

  // Want to solve Ax=b
  Epetra_LinearProblem problem(&A, &x, &b);

  AztecOO solver(problem);

  solver.SetAztecOption(AZ_scaling, AZ_none);  // We want pure GMRES
  solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_max_iter, NumGlobalElements); // Set equal to global dimension
  solver.SetAztecOption(AZ_kspace, NumGlobalElements);
  solver.SetAztecOption(AZ_diagnostics, AZ_none);
  if (!verbose)
      solver.SetAztecOption(AZ_output, AZ_none);

  double single_error = 0.0;
  double double_error = 0.0;
  for (int i=0; i<5; i++) {
    
    if (i==0) {
      solver.SetAztecOption(AZ_orthog, AZ_single_classic);
      solver.Iterate(NumGlobalElements, 1.0E-14);
      single_error = solver.RecursiveResidual();
      
    }
    else if (i==1)  {
      solver.SetAztecOption(AZ_orthog, AZ_double_classic);
      solver.Iterate(NumGlobalElements, 1.0E-14);
      double_error = solver.RecursiveResidual();
      assert(double_error < single_error); // Error from double classic should be less than single
    }
    else if (i==2)  {
      solver.SetAztecOption(AZ_orthog, AZ_single_modified);
      solver.Iterate(NumGlobalElements, 1.0E-14);
      single_error = solver.RecursiveResidual();
    }
    else if (i==3)  {
      solver.SetAztecOption(AZ_orthog, AZ_double_modified);
      solver.Iterate(NumGlobalElements, 1.0E-14);
      double_error = solver.RecursiveResidual();
      assert(double_error < single_error); // Error from double classic should be less than single
    }
    else if (i==4)  {
      solver.SetAztecOption(AZ_solver, AZ_bicgstab);
      solver.Iterate(NumGlobalElements, 1.0E-14);
      assert(solver.RecursiveResidual()>single_error); // BiCGSTAB should always be worse than any GMRES answer
    }

    if (verbose1)
      cout << "Solver performed " << solver.NumIters() << " iterations." << endl
	   << "Norm of Recursive residual = " << solver.RecursiveResidual() << endl << endl;

    assert(solver.NumIters()==NumGlobalElements);
      

  // Print out the result
    if (verbose1) cout << "Computed Solution = " << endl;
    if (verbose) cout << x << endl;

    x.PutScalar(0.0);    
  }

  if (verbose1)
    cout << "Computed solution for 5x5 matrix diag(1, -10, -100, -1000, -10000) using the following unpreconditioned unscaled methods:" << endl
	 << "  1) GMRES with single step classical Gram-Schmidt" << endl
	 << "  2) GMRES with double step classical Gram-Schmidt" << endl
	 << "  3) GMRES with single step modified Gram-Schmidt" << endl
	 << "  4) GMRES with double step modified Gram-Schmidt" << endl
	 << "  5) BiCGSTAB" << endl << endl
	 << "Confirmed that double steps provide superior answer to single step" << endl
	 << "and that GMRES is superior to BiCGSTAB" << endl << endl;
  if (Comm.MyPID()==0)
    cout << "All tests passed." << endl;
    
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  return 0;
}

