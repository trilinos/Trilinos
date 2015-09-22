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
#include "AztecOO_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "AztecOO_Operator.h"
#include "Epetra_InvOperator.h"
#include <string>
// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) { 
  if (fabs((x-y)/x) > 0.01) {
    return(1); 
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
  }
  else {
    if (verbose) cout << message << " check OK." << endl;    
    return(0);
  }
}

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int MyPID = comm.MyPID();

  bool verbose = false;
  bool verbose1 = false; 
  // Check if we should print results to standard out
  if (argc > 1) {
    if ((argv[1][0] == '-') && (argv[1][1] == 'v')) {
      verbose1 = true;
      if (MyPID==0) verbose = true;
    }
  }

  if (verbose1) cout << comm << endl;


  // Uncomment the next three lines to debug in mpi mode
  //int tmp;
  //if (MyPID==0) cin >> tmp;
  //comm.Barrier();

  Epetra_CrsMatrix * A; 
  EPETRA_CHK_ERR(EpetraExt::MatlabFileToCrsMatrix("A.dat", comm, A));

  Epetra_Vector  x(A->OperatorDomainMap()); 
  Epetra_Vector  b(A->OperatorRangeMap());
  x.Random();
  A->Apply(x,b); // Generate RHS from x
  Epetra_Vector xx(x); // Copy x to xx for later use

  Epetra_LinearProblem problem(A, &x, &b);
  // Construct a solver object for this problem

  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond, AZ_none);
  if (!verbose1) solver.SetAztecOption(AZ_output, AZ_none);
  solver.SetAztecOption(AZ_kspace, A->NumGlobalRows());
  AztecOO_Operator AOpInv(&solver, A->NumGlobalRows());
  Epetra_InvOperator AInvOp(&AOpInv);

  EPETRA_CHK_ERR(EpetraExt::OperatorToMatlabFile("Ainv.dat", AInvOp));

  comm.Barrier();

  Epetra_CrsMatrix * AInv; 
  EPETRA_CHK_ERR(EpetraExt::MatlabFileToCrsMatrix("Ainv.dat", comm, AInv));

  EPETRA_CHK_ERR(AInv->Apply(b,x));

  EPETRA_CHK_ERR(x.Update(1.0, xx, -1.0));
  double residual = 0.0;
  EPETRA_CHK_ERR(x.Norm2(&residual));
  if (verbose) cout << "Norm of difference between computed x and exact x = " << residual << endl;
  int ierr = checkValues(residual,0.0,"Norm of difference between computed A1x1 and A1x1 from file", verbose);

  
  delete A;
  delete AInv;


  #ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(ierr);
}
