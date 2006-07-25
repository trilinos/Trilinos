/*@HEADER
// ***********************************************************************
// 
//     AztecOO: An Object-Oriented Aztec Linear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include "AztecOO_ConfigDefs.h"
#ifdef AZTECOO_MPI
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

#ifdef AZTECOO_MPI
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

  Epetra_LinearProblem problem(&A, &x, &bb);
  // Construct a solver object for this problem

  cout << "Building AztecOO solver" << endl;

  AztecOO solver(problem);
  AztecOO_Operator AinvOp(&solver, A->NumGlobalEquations());

  EPETRA_CHK_ERR(EpetraExt::OperatorToMatlabFile("Ainv.dat", AinvOp));

  Epetra_CrsMatrix * Ainv; 
  EPETRA_CHK_ERR(EpetraExt::MatlabFileToCrsMatrix("Ainv.dat", comm, Ainv));

  Ainv->Apply(b,x);

  x.Update(1.0, xx, -1.0);
  double residual = 0.0;
  x.Norm2(&residual);
  if (verbose) cout << "Norm of difference between computed x and exact x = " << residual << endl;
  ierr += checkValues(residual,0.0,"Norm of difference between computed A1x1 and A1x1 from file", verbose);

  
  delete A;
  delete Ainv;


  #ifdef AZTECOO_MPI
  MPI_Finalize() ;
#endif

  return(ierr);
}
