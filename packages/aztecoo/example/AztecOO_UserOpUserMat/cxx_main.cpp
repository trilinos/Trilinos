//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#include "AztecOO.h"
#include "AztecOO_Version.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Poisson2dOperator.h"

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  if (comm.MyPID()==0)
    cout << AztecOO_Version() << endl << endl;

  if (argc!=3) {
    cerr << "Usage: " << argv[0] << " nx ny" << endl;
    exit(1);
  }

  bool vb = (comm.MyPID()==0);

  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  if (vb) cout << "Solving an " << nx << " by " << ny 
               << " grid point Poisson problem. " << endl;

  // Generate nx by ny Poisson operator

  Poisson2dOperator A(nx, ny, comm);

  
  // Generate vectors (xx will be used to generate RHS b)
  Epetra_Vector xx(A.OperatorDomainMap());
  Epetra_Vector x(A.OperatorDomainMap());
  Epetra_Vector b(A.OperatorRangeMap());

  xx.Random();
  A.Apply(xx, b);

  if (vb) cout << "Building Tridiagonal Approximation to Poisson" << endl;
  Epetra_CrsMatrix * PrecMatrix = A.GeneratePrecMatrix();

  if (vb) cout << "Building Epetra_LinearProblem" << endl;
  Epetra_LinearProblem problem(&A, &x, &b);


  if (vb) cout << "Building AztecOO solver" << endl;
  AztecOO solver(problem);

  if (vb) cout << "Setting Preconditioner Matrix" << endl;
  solver.SetPrecMatrix(PrecMatrix);
  //solver.SetAztecOption(AZ_precond, AZ_none);

  int Niters = nx*ny;
  solver.SetAztecOption(AZ_kspace, Niters);

  solver.Iterate(Niters, 1.0E-12);
 
  Epetra_Vector bcomp(A.OperatorRangeMap());
  A.Apply(x, bcomp); 
 
  Epetra_Vector resid(A.OperatorRangeMap());
 
  resid.Update(1.0, b, -1.0, bcomp, 0.0);

  double residual;
  resid.Norm2(&residual);
  if (vb) cout << "Residual    = " << residual << endl;

  resid.Update(1.0, xx, -1.0, x, 0.0);

  resid.Norm2(&residual);
  if (vb)
    cout << "2-norm of difference between computed and exact solution  = " 
         << residual << endl;

  if (residual>1.0e-5) {
    if (vb) cout << "Difference between computed and exact solution is large." 
         << endl      
         << "Computing norm of A times this difference.  "
         << endl      
         << "If this norm is small, then matrix is singular"
         << endl;
    A.Apply(resid, bcomp);
    bcomp.Norm2(&residual);
  if (vb) cout << "2-norm of A times difference between computed and exact "
               << "solution  = " << residual << endl;
  
  }
  else {
   if (vb) cout << "Solver converged" << endl;
  }
  delete PrecMatrix;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0;
}
