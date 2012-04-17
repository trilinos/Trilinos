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
