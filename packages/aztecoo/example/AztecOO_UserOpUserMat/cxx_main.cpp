#include "AztecOO.h"
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
  Epetra_Vector xx(A.DomainMap());
  Epetra_Vector x(A.DomainMap());
  Epetra_Vector b(A.RangeMap());

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
 
  Epetra_Vector bcomp(A.RangeMap());
  A.Apply(x, bcomp); 
 
  Epetra_Vector resid(A.RangeMap());
 
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
  delete PrecMatrix;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0;
}
