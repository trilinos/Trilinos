#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsKundertSparse.h"

// Local prototypes
double computeNorm(const Epetra_Vector & x);
double computeDiffNorm(const Epetra_Vector & x, const Epetra_Vector & y);
void perturbDiag(Epetra_CrsMatrix *A, const Epetra_Vector & diagA, double scaleFactor);
void compareSolutions(const Epetra_CrsMatrix & A, const Epetra_Vector & x, const Epetra_Vector & xx);

int main(int argc, char *argv[])
{
  int    *update;                  /* vector elements updated on this node. */
  int    *indx;   /* MSR format of real and imag parts */
  int    *bindx;
  int    *bpntr;
  int    *rpntr;
  int    *cpntr;
  int    indexBase = 0; 
  double *val;
  double *xguess, *b, *xexact, *xsolve;
  int    n_nonzeros, n_blk_nonzeros, ierr;
  int    N_update;           /* # of block unknowns updated on this node    */
  int    numLocalEquations;
                                 /* Number scalar equations on this node */
  int    numGlobalEquations, numGlobalBlocks; /* Total number of equations */
  int    numLocalBlocks;
  int    *blockSizes, *numNzBlks, *blkColInds;
  int    *numNz, *ColInds;
  int    N_external, N_blk_eqns;
  int    blk_row, *blk_col_inds;
  int    row,     *col_inds, numEntries;
  double *row_vals;

  double *val_msr;
  int *bindx_msr;
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  cout << comm << endl;
  

  if(argc != 2) perror("error: enter name of data file on command line") ; 
  /* Set exact solution to NULL */
  xexact = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

    Trilinos_Util_read_hb(argv[1], comm.MyPID(), &numGlobalEquations, &n_nonzeros,
             &val,  &bindx, &xguess, &b, &xexact);

  Trilinos_Util_distrib_msr_matrix(comm, &numGlobalEquations, &n_nonzeros, &N_update,
		  &update, &val, &bindx, &xguess, &b, &xexact);

  numLocalEquations = N_update;

  /* Make numNzBlks - number of block entries in each block row */

  numNz = new int[numLocalEquations];
  for (int i=0; i<numLocalEquations; i++) numNz[i] = bindx[i+1] - bindx[i] + 1;

  /* Make ColInds - Exactly bindx, offset by diag (just copy pointer) */
  ColInds = bindx+numLocalEquations+1;

  Epetra_Map map(numGlobalEquations, numLocalEquations, update, 0, comm);
 
  Epetra_CrsMatrix A(Copy, map, numNz);
  
  /* Add  rows one-at-a-time */

  for (row=0; row<numLocalEquations; row++) {
    row_vals = val + bindx[row];
    col_inds = bindx + bindx[row];
    numEntries = bindx[row+1] - bindx[row];
    assert(A.InsertGlobalValues(update[row], numEntries, row_vals, col_inds)>=0);
    assert(A.InsertGlobalValues(update[row], 1, val+row, update+row)>=0);
  }  
  assert(A.TransformToLocal()==0);

  Epetra_Vector xx(Copy, map, xexact);

  Epetra_Vector bb(Copy, map, b);


  // Construct a Petra Linear Problem

  Epetra_Vector x(map);
  Epetra_LinearProblem problem(&A, &x, &bb);

  Epetra_Vector diagA(map);

  A.ExtractDiagonalCopy(diagA);

  // Call Kundert's Sparse solver to Solve
  Epetra_Time timer(comm);

  cout << endl << endl << "*************************************************" << endl;
  cout << "One Norm of initial matrix = " << A.NormOne() << endl;
  cout << "Initializing sparse matrix storage" << endl;
  double start = timer.ElapsedTime();
  Epetra_CrsKundertSparse solver(&problem);
  double spOrderAndFactorTime = timer.ElapsedTime() - start;

  cout << "Computing ordering, factoring and solving first system" << endl;
  start = timer.ElapsedTime();
  solver.Solve();
  double spSolveTime = timer.ElapsedTime() - start;

  compareSolutions(A, x, xx);

  cout << endl << endl << "*************************************************" << endl;
  cout << "Factoring and solving using new matrix values but old ordering" << endl;

  perturbDiag(&A, diagA, 1.1); // Perturb diagonal a bit to make matrix different

  A.Multiply(false, xx, bb); // Modify RHS to match perturbed matrix

  cout << "One Norm of second matrix = " << A.NormOne() << endl;
  start = timer.ElapsedTime();
  solver.Solve();
  double spFactorAndSolveTime1 = timer.ElapsedTime() - start;

  compareSolutions(A, x, xx);

  cout << endl << endl << "*************************************************" << endl;
  cout << "Doing the same thing again, should have similar timings" << endl;
  perturbDiag(&A, diagA, 1.2); // Perturb diagonal a bit to make matrix different

  A.Multiply(false, xx, bb); // Modify RHS to match perturbed matrix
  cout << "One Norm of third matrix = " << A.NormOne() << endl;
  start = timer.ElapsedTime();
  solver.Solve();
  double spFactorAndSolveTime2 = timer.ElapsedTime() - start;

  compareSolutions(A, x, xx);

  cout << endl
       << "*************************************************" << endl
       << "Time for matrix initialization                 = " << spOrderAndFactorTime << endl
       << "Time for first ordering/factoring/solving      = " << spSolveTime << endl
       << "Total time for first solve                     = " << spOrderAndFactorTime+spSolveTime 
       << endl <<endl
       << "Time for second factor/solve (no ordering)     = " << spFactorAndSolveTime1 << endl <<endl
       << "Time for third factor/solve (identical to 2nd) = " << spFactorAndSolveTime1 << endl <<endl
       << "*************************************************" << endl;
  free ((void *) xguess);
  free ((void *) b);
  free ((void *) xexact);
  free ((void *) val);
  free ((void *) bindx);
  free ((void *) update);

  delete [] numNz;

				       
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}

double computeNorm(const Epetra_Vector & x) {

  double norm;
  assert(x.Norm2(&norm)==0);
  return(norm);
}
  
double computeDiffNorm(const Epetra_Vector & x, const Epetra_Vector & y) {

  double normdiff, normx;
  assert(x.Norm2(&normx)==0);
  Epetra_Vector resid(x.Map()); 
  assert(resid.Update(1.0, x, -1.0, y, 0.0)==0);
  assert(resid.Norm2(&normdiff)==0);
  return(normdiff/normx);
}
  
void perturbDiag(Epetra_CrsMatrix *A, const Epetra_Vector & diagA, double scaleFactor) {

  int NumMyRows = A->NumMyRows();
  for (int i=0; i< NumMyRows; i++) {
    double newDiagValue = diagA[i]*scaleFactor;
    A->ReplaceMyValues(i, 1, &newDiagValue, &i);
  }
}

void compareSolutions(const Epetra_CrsMatrix & A, const Epetra_Vector & x, const Epetra_Vector & xx) {
  Epetra_Vector Ax(x);
  assert(A.Multiply(false, x, Ax)==0);
  Epetra_Vector Axx(x);
  assert(A.Multiply(false, xx, Axx)==0);

  double normdiff = computeDiffNorm(x,xx);  
  cout << "Norm of (xexact - xcomp)/Norm xcomp      = " << normdiff << endl;
  double normAdiff = computeDiffNorm(Ax,Axx);
  cout << "Norm of (A*xexact - A*xcomp)/Norm Axcomp = " << normAdiff << endl;

  if (normdiff>1e-4 && normAdiff < 1e-4) 
    cout << endl 
	 << "This matrix appears to be singular since two different solutions give a zero residual."
	 << endl;
}
