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
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "PetraSuperLU.h"

// Local prototypes
double computeNorm(const Epetra_Vector & x);
double computeDiffNorm(const Epetra_Vector & x, const Epetra_Vector & y);

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

  // Call SuperLU to Solve

  PetraSuperLU(&problem);

  // Test result that was computed

  bool verbose = (comm.MyPID()==0);

  if (verbose) cout << endl << "*************************************************" << endl << endl;

  double norm = computeNorm(xx);
  if (verbose) cout << "Norm of file-read solution (xexact)      = " << norm << endl;
  norm = computeNorm(x);
  if (verbose) cout << "Norm of computed  solution (xcomp)       = " << norm << endl;

  Epetra_Vector Ax(map);
  assert(A.Multiply(false, x, Ax)==0);
  Epetra_Vector Axx(map);
  assert(A.Multiply(false, xx, Axx)==0);

  double normdiff = computeDiffNorm(x,xx);  
  if (verbose) cout << "Norm of (xexact - xcomp)/Norm xcomp      = " << normdiff << endl;
  double normAdiff = computeDiffNorm(Ax,Axx);
  if (verbose) cout << "Norm of (A*xexact - A*xcomp)/Norm Axcomp = " << normAdiff << endl;

  if (normdiff>1e-4 && normAdiff < 1e-4) 
    cout << endl 
	 << "This matrix appears to be singular since two different solutions give a zero residual."
	 << endl;


  if (verbose) cout << endl << "*************************************************" << endl;
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
  
