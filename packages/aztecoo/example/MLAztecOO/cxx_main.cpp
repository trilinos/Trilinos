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

// Driver to test using ML as a preconditioner to aztecoo.
// Adapted from example in ../AztecOO.
// Reads in a data file from the Matrix Market web site.
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "ml_include.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "ml_epetra_operator.h"
#include "AztecOO.h"
#include "ml_epetra_utils.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }
#define perror1(str,ierr) { fprintf(stderr,"%s %d\n",str,ierr);   exit(-1); }
#define double_quote '"'

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
  
  double norm, d ;

  int has_global_indices, option;
  int i, j, m, mp;
  int ione = 1;

  double time ;
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  printf("proc %d of %d is alive\n",
      comm.MyPID(),comm.NumProc());

  if(argc != 2) perror("error: enter name of data file on command line") ; 
  /* Set exact solution to NULL */
  xexact = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

    Trilinos_Util_read_hb(argv[1], comm.MyPID(), &numGlobalEquations, &n_nonzeros,
             &val,  &bindx, &xguess, &b, &xexact);

  Trilinos_Util_distrib_msr_matrix(comm, &numGlobalEquations, &n_nonzeros, &N_update,
		  &update, &val, &bindx, &xguess, &b, &xexact);

#ifdef DEBUG
  for (i = 0; i<N_update; i++)
    if (val[i] == 0.0 ) printf("Zero diagonal at row %d\n",i);
#endif
  numLocalEquations = N_update;

  /* Make numNzBlks - number of block entries in each block row */

  numNz = new int[numLocalEquations];
  for (i=0; i<numLocalEquations; i++) numNz[i] = bindx[i+1] - bindx[i] + 1;

  /* Make ColInds - Exactly bindx, offset by diag (just copy pointer) */
  ColInds = bindx+numLocalEquations+1;

  Epetra_Map map(numGlobalEquations, numLocalEquations, 
			update, 0, comm);
 
  cout << "Building Epetra_CrsMatrix" << endl;

  Epetra_CrsMatrix A(Copy, map, numNz);
  
  /* Add  rows one-at-a-time */

  for (row=0; row<numLocalEquations; row++)
    {
      row_vals = val + bindx[row];
      col_inds = bindx + bindx[row];
      numEntries = bindx[row+1] - bindx[row];
     if ((ierr = A.InsertGlobalValues(update[row], numEntries, row_vals, col_inds)))
       {
         printf("Row %d:",update[row]);
         perror1("Error putting row:",ierr);
       }
     if ((ierr=(A.InsertGlobalValues(update[row], 1, val+row, update+row)!=0)))
       perror1("Error putting  diagonal",ierr);
    }
  
  if ((ierr=A.FillComplete()))    
    perror1("Error in Epetra_CrsMatrix_FillComplete",ierr);
  
  Epetra_Vector xx(Copy, map, xexact);

  Epetra_Vector bb(Copy, map, b);


  // Make copy of matrix in case it gets scaled by Aztec

  //Epetra_CrsMatrix A_copy(A);

  // Construct a Petra Linear Problem

  Epetra_Vector x(map);

  cout << "Building Epetra_LinearProblem" << endl;

  Epetra_LinearProblem problem(&A, &x, &bb);
  // Construct a solver object for this problem

  cout << "Building AztecOO solver" << endl;

  AztecOO solver(problem);


  // Assert symmetric
  // problem->AssertSymmetric();

  // Set Problem Difficulty Level
  //problem->SetPDL(easy);

  //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  //solver.SetAztecOption(AZ_solver, AZ_gmres);
  //solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  //the following is set within solver.SetPrecOperator()
  //solver.SetAztecOption(AZ_precond, AZ_user_precond);
  solver.SetAztecOption(AZ_solver, AZ_cg);

  // Create and set an ML multilevel preconditioner
  cout << "Creating multigrid hierarchy" << endl;
  ML *ml_handle;
  int N_levels = 10;

  ML_Set_PrintLevel(3);

  ML_Create(&ml_handle,N_levels);
  EpetraMatrix2MLMatrix(ml_handle, 0, &A);

  ML_Aggregate *agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxCoarseSize(agg_object,1);
  N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
                                                  ML_INCREASING, agg_object);
  // Set a symmetric Gauss-Seidel smoother for the MG method
  ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS,
                                  ML_BOTH, 1, ML_DEFAULT);
  ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);
  cout << "Creating Epetra_ML_Operator wrapper" << endl;
  Epetra_ML_Operator  MLop(ml_handle,comm,map,map);
  //note: SetPrecOperator() also sets AZ_user_precond
  solver.SetPrecOperator(&MLop);
 
  //solver.SetAztecOption(AZ_precond, AZ_none);
  //solver.SetAztecOption(AZ_scaling, 8);
  //solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
  //solver.SetAztecOption(AZ_output, 0);
  //solver.SetAztecOption(AZ_graph_fill, 2);
  //solver.SetAztecOption(AZ_overlap, 1);
  //solver.SetAztecOption(AZ_poly_ord, 9);
  //solver.SetAztecParam(AZ_ilut_fill, 1.0);
  //solver.SetAztecParam(AZ_drop, 0.0);
  double rthresh = 1.4;
  //cout << "Rel threshold = " << rthresh << endl;
  solver.SetAztecParam(AZ_rthresh, rthresh);
  double athresh = 10.0;
  //cout << "Abs threshold = " << athresh << endl;
  solver.SetAztecParam(AZ_athresh, athresh);
  solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

  //solver.SetAztecOption(AZ_reorder, 2);

  int Niters = 500;
  solver.SetAztecOption(AZ_kspace, 160);
   
  double norminf = A.NormInf();
  double normone = A.NormOne();
  if (comm.MyPID()==0) 
    cout << "\n Inf-norm of A before scaling = " << norminf 
	 << "\n One-norm of A before scaling = " << normone<< endl << endl;
  solver.Iterate(Niters, 1e-12);
  norminf = A.NormInf();
  normone = A.NormOne(); 
  if (comm.MyPID()==0) 
    cout << "\n Inf-norm of A after  scaling = " << norminf  
	 << "\n One-norm of A after  scaling = " << normone << endl << endl;

  Epetra_Vector bcomp(map);
  assert(A.Multiply(false, x, bcomp)==0);
 
  Epetra_Vector resid(map);
 
  assert(resid.Update(1.0, bb, -1.0, bcomp, 0.0)==0);

  double residual;
  assert(resid.Norm2(&residual)==0);
  if (comm.MyPID()==0) cout << "Residual    = " << residual << endl;

  assert(resid.Update(1.0, xx, -1.0, x, 0.0)==0);

  assert(resid.Norm2(&residual)==0);
  if (comm.MyPID()==0)
    cout << "2-norm of difference between computed and exact solution  = " << residual << endl;

  if (residual>1.0e-5) {
    cout << "Difference between computed and exact solution is large..." << endl      << "Computing norm of A times this difference.  If this norm is small, then matrix is singular"
      << endl;
    assert(A.Multiply(false, resid, bcomp)==0);
    assert(bcomp.Norm2(&residual)==0);
  if (comm.MyPID()==0)
    cout << "2-norm of A times difference between computed and exact solution  = " << residual << endl;
  
  }
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
