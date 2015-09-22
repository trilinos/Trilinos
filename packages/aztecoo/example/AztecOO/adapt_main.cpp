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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
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
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }
#define perror1(str,ierr) { fprintf(stderr,"%s %d\n",str,ierr);   exit(-1); }
#define double_quote '"'

int main(int argc, char *argv[])
{
  int    *update;                  /* vector elements updated on this node. */
  int    *bindx;
  double *val;
  double *xguess, *b, *xexact;
  int    n_nonzeros;
  int    N_update;           /* # of block unknowns updated on this node    */
  int    numLocalEquations;
                                 /* Number scalar equations on this node */
  int    numGlobalEquations; /* Total number of equations */
  int    *numNz, *ColInds;
  int    row,     *col_inds, numEntries;
  double *row_vals;

  int i;

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  cout << comm << endl;

  if(argc != 2) cerr << "error: enter name of data file on command line" << endl; 

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
  for (i=0; i<numLocalEquations; i++) numNz[i] = bindx[i+1] - bindx[i] + 1;

  /* Make ColInds - Exactly bindx, offset by diag (just copy pointer) */
  ColInds = bindx+numLocalEquations+1;

  Epetra_Map map(numGlobalEquations, numLocalEquations, update, 0, comm);
 
  Epetra_CrsMatrix A(Copy, map, numNz);
  
  /* Add  rows one-at-a-time */

  for (row=0; row<numLocalEquations; row++) {
      row_vals = val + bindx[row];
      col_inds = bindx + bindx[row];
      numEntries = bindx[row+1] - bindx[row];
      assert(A.InsertGlobalValues(update[row], numEntries, row_vals, col_inds)==0);
      assert(A.InsertGlobalValues(update[row], 1, val+row, update+row)==0);
    }
  
  assert(A.FillComplete()==0);   
  
  Epetra_Vector xx(Copy, map, xexact);

  Epetra_Vector bb(Copy, map, b);



  // Construct a Petra Linear Problem

  Epetra_Vector x(map);
  Epetra_LinearProblem problem(&A, &x, &bb);
  // Construct a solver object for this problem
  AztecOO solver(problem);


  // Assert symmetric
  // problem->AssertSymmetric();

  // Set Problem Difficulty Level
  //problem->SetPDL(easy);

  //solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  //solver.SetAztecOption(AZ_precond, AZ_ls);
  //solver.SetAztecOption(AZ_scaling, 8);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  //solver.SetAztecOption(AZ_subdomain_solve, AZ_bilu_ifp);
  bool bilu = false;
  //solver.SetAztecOption(AZ_output, 0);
  //solver.SetAztecOption(AZ_graph_fill, 2);
  solver.SetAztecOption(AZ_overlap, 0);
  //solver.SetAztecOption(AZ_reorder, 0);
  //solver.SetAztecOption(AZ_poly_ord, 9);
  solver.SetAztecParam(AZ_ilut_fill, 1.0);
  solver.SetAztecParam(AZ_drop, 0.0);
  //double rthresh = 1.01;
  //cout << "Rel threshold = " << rthresh << endl;
  //solver.SetAztecParam(AZ_rthresh, rthresh);
  //double athresh = 1.0e-2;
  //cout << "Abs threshold = " << athresh << endl;
  //solver.SetAztecParam(AZ_athresh, athresh);

  //solver.SetAztecOption(AZ_/conv, AZ_noscaled);
  //solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e12);


  int Niters = 400;
  solver.SetAztecOption(AZ_kspace, Niters);
  
  double norminf = A.NormInf();
  double normone = A.NormOne();
  if (comm.MyPID()==0) 
    cout << "\n Inf-norm of A before scaling = " << norminf 
	 << "\n One-norm of A before scaling = " << normone<< endl << endl;


  if (bilu) {
    int NumTrials = 3;
    double athresholds[] = {0.0, 1.0E-14, 1.0E-3};
    double rthresholds[] = {0.0, 1.0E-14, 1.0E-3};
    double condestThreshold = 1.0E16;
    double maxFill = 4.0;
    int maxKspace = 4*Niters;
    solver.SetAdaptiveParams(NumTrials, athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
  }
  else {
    int NumTrials = 7;
    double athresholds[] = {0.0, 1.0E-12, 1.0E-12, 1.0E-5, 1.0E-5, 1.0E-2, 1.0E-2};
    double rthresholds[] = {1.0, 1.0,     1.01,    1.0,    1.01,   1.01,   1.1   };
    double condestThreshold = 1.0E16;
    double maxFill = 4.0;
    int maxKspace = 4*Niters;
    solver.SetAdaptiveParams(NumTrials, athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
  }

  Epetra_Time timer(comm);
  solver.AdaptiveIterate(Niters, 1, 5.0e-7);
  double atime = timer.ElapsedTime();
  if (comm.MyPID()==0) cout << "AdaptiveIterate total time = " << atime << endl;

  
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
    cout << "Difference between computed and exact solution is large..." << endl
	 << "Computing norm of A times this difference.  If this norm is small, then matrix is singular" 
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
