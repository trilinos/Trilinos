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
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
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
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }
#define perror1(str,ierr) { fprintf(stderr,"%s %d\n",str,ierr);   exit(-1); }
#define double_quote '"'

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  cout << comm << endl;

  // if (comm.MyPID()==0) {

  //    cout << "Press a key to continue..." << endl;
  //  char aaa;
  //  cin >> aaa;
  //}
  
  comm.Barrier();
  

  if(argc != 3) {
    cerr << "Error: enter name of data and partition file on command line" << endl;
    abort();
  }

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

  int NumGlobalEquations, NumMyNonzeros;
  double *val_msr, *xguess, *b, *xexact;
  int *bindx_msr;

  /* Set exact solution to NULL */
  xexact = NULL;
  Trilinos_Util_read_hb(argv[1], comm.MyPID(), &NumGlobalEquations, &NumMyNonzeros,
			&val_msr,  &bindx_msr, &xguess, &b, &xexact);

  double *val;
  int NumGlobalElements, *indx, *rpntr, *cpntr, *bpntr, *bindx;
  int NumMyBlockEntries, NumMyElements, * MyGlobalElements;
  
  Trilinos_Util_create_vbr(comm, argv[2],
			   &NumGlobalEquations, &NumGlobalElements, 
			   &NumMyNonzeros, &NumMyBlockEntries,
			   &NumMyElements, &MyGlobalElements,
			   bindx_msr, val_msr,
			   &val, &indx, &rpntr, &cpntr,
			   &bpntr, &bindx);
    
    if(comm.MyPID()==0)
      {
	free ((void *) val_msr);
	free ((void *) bindx_msr);
	free ((void *) cpntr);
      }
  
  Trilinos_Util_distrib_vbr_matrix(comm, &NumGlobalEquations, &NumGlobalElements,
				   &NumMyNonzeros, &NumMyBlockEntries, &NumMyElements,
				   &MyGlobalElements, &val, 
				   &indx, &rpntr, &cpntr, &bpntr,
				   &bindx, &xguess, &b, &xexact);
  

  /* Make numNzBlks - number of block entries in each block row */

  int * ElementSizeList = new int[NumMyElements];
  
  for (int i=0; i<NumMyElements; i++) ElementSizeList[i] = rpntr[i+1] - rpntr[i];

  Epetra_BlockMap map(NumGlobalElements, NumMyElements, MyGlobalElements, 
		     ElementSizeList, 0, comm);
 
  Epetra_VbrMatrix A(Copy, map, 0);
  
  /* Add block rows one-at-a-time */

  for (int i=0; i<NumMyElements; i++) {
    int BlockRow = MyGlobalElements[i];
    int NumBlockEntries = bpntr[i+1] - bpntr[i];
    int *BlockIndices = bindx + bpntr[i];
    int ierr = A.BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices);
    if (ierr!=0) {
      cerr << "Error in BeginInsertGlobalValues(GlobalBlockRow = " << BlockRow 
	   << ") = " << ierr << endl; 
      abort();
    }
    int LDA = ElementSizeList[i];
    int NumRows = LDA;
    for (int j=bpntr[i]; j<bpntr[i+1]; j++) {
      int NumCols = (indx[j+1] - indx[j])/LDA;
      double * Values = val + indx[j];
      ierr = A.SubmitBlockEntry(Values, LDA, NumRows, NumCols);
      if (ierr!=0) {
	cerr << "Error in SubmitBlockEntry, GlobalBlockRow = " << BlockRow 
	     << "GlobalBlockCol = " << BlockIndices[j] << "Error = " << ierr << endl; 
	abort();
      }
    }
    ierr = A.EndSubmitEntries();
    if (ierr!=0) {
      cerr << "Error in EndSubmitEntries(GlobalBlockRow = " << BlockRow 
	   << ") = " << ierr << endl; 
      abort();
    }
  }  
  int ierr=A.FillComplete();    
  if (ierr!=0) perror1("Error in Epetra_VbrMatrix FillComplete",ierr);
  
  //cout << A<< endl;
  double * xexactt = xexact;
  Epetra_Vector xx(Copy, map, xexactt);

  double * bt = b;
  Epetra_Vector bb(Copy, map, bt);


  // Make copy of matrix in case it gets scaled by Aztec

  //Epetra_CrsMatrix A_copy(A);

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
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  //solver.SetAztecOption(AZ_precond, AZ_ls);
  //solver.SetAztecOption(AZ_scaling, 8);
  //solver.SetAztecOption(AZ_subdomain_solve, AZ_bilu_ifp); 
  //solver.SetAztecOption(AZ_output, 0);
  solver.SetAztecOption(AZ_graph_fill, 1);
  //solver.SetAztecOption(AZ_overlap, 1);
  //solver.SetAztecOption(AZ_poly_ord, 1);
  //solver.SetAztecParam(AZ_ilut_fill, 2.0);
  //solver.SetAztecParam(AZ_drop, 0.1);
  //solver.SetAztecParam(AZ_rthresh, 1e-3);
  //solver.SetAztecParam(AZ_athresh, 1e-3);

  //solver.SetAztecOption(AZ_reorder, 2);

  int Niters = 500;
  solver.SetAztecOption(AZ_kspace, Niters);
  double norminf = A.NormInf();
  double normone = A.NormOne();
  if (comm.MyPID()==0) 
    cout << "\n Inf-norm of A before scaling = " << norminf 
	 << "\n One-norm of A before scaling = " << normone<< endl << endl;
  bool doRowScaling = false;
  bool doColScaling = false;
  Epetra_Vector rowsum(map);
  A.InvRowSums(rowsum);
  if (doRowScaling) problem.LeftScale(rowsum);
  Epetra_Vector colsum(map);
  A.InvColSums(colsum);
  if (doColScaling) problem.RightScale(colsum);

  if (doRowScaling || doColScaling) {
    norminf = A.NormInf();
    normone = A.NormOne();
    if (comm.MyPID()==0) 
      cout << "\n Inf-norm of A after  scaling = " << norminf  
	   << "\n One-norm of A after  scaling = " << normone << endl << endl;
  } 

  solver.CheckInput();
  solver.Iterate(Niters, 5.0e-14);
  
  if (doRowScaling) {
    Epetra_Vector invrowsum(map);
    invrowsum.Reciprocal(rowsum);
    problem.LeftScale(invrowsum);
  }
  if (doColScaling) {
    Epetra_Vector invcolsum(map);
    invcolsum.Reciprocal(colsum);
    problem.RightScale(invcolsum);
  }
  
  Epetra_Vector bcomp(map);
  if ((ierr=A.Multiply(false, x, bcomp)))
    perror1("Error in matvec",ierr);
 
  Epetra_Vector resid(map); 
 
  if ((ierr=resid.Update(1.0, bb, -1.0, bcomp, 0.0)))  
    perror1("Error in linComb",ierr);

  double residual;
  if ((ierr=resid.Norm2(&residual)))
    perror1("Error in Epetra_Vector_putVector",ierr);
  if (comm.MyPID()==0)
      printf("Residual    = %22.16g\n",residual);

  // Unscale solution

  if ((ierr=resid.Update(1.0, xx, -1.0, x, 0.0)))
    perror1("Error in linComb",ierr);

  if ((ierr=resid.Norm2(&residual)))
    perror1("Error in Epetra_Vector_putVector",ierr);
  if (comm.MyPID()==0)
      printf("2-norm of difference between computed and exact solution  = %22.16g\n",residual);

  free ((void *) xguess);
  free ((void *) b);
  free ((void *) xexact);
  free ((void *) val);
  free ((void *) bindx);
  free ((void *) bpntr);
  free ((void *) rpntr);
  free ((void *) indx);
  free ((void *) MyGlobalElements);

  delete [] ElementSizeList;
				       
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
