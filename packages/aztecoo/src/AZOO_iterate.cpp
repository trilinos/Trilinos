
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

#include "AZOO_iterate.h"
#include <stdlib.h>


void AZOO_iterate(double * xsolve, double * b, 
		  int * options, double * params, 
		  double * status, int *proc_config,
		  AZ_MATRIX * Amat,
		  AZ_PRECOND *precond, struct AZ_SCALING *scaling)
{
  (void)precond;
  (void)scaling;
  bool verbose = (options[AZ_output]!=AZ_none); // Print info unless all output is turned off

  Epetra_Comm * comm;
  Epetra_BlockMap * map;
  Epetra_RowMatrix * A;
  Epetra_Vector * px;
  Epetra_Vector * pb;
  int * global_indices;

  int ierr = Aztec2Petra(proc_config, Amat, xsolve, b, comm, map, A, px, pb, &global_indices);
  if (ierr!=0) {
    cerr << "Error detected in Aztec2Petra. Value = " << ierr << endl;
    exit(1);
  }


  Epetra_LinearProblem problem(A, px, pb);

  Epetra_Vector * leftScaleVec = 0;
  Epetra_Vector * rightScaleVec = 0;
  bool doRowScaling = false;
  bool doColScaling = false;
  
  if ((options[AZ_scaling]==AZ_Jacobi) || options[AZ_scaling]==AZ_BJacobi) {
    doRowScaling = true;
    leftScaleVec = new Epetra_Vector(*map);
    A->ExtractDiagonalCopy(*leftScaleVec); // Extract diagonal of matrix
    leftScaleVec->Reciprocal(*leftScaleVec); // invert it
  }

  else if (options[AZ_scaling]==AZ_row_sum) {
    doRowScaling = true;
    leftScaleVec = new Epetra_Vector(*map);
    A->InvRowSums(*leftScaleVec);
  }
  else if (options[AZ_scaling]==AZ_sym_diag) {
    doRowScaling = true;
    doColScaling = true;
    leftScaleVec = new Epetra_Vector(*map);
    A->ExtractDiagonalCopy(*leftScaleVec); // Extract diagonal of matrix

    int length = leftScaleVec->MyLength();
    for (int i=0; i<length; i++) (*leftScaleVec)[i] = sqrt(fabs((*leftScaleVec)[i])); // Take its sqrt

    rightScaleVec = leftScaleVec; // symmetric, so left and right the same
    leftScaleVec->Reciprocal(*leftScaleVec); // invert it
  }
  else if (options[AZ_scaling]==AZ_sym_row_sum) {
    doRowScaling = true;
    doColScaling = true;
    leftScaleVec = new Epetra_Vector(*map);
    A->InvRowSums(*leftScaleVec);
    int length = leftScaleVec->MyLength();
    for (int i=0; i<length; i++) (*leftScaleVec)[i] = sqrt(fabs((*leftScaleVec)[i])); // Take its sqrt

    rightScaleVec = leftScaleVec; // symmetric, so left and right the same
  }
  if ((doRowScaling || doColScaling) && verbose) {
    double norminf = A->NormInf();
    double normone = A->NormOne();
    if (comm->MyPID()==0) 
      cout << "\n Inf-norm of A before scaling = " << norminf 
	   << "\n One-norm of A before scaling = " << normone<< endl << endl;
  }
  if (doRowScaling) problem.LeftScale(*leftScaleVec);
  if (doColScaling) problem.RightScale(*rightScaleVec);

  if ((doRowScaling || doColScaling) && verbose) {
    double norminf = A->NormInf();
    double normone = A->NormOne();
    if (comm->MyPID()==0) 
      cout << "\n Inf-norm of A after  scaling = " << norminf  
	   << "\n One-norm of A after  scaling = " << normone << endl << endl;
  }



  AztecOO solver(problem);

  solver.SetAllAztecParams(params); // set all AztecOO params with user-provided params
  solver.SetAllAztecOptions(options); // set all AztecOO options with user-provided options

  solver.CheckInput();
  solver.SetAztecOption(AZ_scaling, AZ_none); // Always must have scaling off
  solver.Iterate(options[AZ_max_iter], params[AZ_tol]);
  solver.GetAllAztecStatus(status);
  
  if (doColScaling) {
    rightScaleVec->Reciprocal(*rightScaleVec);
    problem.RightScale(*rightScaleVec);
  }
  if (doRowScaling) {
    leftScaleVec->Reciprocal(*leftScaleVec);
    problem.LeftScale(*leftScaleVec);
  }

  if ((rightScaleVec!=0) && (rightScaleVec!=leftScaleVec)) delete rightScaleVec;
  if (leftScaleVec!=0) delete leftScaleVec;

  delete pb; // These are all objects created here and we have to delete them
  delete px;
  delete A;
  delete map;
  delete comm;
  if (global_indices!=0) AZ_free((void *) global_indices); // Note: we used a special version of free here

  return;
}
