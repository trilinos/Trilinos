
/*@HEADER
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
