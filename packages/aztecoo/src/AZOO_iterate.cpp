
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "AZOO_iterate.h"


void AZOO_iterate(double * xsolve, double * b, 
		  int * options, double * params, 
		  double * status, int *proc_config,
		  AZ_MATRIX * Amat,
		  AZ_PRECOND *precond, struct AZ_SCALING *scaling) {

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
