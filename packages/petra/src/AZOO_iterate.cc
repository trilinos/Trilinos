#include "AZOO_iterate.h"


void AZOO_iterate(double * xsolve, double * b, 
		  int * options, double * params, 
		  double * status, int *proc_config,
		  AZ_MATRIX * Amat,
		  AZ_PRECOND *precond, struct AZ_SCALING *scaling) {

  bool verbose = (options[AZ_output]!=AZ_none); // Print info unless all output is turned off

  Petra_Comm * comm;
  Petra_BlockMap * map;
  Petra_RDP_RowMatrix * A;
  Petra_RDP_Vector * px;
  Petra_RDP_Vector * pb;

  int ierr = Aztec2Petra(proc_config, Amat, xsolve, b, comm, map, A, px, pb);


  Petra_RDP_LinearProblem problem(A, px, pb);

  Petra_RDP_Vector * leftScaleVec = 0;
  Petra_RDP_Vector * rightScaleVec = 0;
  bool doRowScaling = false;
  bool doColScaling = false;
  
  if ((options[AZ_scaling]==AZ_Jacobi) || options[AZ_scaling]==AZ_BJacobi) {
    doRowScaling = true;
    leftScaleVec = new Petra_RDP_Vector(*map);
    A->ExtractDiagonalCopy(*leftScaleVec); // Extract diagonal of matrix
    leftScaleVec->Reciprocal(*leftScaleVec); // invert it
  }

  else if (options[AZ_scaling]==AZ_row_sum) {
    doRowScaling = true;
    leftScaleVec = new Petra_RDP_Vector(*map);
    A->InvRowSums(*leftScaleVec);
  }
  else if (options[AZ_scaling]==AZ_sym_diag) {
    doRowScaling = true;
    doColScaling = true;
    leftScaleVec = new Petra_RDP_Vector(*map);
    A->ExtractDiagonalCopy(*leftScaleVec); // Extract diagonal of matrix

    int length = leftScaleVec->MyLength();
    for (int i=0; i<length; i++) (*leftScaleVec)[i] = sqrt(fabs((*leftScaleVec)[i])); // Take its sqrt

    rightScaleVec = leftScaleVec; // symmetric, so left and right the same
    leftScaleVec->Reciprocal(*leftScaleVec); // invert it
  }
  else if (options[AZ_scaling]==AZ_sym_row_sum) {
    doRowScaling = true;
    doColScaling = true;
    leftScaleVec = new Petra_RDP_Vector(*map);
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



  Aztec_OO solver(problem);

  solver.SetAllAztecParams(params); // set all Aztec_OO params with user-provided params
  solver.SetAllAztecOptions(options); // set all Aztec_OO options with user-provided options

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
  
  return;
}
