#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Condest.h"
#include "Ifpack_CondestType.h"
#include "Ifpack_Preconditioner.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#ifdef HAVE_IFPACK_AZTECOO
#include "AztecOO.h"
#endif

double Ifpack_Condest(const Ifpack_Preconditioner& IFP,
		      const Ifpack_CondestType CT,
		      const int MaxIters,
		      const double Tol,
		      Epetra_RowMatrix* Matrix)
{
  double ConditionNumberEstimate = -1.0;

  if (CT == Ifpack_Cheap) {

    // Create a vector with all values equal to one
    Epetra_Vector Ones(IFP.OperatorDomainMap());
    Ones.PutScalar(1.0);
    // Create the vector of results
    Epetra_Vector OnesResult(IFP.OperatorRangeMap());
    // Compute the effect of the solve on the vector of ones
    IFPACK_CHK_ERR(IFP.ApplyInverse(Ones, OnesResult)); 
    // Make all values non-negative
    IFPACK_CHK_ERR(OnesResult.Abs(OnesResult)); 
    // Get the maximum value across all processors
    IFPACK_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate)); 

  }
  else if (CT == Ifpack_CG) {

#ifdef HAVE_IFPACK_AZTECOO
    if (Matrix == 0)
      Matrix = (Epetra_RowMatrix*)&(IFP.Matrix());

    Epetra_Vector LHS(IFP.OperatorDomainMap());
    LHS.PutScalar(0.0);
    Epetra_Vector RHS(IFP.OperatorRangeMap());
    RHS.Random();
    Epetra_LinearProblem Problem;
    Problem.SetOperator(Matrix);
    Problem.SetLHS(&LHS);
    Problem.SetRHS(&RHS);

    AztecOO Solver(Problem);
    Solver.SetAztecOption(AZ_output,AZ_none);
    Solver.SetAztecOption(AZ_solver,AZ_cg_condnum);
    Solver.Iterate(MaxIters,Tol);

    const double* status = Solver.GetAztecStatus();
    ConditionNumberEstimate = status[AZ_condnum];
#endif

  } else if (CT == Ifpack_GMRES) {

#ifdef HAVE_IFPACK_AZTECOO
    if (Matrix == 0)
      Matrix = (Epetra_RowMatrix*)&(IFP.Matrix());

    Epetra_Vector LHS(IFP.OperatorDomainMap());
    LHS.PutScalar(0.0);
    Epetra_Vector RHS(IFP.OperatorRangeMap());
    RHS.Random();
    Epetra_LinearProblem Problem;
    Problem.SetOperator(Matrix);
    Problem.SetLHS(&LHS);
    Problem.SetRHS(&RHS);

    AztecOO Solver(Problem);
    Solver.SetAztecOption(AZ_solver,AZ_gmres_condnum);
    Solver.SetAztecOption(AZ_output,AZ_none);
    // the following can be problematic for large problems,
    // but any restart would destroy useful information about
    // the condition number.
    Solver.SetAztecOption(AZ_kspace,MaxIters);
    Solver.Iterate(MaxIters,Tol);

    const double* status = Solver.GetAztecStatus();
    ConditionNumberEstimate = status[AZ_condnum];
#endif
  }

  return(ConditionNumberEstimate);

}
