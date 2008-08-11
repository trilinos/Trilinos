/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
