//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to provide for residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface.H"

#include "GenericProblem.H"

class  Problem_Interface : public NOX::Epetra::Interface
{
public:
  Problem_Interface(GenericProblem& Problem) :
  problem(Problem)
  {};

  ~Problem_Interface()
  {};

  //! Compute and return F.  Returns true if computation was successful.
  bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag = F)
  {
    return problem.evaluate(F_ONLY, &x, &FVec, NULL);
  }

  //! Compute an explicit Jacobian.  Returns true if computation was successful.
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
  {
    Epetra_RowMatrix* Jacobian = dynamic_cast<Epetra_RowMatrix*>(&Jac);
    if (Jacobian == NULL) {
      cout << "ERROR: Problem_Interface::computeJacobian() - The supplied"
  	 << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
      throw;
    }
    return problem.evaluate(MATRIX_ONLY, &x, NULL, Jacobian);
  }
  
  //! Compute the Epetra_RowMatrix M, that will be used by the Aztec preconditioner instead of the Jacobian.  This is used when there is no explicit Jacobian present (i.e. Matrix-Free Newton-Krylov).  This MUST BE an Epetra_RowMatrix since the Aztec preconditioners need to know the sparsity pattern of the matrix.  Returns true if computation was successful.
  bool computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
  {
    return problem.computePrecMatrix(x, M);
  }

  //! Computes a user supplied preconditioner based on input vector x.  Returns true if computation was successful.
  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& precOp)
  {
    return problem.computePreconditioner(x, precOp);
  };

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 
  GenericProblem& problem;
};

#endif
