//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NLS_PetraGroup to provide for 
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface.H"

// ---------- Forward Declarations ----------
class DennisSchnabel;

class  Problem_Interface : public NOX::Epetra::Interface
{
public:
  Problem_Interface(DennisSchnabel& Problem);
  ~Problem_Interface();

  //! Compute and return RHS.  Returns true if computation was successful.
  bool computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS);

  //! Compute an explicit Jacobian.  Returns true if computation was successful.
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);
  
  //! Compute the Epetra_RowMatrix M, that will be used by the Aztec preconditioner instead of the Jacobian.  This is used when there is no explicit Jacobian present (i.e. Matrix-Free Newton-Krylov).  This MUST BE an Epetra_RowMatrix since the Aztec preconditioners need to know the sparsity pattern of the matrix.  Returns true if computation was successful.
  bool computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M);

  //! Computes a user supplied preconditioner based on input vector x.  Returns true if computation was successful.
  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M);

  //! Set a parameter in the user's code.
  bool setParameter(string param, double value);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 
  DennisSchnabel& problem;
};

#endif

