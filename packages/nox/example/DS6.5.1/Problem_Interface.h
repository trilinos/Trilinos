//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NOX::Petra::Group to provide for 
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#ifdef TFLOP
#include <iostream.h>
#else
#include <iostream>
#endif

#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface.H"

// ---------- Forward Declarations ----------
class Fill;

class  Problem_Interface : public NOX::Epetra::Interface 
{
public:
  Problem_Interface ();
  ~Problem_Interface ();

  //! Compute and return RHS
  void computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS);

  //! Compute RHS
  void computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac);

  //! Register objects from the main code required to call the fill routines
  void registerFill(Fill * tmpFill);

  //! Compute the matrix M that will be used as the preconditioner.
  void computePreconditioner(Epetra_RowMatrix& M){};

  //! Return the action of the preconditioner on a vector.
  void preconditionVector(Epetra_Vector& y){};

  //! Object to store the auxiliary data from user's code required 
  //! for RHS and Jacobian fills.
  Fill* LO;
};

#endif

