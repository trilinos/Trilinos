//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NLS_PetraGroup to provide for 
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

  //! Object to store the auxiliary data from user's code required 
  //! for RHS and Jacobian fills.
  Fill* LO;
};

#endif

