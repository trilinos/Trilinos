//-----------------------------------------------------------------------------
#ifndef NLS_Interface_H
#define NLS_Interface_H

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
#include "NLS_PetraGroupInterface.H"

// ---------- Forward Declarations ----------
class Fill;

class  NLS_Interface : public NLS_PetraGroupInterface 
{
public:
  NLS_Interface ();
  ~NLS_Interface ();

  //! Compute and return RHS
  void computeRHS(Epetra_Vector& x, Epetra_Vector& RHS);

  //! Compute RHS
  void computeJacobian(Epetra_Vector& x, Epetra_RowMatrix& Jac);

  //! Register objects from the main code required to call the fill routines
  void registerFill(Fill * tmpFill);

  //! Object to store the auxiliary data from user's code required 
  //! for RHS and Jacobian fills.
  Fill* LO;
};

#endif

