#ifndef IFPACK_PRECONDITIONER_H
#define IFPACK_PRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Operator.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Teuchos_ParameterList.hpp"
#endif
class Epetra_RowMatrix;

//! Ifpack_Preconditioner: basic class for preconditioning in Ifpack

/*!
  Class Ifpack_Preconditioner is a pure virtual class, and it defines
  the structure of all Ifpack preconditioners.

  This class is a simple extension to Epetra_Operator. It provides 
  the following additional methods:
  - Compute(), that should compute all is required to apply the
    preconditioner;
  - IsComputed(), that should return true is the preconditioner
    has been successfully computed, false otherwise.
  - Condest() returns an estimation of the condition number, or -1.0
    if not available
  - Matrix() returns a reference to the matrix to be preconditioned.

  If IFPACK is configure with Teuchos support, method SetParameters()
  should be adopted. Otherwise, users can set parameters (one at-a-time),
  using methods SetParameter(), for integers and doubles.
  \date Sep-04

  Ifpack_Preconditioner objects overload the << operator. Derived
  classes should specify a Print() method, that will be used in
  operator <<.

*/

class Ifpack_Preconditioner : public Epetra_Operator {

public:

#ifdef HAVE_IFPACK_TEUCHOS
  //! Sets all parameters for the preconditioner.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;
#endif

  //! Sets integer parameters `Name' for the preconditioner.
  virtual int SetParameter(const string Name, const int Value) = 0;

  //! Sets double parameters `Name' for the preconditioner.
  virtual int SetParameter(const string Name, const double Value) = 0;

  //! Computes all it is necessary to apply the preconditioner.
  virtual int Compute() = 0;

  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  virtual bool IsComputed() const = 0;

  //! Returns the condition number estimate, computes it if necessary.
  virtual double Condest() const = 0;

  //! Applies the preconditioner to vector X, returns the result in Y.
  virtual int ApplyInverse(const Epetra_MultiVector& X,
			   Epetra_MultiVector& Y) const = 0;

  //! Returns a pointer to the matrix to be preconditioned.
  virtual const Epetra_RowMatrix& Matrix() const = 0;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const = 0;

};

inline ostream& operator<<(ostream& os, const Ifpack_Preconditioner& obj)
{
  return(obj.Print(os));
}

#endif // IFPACK_PRECONDITIONER_H
