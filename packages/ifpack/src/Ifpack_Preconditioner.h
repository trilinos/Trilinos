#ifndef IFPACK_PRECONDITIONER_H
#define IFPACK_PRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
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
  only two additional methods:
  - Compute(), that should compute all is required to apply the
    preconditioner;
  - IsComputed(), that should return true is the preconditioner
    has been successfully computed, false otherwise.

  \date Sep-04
*/

class Ifpack_Preconditioner : public Epetra_Operator {

public:

#ifdef HAVE_IFPACK_TEUCHOS
  //! Sets all parameters for the preconditioner.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;
#endif

  //! Computes all it is necessary to apply the preconditioner.
  virtual int Compute() = 0;

  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  virtual bool IsComputed() const = 0;

  //! Returns a pointer to the matrix to be preconditioned.
  virtual const Epetra_RowMatrix* Matrix() const = 0;

  // FIXME: copy constructors...
  // FIXME cout << 
};

#endif // IFPACK_PRECONDITIONER_H
