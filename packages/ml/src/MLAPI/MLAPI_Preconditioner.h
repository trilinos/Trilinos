#ifndef MLAPI_PRECONDITIONER_H
#define MLAPI_PRECONDITIONER_H
#include "MLAPI_BaseObject.h"

namespace MLAPI {

class DoubleVector;
class Space;

/*!
\class Preconditioner

\brief Pure virtual class for MLAPI preconditioners.

This class contains the basic methods that must be defined by all
preconditioners in the MLAPI workspace. By deriving from this class, any
preconditioner can be used for AztecOO's Krylov solvers through the
MLAPI::EpetraPreconditioner class.

\author Marzio Sala, SNL 9214.

\date Last modified on 07-Jan-05.

*/

class Preconditioner : public BaseObject{

public:

  //! Virtual destructor.
  virtual ~Preconditioner() {}

  //! Applies the preconditioner to LHS, using RHS as starting solution.
  virtual int Solve(const DoubleVector& LHS, DoubleVector& RHS) const = 0;

  //! Returns a copy of the domain space of \c this object.
  virtual const Space DomainSpace() const = 0;

  //! Returns a copy of the range space of \c this object.
  virtual const Space RangeSpace() const = 0;

};
} // namespace MLAPI

#endif
