// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PETRA_SHARED_JACOBIAN_H
#define _NLS_PETRA_SHARED_JACOBIAN_H

#include "Epetra_CrsMatrix.h"

namespace NOX {
namespace Epetra {

class Group;			// forward declaration

  /*! \brief Shared Jacobian for NOX Epetra Group objects.
   * 
   * Due to the size of the Jacobian, we cannot afford to have
   * multiple copies. Instead we implement a shared Jacobian. */

class SharedJacobian {

public:

  //! Constructor
  SharedJacobian(Epetra_RowMatrix& j);

  //! Destructor
  ~SharedJacobian();

  //! Get a reference to the Jacobian AND take ownership
  Epetra_RowMatrix& getJacobian(const Group* newowner);

  //! Return true if grp is the owner of the Jacobian
  bool isOwner(const Group* grp) const;

  //! Return a const pointer to the Jacobian
  const Epetra_RowMatrix& getJacobian() const;

protected:

  //! Pointer to Epetra Jacobian matrix 
  Epetra_RowMatrix* jacobian;

  //! Pointer to const group that owns Jacobian
  const Group* owner;
};
}
}

#endif
