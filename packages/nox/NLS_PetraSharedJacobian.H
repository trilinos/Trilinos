// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PETRA_SHARED_JACOBIAN_H
#define _NLS_PETRA_SHARED_JACOBIAN_H

//#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

class NLS_PetraGroup;		// forward declaration

//! shared Jacobian for NLS_PetraGroup objects

/*! Due to the size of the Jacobian, we cannot afford to have multiple
  copies. Instead we implement a shared Jacobian. */

class NLS_PetraSharedJacobian {

public:

  //! Constructor
  NLS_PetraSharedJacobian(Epetra_RowMatrix& j);

  //! Deconstructor
  ~NLS_PetraSharedJacobian();

  //! Get a reference to the Jacobian AND take ownership
  Epetra_RowMatrix& getJacobian(const NLS_PetraGroup* newowner);

  //! Give ownership to newowner
  void takeOwnership(const NLS_PetraGroup* newowner);

  //! Return true is grp is the owner of the Jacobian
  bool isOwner(const NLS_PetraGroup* grp) const;

  //! Return a const pointer to the Jacobian
  const Epetra_RowMatrix& getJacobian() const;

protected:

  //! Pointer to Epetra Jacobian Matrix 
  Epetra_RowMatrix* jacobian;

  //! Pointer to const NLS_PetraGroup that owns Jacobian
  const NLS_PetraGroup* owner;
};


#endif
