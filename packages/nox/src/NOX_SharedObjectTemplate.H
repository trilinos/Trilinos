// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SHARED_OBJECT_TEMPLATE_H
#define NOX_SHARED_OBJECT_TEMPLATE_H

#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"

namespace NOX {

  /*! \brief Holder for objects that are shared between NOX::Abstract::Groups.
   *
   * Due to the size of certain operators (i.e. the Jacobian and
   * Preconditioning Matrices), we cannot afford to have
   * multiple copies. Instead we implement a shared object class that
   * that all groups use.
   */
template <class Object, class Owner>
class SharedObject {

public:

  //! Constructor
  SharedObject(const Teuchos::RCP<Object>& newObject) :
    object(newObject), owner(NULL) {};

  //! Destructor
  ~SharedObject(){};

  //! Get a non-const reference to the underlying object
  Teuchos::RCP<Object>
  getObject(const Owner* newOwner)
  {
    owner = newOwner;
    return object;
  };

  //! Return a const reference to the underlying object
  Teuchos::RCP<const Object> getObject() const
  { return object;};

  //! Return true if testOwner is the owner of the Jacobian
  bool isOwner(const Owner* checkOwner) const
  {return (owner == checkOwner);};

protected:

  //! Pointer to the shared object.
  Teuchos::RCP<Object> object;

  //! Pointer to const object that owns the shared object.
  const Owner* owner;
};

} // namespace NOX

#endif
