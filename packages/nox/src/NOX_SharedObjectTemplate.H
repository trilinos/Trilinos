// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
