// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
