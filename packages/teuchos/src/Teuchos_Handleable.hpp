// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_HANDLEABLE_HPP
#define TEUCHOS_HANDLEABLE_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Teuchos 
{
  /** 
   * Class ConstHandleable provides an abstract interface for polymorphic
   * conversion from raw pointers to const smart pointers. Recall from the
   * Teuchos RefCountPtr documentation that one should never create
   * directly a smart pointer from a raw pointer; rather, smart pointers
   * should be created through a call to rcp(). The type of the argument
   * to rcp() must be known at compile time. This makes the syntax
   * \code
   * ConstHandle h = new Derived();
   * \endcode
   * impossible with the straightforward implementation in which 
   * ConstHandle takes
   * a raw pointer to a Base. In order to preserve this clean syntax, we
   * require any handles supporting this syntax to take a raw
   * pointer to a ConstHandleable<Base>, where ConstHandleable<Base> 
   * provides a getConstRcp() method which returns the result of 
   a call to rcp() on this. 
   */
  template <typename Base> 
  class ConstHandleable 
  {
  public:
    /** Virtual dtor */
    virtual ~ConstHandleable(){;}

    /** Return a safely-created RefCountPtr to the base type */
    virtual RefCountPtr<const Base> getConstRcp() const = 0 ;
  };

  /** 
   * Class Handleable provides an abstract interface for polymorphic
   * conversion from raw pointers to smart pointers. Recall from the
   * Teuchos RefCountPtr documentation that one should never create
   * directly a smart pointer from a raw pointer; rather, smart pointers
   * should be created through a call to rcp(). The type of the argument
   * to rcp() must be known at compile time. This makes the syntax
   * \code
   * Handle h = new Derived();
   * \endcode
   * impossible with the straightforward implementation in which Handle takes
   * a raw pointer to a Base. In order to preserve this clean syntax, we
   * require any handles supporting this syntax to take a raw
   * pointer to a Handleable<Base>, where Handleable<Base> provides a 
   * getRcp() method which returns the result of a call to rcp() on this.
   */
  template <typename Base> 
  class Handleable : public virtual ConstHandleable<Base>
  {
  public:
    /** Virtual dtor */
    virtual ~Handleable(){;}

    /** Return a safely-created RefCountPtr to the base type */
    virtual RefCountPtr<Base> getRcp() = 0 ;
  };
}


/** \def
 * Use this macro as an easy way to implement the Handleable interface
 * in a derived class. For example,
 *
 * \code
 * class Derived : public Handleable<Base>
 * {
 * public:
 * TEUCHOS_GET_RCP(Base);
 * };
 * \endcode
 */
#define TEUCHOS_GET_RCP(Base)                                           \
  /** Teuchos::ConstHandleable<##Base> interface */                     \
  virtual Teuchos::RefCountPtr<const Base> getConstRcp() const {return rcp(this);} \
  /** Teuchos::Handleable<##Base> interface */                          \
  virtual Teuchos::RefCountPtr<Base > getRcp() {return rcp(this);} 

/** \def 
 * Use this macro as an easy way to implement the ConstHandleable interface
 * in a derived class. For example,
 *
 * \code
 * class Derived : public ConstHandleable<Base>
 * {
 * public:
 * TEUCHOS_GET_CONST_RCP(Base);
 * };
 * \endcode
 */
#define TEUCHOS_GET_CONST_RCP(Base) \
/** Teuchos::ConstHandleable<##Base> interface */                            \
virtual Teuchos::RefCountPtr<const Base> getConstRcp() const {return rcp(this);}




#endif // TEUCHOS_HANDLEABLE_HPP
