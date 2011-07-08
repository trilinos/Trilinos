// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_HANDLEABLE_HPP
#define TEUCHOS_HANDLEABLE_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"


namespace Teuchos 
{
  /** \brief Class ConstHandleable provides an abstract interface for
   * polymorphic conversion from raw pointers to const smart pointers.
   *
   * Recall from the Teuchos RCP documentation that one should never
   * create directly a smart pointer from a raw pointer; rather, smart
   * pointers should be created through a call to rcp(). The type of the
   * argument to rcp() must be known at compile time. This makes the syntax
   * \code ConstHandle h = new Derived(); \endcode impossible with the
   * straightforward implementation in which ConstHandle takes a raw pointer
   * to a Base. In order to preserve this clean syntax, we require any handles
   * supporting this syntax to take a raw pointer to a ConstHandleable<Base>,
   * where ConstHandleable<Base> provides a getConstRcp() method which returns
   * the result of a call to rcp() on this.
   */
  template <typename Base> 
  class ConstHandleable 
  {
  public:
    /** \brief . */
    virtual ~ConstHandleable(){}

    /** \brief Virtual dtorReturn a safely-created RCP to the base
     * type */
    virtual RCP<const Base> getConstRcp() const = 0 ;
  };

  /** \brief Class Handleable provides an abstract interface for polymorphic
   * conversion from raw pointers to smart pointers.
   *
   * Recall from the Teuchos RCP documentation that one should never
   * create directly a smart pointer from a raw pointer; rather, smart
   * pointers should be created through a call to rcp(). The type of the
   * argument to rcp() must be known at compile time. This makes the syntax
   * \code Handle h = new Derived(); \endcode impossible with the
   * straightforward implementation in which Handle takes a raw pointer to a
   * Base. In order to preserve this clean syntax, we require any handles
   * supporting this syntax to take a raw pointer to a Handleable<Base>, where
   * Handleable<Base> provides a getRcp() method which returns the result of a
   * call to rcp() on this.
   */
  template <typename Base> 
  class Handleable : public virtual ConstHandleable<Base>
  {
  public:

    /** \brief . */
    virtual ~Handleable(){;}

    /** \brief Return a safely-created RCP to the base type */
    virtual RCP<Base> getRcp() = 0 ;

  };
}


/** \brief Use this macro as an easy way to implement the Handleable interface
 * in a derived class.
 *
 * For example,
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
  virtual Teuchos::RCP<const Base > getConstRcp() const {return rcp(this);} \
  virtual Teuchos::RCP<Base > getRcp() {return rcp(this);} 

/** \brief Use this macro as an easy way to implement the ConstHandleable
 * interface in a derived class. For example,
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
virtual Teuchos::RCP<const Base > getConstRcp() const {return rcp(this);}




#endif // TEUCHOS_HANDLEABLE_HPP
