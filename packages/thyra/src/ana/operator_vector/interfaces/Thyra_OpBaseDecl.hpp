// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_OP_BASE_DECL_HPP
#define THYRA_OP_BASE_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Describable.hpp"

namespace Thyra {

/** \brief Base class for all operators.
 *
 * It is not expected that clients will manipulate objects directly
 * through this interface.  This interface is just meant to provide
 * common declarations the functions <tt>domain()</tt>,
 * <tt>range()</tt> and <tt>opSupported()</tt>.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class RangeScalar, class DomainScalar = RangeScalar>
class OpBase : virtual public Teuchos::Describable {
public:

  /** \brief . */
  using Teuchos::Describable::describe;

  /** \brief . */
  virtual ~OpBase() {}

  /** @name Pure virtual methods (must be overridden by subclass) */

  /** \brief Return a smart pointer for the range space for <tt>this</tt> operator.
   *
   * Note that a return value of <tt>return.get()==NULL</tt> is a flag
   * that <tt>*this</tt> is not fully initialized.
   *
   * If <tt>return.get()!=NULL</tt>, it is required that the object
   * referenced by <tt>*return.get()</tt> must have lifetime that
   * extends past the lifetime of the returned smart pointer object.
   * However, the object referenced by <tt>*return.get()</tt> may
   * change if <tt>*this</tt> modified so this reference should not
   * be maintained for too long.
   *
   * Once more, the client should not expect the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> to be valid past the lifetime
   * of <tt>*this</tt>.
   */
  virtual Teuchos::RefCountPtr< const VectorSpaceBase<RangeScalar> > range() const = 0;

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   *
   * Note that a return value of <tt>return.get()==NULL</tt> is a flag
   * that <tt>*this</tt> is not fully initialized.
   *
   * If <tt>return.get()!=NULL</tt>, it is required that the object
   * referenced by <tt>*return.get()</tt> must have lifetime that
   * extends past the lifetime of the returned smart pointer object.
   * However, the object referenced by <tt>*return.get()</tt> may
   * change if <tt>*this</tt> modified so this reference should not
   * be maintained for too long.
   *
   * Once more, the client should not expect the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> to be valid past the lifetime
   * of <tt>*this</tt>.
   */
  virtual Teuchos::RefCountPtr< const VectorSpaceBase<DomainScalar> > domain() const = 0;

  //@}

};	// end class OpBase

}	// end namespace Thyra

#endif	// THYRA_OP_BASE_DECL_HPP
