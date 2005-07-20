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

#ifndef THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP
#define THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP

#include "Thyra_VectorSpaceBaseDecl.hpp"

namespace Thyra {

/** \brief Node subclass that provides a default implementation for
 * multi-vectors.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
class VectorSpaceDefaultBase : virtual public VectorSpaceBase<Scalar> {
public:

  /** @name Public functions overridden from VectorSpaceBase */ 
  //@{

  /** \brief .
   *
   * The default implementation returns
   * <tt>dynamic_cast<SerialVectorSpaceFactory>(return.get())!=NULL</tt>.
   * Note that if a subclass overrides <tt>createMembers()</tt> then
   * it may also need to override this method as well.
   */
  Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const;

  //@}

protected:

  /** @name Protected functions overridden from VectorSpaceBase */
  //@{

  /** \brief .
   *
   * The default implementation returns
   * <tt>dynamic_cast<MultiVectorCols>(return.get())!=NULL</tt>.
   */
  Teuchos::RefCountPtr< MultiVectorBase<Scalar> > createMembers(int numMembers) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements and then only when
   * the vector is destroyed is the data copied out of the vector and
   * back into the elements pointed to by
   * <tt>raw_v.values()</tt>.
   */
  Teuchos::RefCountPtr<VectorBase<Scalar> > createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements.
   */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>createMembers(raw_mv.numSubCols())</tt> to create a
   * multi-vector then uses the explicit element access functions to
   * set the elements and then only when the multi-vector is destroyed
   * is the data copied out of the multi-vector and back into the
   * elements pointed to by <tt>raw_mv.values()</tt>.
   */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>createMembers()</tt> to create a multi-vector then uses the explicit
   * element access functions to set the elements.
   */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const;

  //@}

}; // end class VectorSpaceBase

} // end namespace Thyra

#endif  // THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP
