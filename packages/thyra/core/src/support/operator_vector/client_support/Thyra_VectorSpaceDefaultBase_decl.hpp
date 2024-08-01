// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP
#define THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP

#include "Thyra_VectorSpaceBase_decl.hpp"


namespace Thyra {


/** \brief Node <tt>VectorSpaceBase</tt> subclass that provides default
 * implementations for many functions using a default multi-vectors
 * implementation.
 *
 * <b>Notes to Subclass Developers</b>
 *
 * Because of the default multi-vector implementation used in this node
 * subclasses, a concrete subclass is only required to override four
 * functions: <tt>dim()</tt>, <tt>isCompatible()</tt>, <tt>createMember()</tt>
 * and <tt>smallVecSpcFcty()</tt>.  Note that implementing the
 * <tt>createMember()</tt> method also entails defining a concrete
 * <tt>VectorBase</tt> subclass and defining <tt>smallVecSpcFcty()</tt>
 * entails defining a concrete <tt>VectorSpaceFactoryBase</tt> subclass.
 *
 * If a subclass can support specialized multi-vectors, then the
 * <tt>createMembers()</tt> function should be overridden as well.  Note that
 * implementing <tt>createMembers()</tt> also entails defining a concrete
 * <tt>MultiVectorBase</tt> subclass.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class VectorSpaceDefaultBase : virtual public VectorSpaceBase<Scalar> {
protected:

  /** @name Protected functions overridden from VectorSpaceBase */
  //@{

  /** \brief .
   *
   * The default implementation returns <tt>dynamic_cast<
   * DefaultColumnwiseMultiVector<Scalar> >(return.get())!=NULL</tt>.
   */
  RCP< MultiVectorBase<Scalar> > createMembers(int numMembers) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements and then only when the
   * vector is destroyed is the data copied out of the vector and back into
   * the elements pointed to by <tt>raw_v.values()</tt>.
   */
  RCP<VectorBase<Scalar> >
  createMemberView( const RTOpPack::SubVectorView<Scalar> &raw_v ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements from
   * <tt>raw_v.values()</tt>.
   */
  RCP<const VectorBase<Scalar> >
  createMemberView( const RTOpPack::ConstSubVectorView<Scalar> &raw_v ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMembers(raw_mv.numSubCols())</tt> to create a
   * multi-vector then uses the explicit element access functions to set the
   * elements and then only when the multi-vector is destroyed is the data
   * copied out of the multi-vector and back into the elements pointed to by
   * <tt>raw_mv.values()</tt>.
   */
  RCP<MultiVectorBase<Scalar> > 
  createMembersView( const RTOpPack::SubMultiVectorView<Scalar> &raw_mv ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMembers()</tt> to create a multi-vector then uses the
   * explicit element access functions to set the elements.
   */
  RCP<const MultiVectorBase<Scalar> >
  createMembersView( const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv ) const;

  //@}

};


} // end namespace Thyra


#endif  // THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP
