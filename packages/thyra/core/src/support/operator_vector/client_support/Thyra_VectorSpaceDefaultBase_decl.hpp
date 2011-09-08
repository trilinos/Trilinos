// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
