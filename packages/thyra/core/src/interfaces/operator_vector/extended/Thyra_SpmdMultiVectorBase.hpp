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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SPMD_MULTI_VECTOR_BASE_DECL_HPP
#define THYRA_SPMD_MULTI_VECTOR_BASE_DECL_HPP


#include "Thyra_MultiVectorBase.hpp"


namespace Thyra {


template<class Scalar> class SpmdVectorSpaceBase;


/** \brief Base interface class for SPMD multi-vectors.
 *
 * By inheriting from this base class, multi-vector implementations allow
 * their multi-vector objects to be seamlessly combined with other SPMD
 * multi-vector objects (of different concrete types) in <tt>applyOp()</tt>
 * and <tt>apply()</tt>.  A big part of this protocol is that every
 * multi-vector object can expose an <tt>SpmdVectorSpaceBase</tt> object
 * through the function <tt>spmdSpace()</tt>.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class SpmdMultiVectorBase : virtual public MultiVectorBase<Scalar>
{
public:

  /** @name Public non-virtual interface functions */
  //@{

  /** \brief Returns the SPMD vector space object for the range of
   * <tt>*this</tt> multi-vector.
   */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const
    { return spmdSpaceImpl(); }

  /** \brief Get a non-const generalized view of local multi-vector data.
   */
  RTOpPack::SubMultiVectorView<Scalar> getNonconstLocalSubMultiVector()
    { return getNonconstLocalSubMultiVectorImpl(); }

  /** \brief Get a const generalized view of local multi-vector data.
   */
  RTOpPack::ConstSubMultiVectorView<Scalar> getLocalSubMultiVector() const
    {  return getLocalSubMultiVectorImpl(); }

  /** \brief Returns a non-<tt>const</tt> pointer to a Fortran-style view of
   * the local multi-vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to
   * the first element in the first column of the local multi-vector stored as
   * a column-major dense Fortran-style matrix.
   *
   * \param leadingDim [out] On output <tt>*leadingDim</tt> gives the leading
   * dimension of the Fortran-style local multi-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>leadingDim!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*leadingDim!=0</tt>
   * </ul>
   */
  void getNonconstLocalData(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    )
    { getNonconstLocalMultiVectorDataImpl(localValues, leadingDim); }

  /** \brief Returns a <tt>const</tt> pointer to a Fortran-style view of the
   * local multi-vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to
   * the first element in the first column of the local multi-vector stored as
   * a column-major dense Fortran-style matrix.
   *
   * \param leadingDim [out] On output <tt>*leadingDim</tt> gives the leading
   * dimension of the Fortran-style local multi-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>leadingDim!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*leadingDim!=0</tt>
   * </ul>
   */
  void getLocalData(
    const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const
    { getLocalMultiVectorDataImpl(localValues, leadingDim); }

  //@}

protected:

  /** @name Virtual functions to be overridden by sublcasses. */
  //@{

  /** \brief Virtual implementation for spmdSpace(). */
  virtual RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpaceImpl() const = 0;

  /** \brief Virtual implementation for getNonconstLocalSubMultiVector(). */
  virtual RTOpPack::SubMultiVectorView<Scalar>
  getNonconstLocalSubMultiVectorImpl() = 0;

  /** \brief Virtual implementation for getLocalSubMultiVector(). */
  virtual RTOpPack::ConstSubMultiVectorView<Scalar>
  getLocalSubMultiVectorImpl() const = 0;

  /** \brief Virtual implementation for getNonconstLocalData(). */
  virtual void getNonconstLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) = 0;

  /** \brief Virtual implementation for getLocalData(). */
  virtual void getLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const = 0;

  //@}
  
}; // end class SpmdMultiVectorBase


} // end namespace Thyra


#endif // THYRA_SPMD_MULTI_VECTOR_BASE_DECL_HPP
