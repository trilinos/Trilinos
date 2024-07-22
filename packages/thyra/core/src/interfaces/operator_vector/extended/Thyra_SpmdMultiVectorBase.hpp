// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
