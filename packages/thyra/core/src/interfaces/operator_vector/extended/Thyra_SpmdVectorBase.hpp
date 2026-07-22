// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SPMD_VECTOR_BASE_DECL_HPP
#define THYRA_SPMD_VECTOR_BASE_DECL_HPP


#include "Thyra_VectorBase.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"


namespace Thyra {


/** \brief Base class for SPMD vectors that can provide views of contiguous
 * elements in a process.
 *
 * By inheriting from this base class, vector implementations allow their
 * vector objects to be seamlessly combined with other SPMD vector objects (of
 * potentially different concrete types) in <tt>applyOp()</tt>.  A big part of
 * this protocol is that every vector object can expose an
 * <tt>SpmdVectorSpaceBase</tt> object through the function
 * <tt>spmdSpace()</tt>.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class SpmdVectorBase :
    virtual public VectorBase<Scalar>,
    virtual public SpmdMultiVectorBase<Scalar>
{
public:

  /** \brief . */
  using SpmdMultiVectorBase<Scalar>::getNonconstLocalData;
  /** \brief . */
  using SpmdMultiVectorBase<Scalar>::getLocalData;

  /** @name Public non-virtual interface functions */
  //@{

  /** \brief Get a non-const generalized view of local vector data.
   */
  RTOpPack::SubVectorView<Scalar> getNonconstLocalSubVector()
    { return getNonconstLocalSubVectorImpl(); }

  /** \brief Get a const generalized view of local vector data.
   */
  RTOpPack::ConstSubVectorView<Scalar> getLocalSubVector() const
    { return getLocalSubVectorImpl(); }

  /** \brief Returns a non-<tt>const</tt> pointer to the beginning of the
   * local vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to an
   * array of the local values.
   *
   * Preconditions:<ul>
   * <li> <tt>nonnull(localValues)==true</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>nonnull(*localValues)==true</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be freed by removing
   * all of the <tt>ArrayRCP</tt> objects (or setting them to null).
   */
  void getNonconstLocalData(const Ptr<ArrayRCP<Scalar> > &localValues)
    { this->getNonconstLocalVectorDataImpl(localValues); }

  /** \brief Returns a <tt>const</tt> pointer to the beginning of the local
   * vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to an
   * array of the local values.
   *
   * Preconditions:<ul>
   * <li> <tt>nonnull(localValues)==true</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>nonnull(*localValues)==true</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be freed by removing
   * all of the <tt>ArrayRCP</tt> objects (or setting them to null).
   */
  void getLocalData(const Ptr<ArrayRCP<const Scalar> > &localValues) const
    { this->getLocalVectorDataImpl(localValues); }

  //@}

protected:

  /** \name Protected pure virtual functions to be overridden */
  //@{

  /** \brief Virtual implementation for getNonconstLocalSubVector(). */
  virtual RTOpPack::SubVectorView<Scalar>
  getNonconstLocalSubVectorImpl() = 0;

  /** \brief Virtual implementation for getLocalSubVector(). */
  virtual RTOpPack::ConstSubVectorView<Scalar>
  getLocalSubVectorImpl() const = 0;

  /** \brief Implementation of getNonconstLocalData() */
  virtual void getNonconstLocalVectorDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues) = 0;

  /** \brief Implementation of getLocalData() */
  virtual void getLocalVectorDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues) const = 0;

  //@}


}; // end class SpmdVectorBase


} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_BASE_DECL_HPP
