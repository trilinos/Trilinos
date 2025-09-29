// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SPMD_VECTOR_DECL_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_DECL_HPP

#include "Thyra_SpmdVectorDefaultBase_decl.hpp"

namespace Thyra {

/** \brief Efficient concrete implementation subclass for SPMD vectors.
 *
 * This subclass provides a very efficient and very general concrete
 * implementation of a <tt>Thyra::VectorBase</tt> object for any SPMD
 * platform.
 *
 * Objects of this type generally should not be constructed directly by a
 * client but instead by using the concrete vector space subclass
 * <tt>Thyra::DefaultSpmdVectorSpace</tt> and using the function
 * <tt>Thyra::createMember()</tt>.
 *
 * The storage type can be anything since an <tt>ArrayRCP</tt> is used to pass
 * in the local values pointer into the constructor and <tt>initialize()</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template <class Scalar>
class DefaultSpmdVector : virtual public SpmdVectorDefaultBase<Scalar> {
 public:
  /** @name Constructors/initializers */
  //@{

  /** \brief Construct to uninitialized. */
  DefaultSpmdVector();

  /** \brief Calls <tt>initialize()</tt>. */
  DefaultSpmdVector(
      const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdSpace, const ArrayRCP<Scalar> &localValues, const Ordinal stride);

  /** \brief Initialize.
   *
   * \param  spmdSpace
   *           [in] Smart pointer to <tt>SpmdVectorSpaceBase</tt> object
   *           that defines the data distribution for <tt>spmdSpace()</tt> and <tt>space()</tt>.
   * \param  localValues
   *           [in] Smart pointer to beginning of local strided vector data.
   *           This array must be at least of dimension <tt>mpiRangeSpace->localDim()*stride</tt>
   *           and <tt>(&*localValues)[ i*stride ]</tt> gives the local value
   *           of the zero-based entry <tt>(i)</tt> where <tt>i=0...spmdSpace()->localSubDim()-1</tt>.
   * \param  stride
   *           [in] Stride between local vector elements.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>spmdSpace.get()!=NULL</tt>
   * <li><tt>localValues.get()!=NULL</tt>
   * <li><tt>stride != 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>this->getRCptr().get() == localValues.get()</tt>
   * <li> <tt>this->getPtr() == &*localValues</tt>
   * <li> <tt>this->getStride() == stride</tt>
   * </ul>
   */
  void initialize(
      const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdSpace, const ArrayRCP<Scalar> &localValues, const Ordinal stride);

  /** \brief Set to an uninitialized state.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->spmdSpace().get() == NULL</tt>.
   */
  void uninitialize(
      RCP<const SpmdVectorSpaceBase<Scalar> > *spmdSpace = NULL, ArrayRCP<Scalar> *localValues = NULL, Ordinal *stride = NULL);

  //@}

  /** @name Accessors (inlined for minimal overhead) */
  //@{

  /** \brief . */
  ArrayRCP<Scalar> getRCPtr();
  /** \brief . */
  ArrayRCP<const Scalar> getRCPtr() const;
  /** \brief . */
  Scalar *getPtr();
  /** \brief . */
  const Scalar *getPtr() const;
  /** \brief . */
  Ordinal getStride() const;

  //@}

  /** @name Overridden from SpmdMultiVectorBase */
  //@{

  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpaceImpl() const;

  //@}

  /** @name Overridden from SpmdVectorBase */
  //@{
  /** \brief . */
  void getNonconstLocalVectorDataImpl(const Ptr<ArrayRCP<Scalar> > &localValues);
  /** \brief . */
  void getLocalVectorDataImpl(const Ptr<ArrayRCP<const Scalar> > &localValues) const;

  //@}

 private:
  // ///////////////////////////////////////
  // Private data members

  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace_;
  ArrayRCP<Scalar> localValues_;
  Ordinal stride_;
};

// /////////////////////////////////////////////////////
// Inline members

template <class Scalar>
inline ArrayRCP<Scalar>
DefaultSpmdVector<Scalar>::getRCPtr() {
  return localValues_;
}

template <class Scalar>
inline ArrayRCP<const Scalar>
DefaultSpmdVector<Scalar>::getRCPtr() const {
  return localValues_;
}

template <class Scalar>
inline Scalar *DefaultSpmdVector<Scalar>::getPtr() {
  return localValues_.get();
}

template <class Scalar>
inline const Scalar *DefaultSpmdVector<Scalar>::getPtr() const {
  return localValues_.get();
}

template <class Scalar>
inline Ordinal DefaultSpmdVector<Scalar>::getStride() const {
  return stride_;
}

}  // end namespace Thyra

#endif  // THYRA_DEFAULT_SPMD_VECTOR_DECL_HPP
