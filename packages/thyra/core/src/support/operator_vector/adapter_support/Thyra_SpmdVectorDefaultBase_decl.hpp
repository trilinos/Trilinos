// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SPMD_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_SPMD_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_VectorDefaultBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase_decl.hpp"

//#define THYRA_SPMD_VECTOR_BASE_DUMP

namespace Thyra {

/** \brief Base class for SPMD vectors that can provide views of contiguous
 * elements in a process.
 *
 * By inheriting from this base class, vector implementations allow their
 * vector objects to be seamlessly combined with other SPMD vector objects (of
 * potentially different concrete types) in <tt>applyOp()</tt>.  A big part of
 * this protocol is that every vector object can expose an
 * <tt>SpmdVectorSpaceBase</tt> object through the virtual function
 * <tt>spmdSpace()</tt>.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must override only two functions:
 * <tt>spmdSpaceImpl()</tt> and <tt>getLocalData(Scalar**,Ordinal*)</tt>.  The
 * default implementation of <tt>getLocalData(cons Scalar**,Ordinal*)</tt>
 * should rarely need to be overridden as it just calls the pure-virtual
 * non-<tt>const</tt> version.  Note that providing an implementation for
 * <tt>spmdSpaceImpl()</tt> of course means having to implement or use a
 * pre-implemented <tt>SpmdVectorSpaceBase</tt> subclass.
 *
 * Vector subclasses must also call the protected function
 * <tt>updateSpmdState()</tt> whenever the state of
 * <tt>*this->spmdSpaceImpl()</tt> vector space changes.  This function
 * gathers some cached data that makes the rest of the class more efficient.
 * This function must be called in a constructor or any other function that
 * changes the state of the vector space.
 *
 * Note that vector subclasses derived from this node interface class must
 * only be directly used in SPMD mode to work properly.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_support_grp
 */
template <class Scalar>
class SpmdVectorDefaultBase
  : virtual public SpmdVectorBase<Scalar>,
    virtual public VectorDefaultBase<Scalar> {
 public:
  /** @name Public interface functions */
  //@{

  /** \brief . */
  SpmdVectorDefaultBase();

  //@}

  /** \brief Implementation of applyOpImpl(...) that uses an input Comm.
   *
   * \param comm [in] The Spmd communicator to use in the global reduction.
   * If <tt>comm==NULL</tt>, then the local communicator will be used instead.
   */
  virtual void applyOpImplWithComm(
      const Ptr<const Teuchos::Comm<Ordinal> > &comm,
      const RTOpPack::RTOpT<Scalar> &op,
      const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
      const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
      const Ptr<RTOpPack::ReductTarget> &reduct_obj,
      const Ordinal global_offset) const;

  //@}

  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Overridden public functions from VectorBase */
  //@{

  /** \brief Returns <tt>this->spmdSpace()</tt>. */
  Teuchos::RCP<const VectorSpaceBase<Scalar> > space() const;

  //@}

 protected:
  /** \name Overridden protected functions from VectorBase */
  //@{

  /** \brief Calls applyOpImplWithComm(null,op,...). */
  void applyOpImpl(
      const RTOpPack::RTOpT<Scalar> &op,
      const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
      const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
      const Ptr<RTOpPack::ReductTarget> &reduct_obj,
      const Ordinal global_offset) const;
  /** \brief Implemented through <tt>this->getLocalData()</tt> */
  void acquireDetachedVectorViewImpl(
      const Range1D &rng, RTOpPack::ConstSubVectorView<Scalar> *sub_vec) const;
  /** \brief Implemented through <tt>this->freeLocalData()</tt> */
  void releaseDetachedVectorViewImpl(
      RTOpPack::ConstSubVectorView<Scalar> *sub_vec) const;
  /** \brief Implemented through <tt>this->getLocalData()</tt> */
  void acquireNonconstDetachedVectorViewImpl(
      const Range1D &rng, RTOpPack::SubVectorView<Scalar> *sub_vec);
  /** \brief Implemented through <tt>this->commitLocalData()</tt> */
  void commitNonconstDetachedVectorViewImpl(
      RTOpPack::SubVectorView<Scalar> *sub_vec);

  //@}

  /** \name Overridden Protected functions from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RTOpPack::SubMultiVectorView<Scalar> getNonconstLocalSubMultiVectorImpl();
  /** \brief . */
  RTOpPack::ConstSubMultiVectorView<Scalar> getLocalSubMultiVectorImpl() const;
  /** \brief . */
  void getNonconstLocalMultiVectorDataImpl(
      const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim);
  /** \brief . */
  void getLocalMultiVectorDataImpl(
      const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim) const;
  //@}

  /** \name Overridden Protected functions from SpmdVectorBase */
  //@{
  /** \brief Virtual implementation for getNonconstLocalSubVector(). */
  RTOpPack::SubVectorView<Scalar> getNonconstLocalSubVectorImpl();
  /** \brief Virtual implementation for getLocalSubVector(). */
  RTOpPack::ConstSubVectorView<Scalar> getLocalSubVectorImpl() const;
  //@}

  /** \name Protected functions to be used by subclasses */
  //@{

  /** \brief Subclasses must call this function whenever the structure of the
   * VectorSpaceBase changes.
   */
  virtual void updateSpmdSpace();

  //@}

 private:
  // ///////////////////////////////////////
  // Private data members

  // Cached (only on vector space!)
  mutable Ordinal globalDim_;
  mutable Ordinal localOffset_;
  mutable Ordinal localSubDim_;

  // /////////////////////////////////////
  // Private member functions

  Range1D validateRange(const Range1D &rng_in) const;

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
 public:
  static bool show_dump;
#endif  // THYRA_SPMD_VECTOR_BASE_DUMP

};  // end class SpmdVectorDefaultBase

}  // end namespace Thyra

#endif  // THYRA_SPMD_VECTOR_DEFAULT_BASE_DECL_HPP
