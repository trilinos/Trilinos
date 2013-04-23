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

#ifndef THYRA_SPMD_VECTOR_BASE_DECL_HPP
#define THYRA_SPMD_VECTOR_BASE_DECL_HPP


#include "Thyra_VectorDefaultBase.hpp"
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
 * This base class contains an implementation of <tt>applyOp()</tt> that
 * relies on implementations of the <tt>const</tt> functions
 * <tt>acquireDetachedView()</tt> and <tt>releaseDetachedView()</tt>, and the
 * non-<tt>const</tt> functions <tt>acquireDetachedView()</tt> and
 * <tt>commitDetachedView()</tt> (which all have default implementations in
 * this subclass).  In essence, this implementation will only call the
 * <tt>acquireDetachedView()</tt> functions using a range of (global) indexes
 * for elements that exist in the local process.  As long as the number of
 * local elements in each process is fairly large, the virtual function call
 * overhead will be minimal and this will result in a near optimal
 * implementation.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must override only two functions: <tt>spmdSpace()</tt>
 * and <tt>getLocalData(Scalar**,Ordinal*)</tt>.  The default implementation of
 * <tt>getLocalData(cons Scalar**,Ordinal*)</tt> should rarely need to be
 * overridden as it just calls the pure-virtual non-<tt>const</tt> version.
 * Note that providing an implementation for <tt>spmdSpace()</tt> of course
 * means having to implement or use a pre-implemented
 * <tt>SpmdVectorSpaceBase</tt> subclass.
 *
 * Vector subclasses must also call the protected function
 * <tt>updateSpmdState()</tt> whenever the state of <tt>*this->spmdSpace()</tt>
 * vector space changes.  This function gathers some cached data that makes
 * the rest of the class more efficient.  This function must be called in a
 * constructor or any other function that changes the state of the vector
 * space.
 *
 * If the <tt>acquireDetachedView()</tt> functions are ever called with index
 * ranges outside of those of the local process, then the default
 * implementations in <tt>VectorBase</tt> of all of the functions
 * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt>, and
 * <tt>commitDetachedView()</tt> are called instead.  Alternatively, a
 * subclass could provide more specialized implementations of these functions
 * (for more efficient gather/scatter operations) if desired but this should
 * not be needed for most use cases.
 *
 * It is interesting to note that in the above use case when the explicit
 * subvector access functions call on its default implementation defined in
 * <tt>VectorBase</tt> (which calls on <tt>applyOpImpl()</tt>), the operator will
 * be properly applied since the version of <tt>applyOpImpl()</tt> implemented in
 * this class will only request local vector data and hence there will only be
 * two levels of recursion for any call to an explicit subvector access
 * function.  This is a truly elegant result.
 *
 * Note that vector subclasses derived from this node interface class must
 * only be directly used in SPMD mode to work properly.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_support_grp
 */
template<class Scalar>
class SpmdVectorBase : virtual public VectorDefaultBase<Scalar> {
public:

  /** @name Public interface functions */
  //@{

  /** \brief . */
  SpmdVectorBase();

  /** \brief Returns the Spmd-based vector space object for <tt>*this</tt> vector.
   */
  virtual Teuchos::RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const = 0;

  /** \brief Get a non-const view of the local data.
   *
   * ToDo: Refactor this interface to use iterator access.
   */
  RTOpPack::SubVectorView<Scalar> getNonconstLocalSubVector();

  /** \brief Get a const view of the local data.
   *
   * ToDo: Refactor this interface to use iterator access.
   */
  RTOpPack::ConstSubVectorView<Scalar> getLocalSubVector() const;

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
    { this->getNonconstLocalDataImpl(localValues); }

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
    { this->getLocalDataImpl(localValues); }

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
    const Ordinal global_offset
    ) const;

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

  /** \name Protected pure virtual functions to be overridden */
  //@{

  /** \brief Implementation of getNonconstLocalData() */
  virtual void getNonconstLocalDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues) = 0;

  /** \brief Implementation of getLocalData() */
  virtual void getLocalDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues) const = 0;

  //@}

  /** \name Overridden protected functions from VectorBase */
  //@{

  /** \brief Calls applyOpImplWithComm(null,op,...). */
  void applyOpImpl(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const;
  /** \brief Implemented through <tt>this->getLocalData()</tt> */
  void acquireDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief Implemented through <tt>this->freeLocalData()</tt> */
  void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief Implemented through <tt>this->getLocalData()</tt> */
  void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief Implemented through <tt>this->commitLocalData()</tt> */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    );

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
  
  mutable bool in_applyOpImpl_;

  // Cached (only on vector space!)
  mutable Ordinal  globalDim_;
  mutable Ordinal  localOffset_;
  mutable Ordinal  localSubDim_;

  // /////////////////////////////////////
  // Private member functions

  Range1D validateRange( const Range1D& rng_in ) const;

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
public:
  static bool show_dump;
#endif // THYRA_SPMD_VECTOR_BASE_DUMP

}; // end class SpmdVectorBase


} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_BASE_DECL_HPP
