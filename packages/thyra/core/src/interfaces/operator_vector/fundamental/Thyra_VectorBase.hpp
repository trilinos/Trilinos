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

#ifndef THYRA_VECTOR_BASE_DECL_HPP
#define THYRA_VECTOR_BASE_DECL_HPP


#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_MultiVectorBase_decl.hpp"
#include "RTOpPack_RTOpT.hpp"
#include "RTOpPack_SparseSubVectorT.hpp"


namespace Thyra {


/** \brief Abstract interface for finite-dimensional dense vectors.
 *
 * This interface contains the minimal set of operations needed to define an
 * abstract vector.
 *
 * \section Thyra_VB_outline_sec Outline
 *
 * <ul>
 * <li>\ref Thyra_VB_rtop_sec
 * <li>\ref Thyra_VB_rtop_collection_sec
 * <li>\ref Thyra_VB_expl_access_sec
 * <li>\ref Thyra_VB_expl_access_utils_sec
 * <li>\ref Thyra_VB_expl_access_assign_sec
 * <li>\ref Thyra_VB_is_lin_op_sec
 * <li>\ref Thyra_VB_dev_notes_sec
 * </ul>
 *
 * \section Thyra_VB_rtop_sec Reduction/transformation operator (RTOp) support
 *
 * The main feature of this interface is the function <tt>applyOp()</tt>
 * which is used to implement all types of vector reduction and transformation
 * operations (RTOp) through RTOp operators .  Every standard (i.e. BLAS) and
 * nearly every non-standard element-wise operation that can be performed on a
 * set of vectors can be performed efficiently through
 * reduction/transformation operators.  More standard vector operations could
 * be included in this interface and allow for specialized implementations
 * but, in general, assuming the sub-vectors are large enough, such
 * implementations would not be significantly faster than those implemented
 * through reduction/transformation operators.  There are some operations
 * however that can not always be efficiently implemented with
 * reduction/transformation operators and a few of these important operations
 * are included in this interface.  The <tt>applyOp()</tt> function allows to
 * client to specify a sub-set of the vector elements to include in
 * reduction/transformation operation.  This greatly increases the generality
 * of this vector interface as vector objects can be used as sub objects in
 * larger composite vectors and sub-views of a vector can be created.
 *
 * \section Thyra_VB_rtop_collection_sec Collection of pre-written RTOps and wrapper functions
 *
 * There already exists RTOp-based implementations of several standard vector
 * operations and some convenience functions that wrap these operators and
 * call <tt>applyOp()</tt>.  These wrapper functions can be found
 * <a href="../../../../../../support/operator_vector/doc/html/group__Thyra__Op__Vec__VectorStdOps__grp.html">here</a>
 *
 * \section Thyra_VB_expl_access_sec Explicit vector coefficient access
 *
 * This interface also allows a client to extract a sub-set of vector
 * coefficients in an explicit form as non-mutable
 * <tt>RTOpPack::ConstSubVectorView</tt> or mutable
 * <tt>RTOpPack::SubVectorView</tt> objects using the
 * <tt>acquireDetachedView()</tt> functions.  In general, this is a very
 * inefficient thing to do and should be avoided.  However, there are some
 * situations where getting explicit access to the coefficients of a vector is
 * a very reasonable and efficient thing to do (i.e. for vectors in the domain
 * of a multi-vector for instance) and therefore this functionality is
 * supported.  These views and the parent vector follow the state behavior
 * outlined \ref Thyra_Op_Vec_Behavior_Of_Views_grp "here".
 *
 * \section Thyra_VB_expl_access_utils_sec Explicit vector coefficient access utilities
 *
 * Note that client code in general should not directly call the above
 * explicit sub-vector access functions but should use the utility classes
 * <tt>ConstDetachedVectorView</tt> and <tt>DetachedVectorView</tt> instead
 * since these are easier an safer in the event that an exception is thrown.
 *
 * \section Thyra_VB_expl_access_assign_sec Explicit vector coefficient assignment
 *
 * In addition to being able to extract an explicit non-mutable and
 * mutable views of some (small?) sub-set of elements, this interface
 * allows a client to set sub-vectors using <tt>setSubVector()</tt>.
 *
 * \section Thyra_VB_is_lin_op_sec Vector is a MultiVectorBase is a LinearOpBase
 *
 * It is also worth mentioning that that this <tt>%VectorBase</tt> interface
 * class also inherits from <tt>MultiVectorBase</tt> so that every
 * <tt>%VectorBase</tt> object is also a <tt>%MultiVectorBase</tt> object.
 * This allows any piece of code that accepts <tt>%MultiVectorBase</tt>
 * objects to automatically accept <tt>%VectorBase</tt> objects as well.  In
 * addition, since <tt>MultiVectorBase</tt> inherits from
 * <tt>LinearOpBase</tt>, then this means that every vector is also a linear
 * operator.
 *
 * \section Thyra_VB_dev_notes_sec Notes for subclass developers
 *
 * The support subclass <tt>VectorDefaultBase</tt> provides default
 * implementations for as many functions as possible and should be considered
 * a first choice for creating concrete subclasses.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class VectorBase : virtual public MultiVectorBase<Scalar>
{
public:

#ifdef THYRA_INJECT_USING_DECLARATIONS
  using MultiVectorBase<Scalar>::apply;
#endif

  /** @name Space membership */
  //@{

  /** \brief Return a smart pointer to the vector space that this vector
   * belongs to.
   *
   * A return value of <tt>space().get()==NULL</tt> is a flag that
   * <tt>*this</tt> is not fully initialized.
   *
   * If <tt>return.get()!=NULL</tt>, then it is required that the object
   * referenced by <tt>*return.get()</tt> must have lifetime that extends past
   * the lifetime of the returned smart pointer object.  However, the object
   * referenced by <tt>*return.get()</tt> may change if <tt>*this</tt> is
   * modified so this reference should not be maintained for too long.
   *
   * <b>New Behavior!</b> It is required that the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> must be valid past the lifetime of
   * <tt>*this</tt> vector object.
   */
  virtual RCP< const VectorSpaceBase<Scalar> > space() const = 0;

  //@}

  /** @name Reduction/Transformation operator support */
  //@{

  /** \brief Calls applyOpImpl().
   *
   * Temporary NVI function.
   */
  void applyOp(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const
    {
      applyOpImpl(op, vecs, targ_vecs, reduct_obj, global_offset);
    }

  //@}

  /** @name Vector Cloning */
  //@{

  /** \brief Returns a cloned copy of <tt>*this</tt> vector.
   *
   * This function exists to be consistent with the clone functions
   * <tt>clone()</tt> which creates a <tt>LinearOpBase</tt> object and
   * <tt>clone_mv()</tt> which creates a <tt>MultiVectorBase</tt> object.
   * However, this function is not really necessary because this capability is
   * already present by using the <tt>VectorSpaceBase</tt> returned from
   * <tt>this->space()</tt>.
   *
   * Subclasses should only consider overriding this function if there they
   * want to be very sophisticated and implement some form of lazy evaluation
   * in case the created object might not actually be modified before it is
   * destroyed.  However, this is not advised.
   */
  virtual RCP<VectorBase<Scalar> > clone_v() const = 0;

  //@}

  /** @name Explicit sub-vector access */
  //@{

  /** \brief Calls acquireDetachedVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void acquireDetachedView(
    const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const
    { acquireDetachedVectorViewImpl(rng,sub_vec); }

  /** \brief Calls releaseDetachedVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void releaseDetachedView(
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const
    { releaseDetachedVectorViewImpl(sub_vec); }

  /** \brief Calls acquireNonconstDetachedVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void acquireDetachedView(
    const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec
    )
    { acquireNonconstDetachedVectorViewImpl(rng,sub_vec); }

  /** \brief Calls commitDetachedView().
   *
   * Temporary NVI function.
   */
  void commitDetachedView(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    )
    { commitNonconstDetachedVectorViewImpl(sub_vec); }

  /** \brief Calls setSubVectorImpl().
   *
   * Temporary NVI function.
   */
  void setSubVector(
    const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
    )
    { setSubVectorImpl(sub_vec); }

  //@}

protected:

  /** @name Protected virtual functions to be overridden by subclasses */
  //@{

  /** \brief Apply a reduction/transformation operator over a set of vectors:
   * <tt>op(op(v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) ->
   * z[0]...z[nz-1],(*reduct_obj)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * </ul>
   *
   * The vector <tt>*this</tt> that this function is called on is assumed to
   * be one of the vectors in <tt>v[0]...v[nv-1],z[0]...z[nz-1]</tt>.  This
   * function is generally should not called directly by a client but instead
   * the client should call the nonmember function <tt>Thyra::applyOp()</tt>.
   *
   * See the documentation for the nonmember function <tt>Thyra::applyOp()</tt>
   * for a description of what this function does.
   */
  virtual void applyOpImpl(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const = 0;
  
  /** \brief Get a non-mutable explicit view of a sub-vector.
   *
   * \param rng [in] The range of the elements to extract the sub-vector view.
   *
   * \param sub_vec [in/out] View of the sub-vector.  Prior to the first call
   * to this function, <tt>sub_vec->set_uninitialized()</tt> must be called.
   * Technically <tt>*sub_vec</tt> owns the memory but this memory can be
   * freed only by calling <tt>this->releaseDetachedView(sub_vec)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() < this->space()->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>*sub_vec</tt> contains an explicit non-mutable view to the elements
   *      in the range <tt>full_range(rng,0,this->space()->dim()-1)</tt>
   * </ul>
   *
   * This is only a transient view of a sub-vector that is to be immediately
   * used and then released with a call to <tt>releaseDetachedView()</tt>.
   *
   * Note that calling this function might require some dynamic memory
   * allocations and temporary memory.  Therefore, it is critical that
   * <tt>this->releaseDetachedView(sub_vec)</tt> is called to clean up memory and
   * avoid memory leaks after the sub-vector is used.
   *
   * <b>Heads Up!</b> Note that client code in general should not directly
   * call this function but should instead use the utility class
   * <tt>ConstDetachedVectorView</tt> which will also take care of calling
   * <tt>releaseDetachedView()</tt>.
   *
   * If <tt>this->acquireDetachedView(...,sub_vec)</tt> was previously called on
   * <tt>sub_vec</tt> then it may be possible to reuse this memory if it is
   * sufficiently sized.  The user is encouraged to make multiple calls to
   * <tt>this->acquireDetachedView(...,sub_vec)</tt> before
   * <tt>this->releaseDetachedView(sub_vec)</tt> to finally clean up all of the
   * memory.  Of course, the same <tt>sub_vec</tt> object must be passed to
   * the same vector object for this to work correctly.
   */
  virtual void acquireDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const = 0;

  /** \brief Free an explicit view of a sub-vector.
   *
   * \param sub_vec [in/out] The memory referred to by
   * <tt>sub_vec->values()</tt> will be released if it was allocated and
   * <tt>*sub_vec</tt> will be zeroed out using
   * <tt>sub_vec->set_uninitialized()</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_vec</tt> must have been passed through a call to 
   *      <tt>this->acquireDetachedView(...,sub_vec)</tt>
   * </ul>
    *
   * <b>Postconditions:</b><ul>
   * <li> See <tt>RTOpPack::ConstSubVectorView::set_uninitialized()</tt> for <tt>sub_vec</tt>
   * </ul>
   *
   * The sub-vector view must have been allocated by <tt>this->acquireDetachedView()</tt> first.
   */
  virtual void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const = 0;

  /** \brief Get a mutable explicit view of a sub-vector.
   *
   * \param rng [in] The range of the elements to extract the sub-vector view.
   *
   * \param sub_vec [in/out] Mutable view of the sub-vector.  Prior to the
   * first call to this function <tt>sub_vec->set_uninitialized()</tt> must
   * have been called for the correct behavior.  Technically <tt>*sub_vec</tt>
   * owns the memory but this memory must be committed and freed by calling
   * <tt>this->commitDetachedView(sub_vec)</tt> after the client is finished
   * modifying the view.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!rng.full_range()</tt>] <tt>rng.ubound() < this->space()->dim()</tt>
   *      (throw <tt>std::out_of_range</tt>)
   * </ul>
    *
   * <b>Postconditions:</b><ul>
   * <li> <tt>*sub_vec</tt> contains an explicit mutable view to the elements
   *      in the range <tt>\ref Thyra::full_range() "full_range"(rng,0,this->space()->dim()-1)</tt>
   * </ul>
   *
   * This is only a transient view of a sub-vector that is to be immediately
   * used and then committed back with a call to <tt>commitDetachedView()</tt>.
   *
   * Note that calling this function might require some internal allocations
   * and temporary memory.  Therefore, it is critical that
   * <tt>this->commitDetachedView(sub_vec)</tt> is called to commit the changed
   * entries, clean up memory, and avoid memory leaks after the sub-vector is
   * modified.
   *
   * <b>Heads Up!</b> Note that client code in general should not directly
   * call this function but should instead use the utility class
   * <tt>DetachedVectorView</tt> which will also take care of calling
   * <tt>commitDetachedView()</tt>.
   *
   * If <tt>this->acquireDetachedView(...,sub_vec)</tt> was previously called on
   * <tt>sub_vec</tt> then it may be possible to reuse this memory if it is
   * sufficiently sized.  The user is encouraged to make multiple calls to
   * <tt>this->acquireDetachedView(...,sub_vec)</tt> before
   * <tt>this->commitDetachedView(sub_vec)</tt> is called to finally clean up all
   * of the memory.  Of course the same <tt>sub_vec</tt> object must be passed
   * to the same vector object for this to work correctly.
   *
   * Changes to the underlying sub-vector are not guaranteed to become
   * permanent until <tt>this->acquireDetachedView(...,sub_vec)</tt> is called again,
   * or <tt>this->commitDetachedView(sub_vec)</tt> is called.
   */
  virtual void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec
    ) = 0;

  /** \brief Commit changes for a mutable explicit view of a sub-vector.
   *
   * \param sub_vec [in/out] The data in <tt>sub_vec->values()</tt> will be
   * written back internal storage and the memory referred to by
   * <tt>sub_vec->values()</tt> will be released if it was allocated and
   * <tt>*sub_vec</tt> will be zeroed out using
   * <tt>sub_vec->set_uninitialized()</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_vec</tt> must have been passed through a call to 
   *      <tt>this->acquireDetachedView(...,sub_vec)</tt>
   * </ul>
    *
   * <b>Postconditions:</b><ul>
   * <li> See <tt>RTOpPack::SubVectorView::set_uninitialized()</tt> for <tt>sub_vec</tt>
   * <li> <tt>*this</tt> will be updated according the the changes made to <tt>sub_vec</tt>
   * </ul>
   *
   * The sub-vector view must have been allocated by
   * <tt>this->acquireDetachedView()</tt> first.
   */
  virtual void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    ) = 0;

  /** \brief Set a specific sub-vector.
   *
   * \param sub_vec [in] Represents the elements in the sub-vector to be set.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_vec.globalOffset() + sub_vec.subDim() < this->space()->dim()</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> All of the elements in the range
   *      <tt>[sub_vec.globalOffset(),sub_vec.globalOffset()+sub_vec.subDim()-1]</tt>
   *      in <tt>*this</tt> are set to 0.0 except for those that have that
   *      have entries in <tt>sub_vec</tt> which are set to the values specified
   *      by <tt>(*this)(sub_vec.globalOffset()+vec.localOffset()+sub_vec.indices()[sub_vec.indicesStride()*k])
   *      = vec.values[vec.valueStride()*k]</tt>, for <tt>k = 0..sub_vec.subNz()-1</tt>
   * </ul>
   *
   * After this function returns, the corresponding elements in <tt>*this</tt>
   * vector object will be set equal to those in the input view
   * <tt>sub_vec</tt>.
   */
  virtual void setSubVectorImpl(
    const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
    ) = 0;

  //@}

private:
  
  // Not defined and not to be called
  VectorBase<Scalar>&
  operator=(const VectorBase<Scalar>&);

};


/** \brief Apply a reduction/transformation operator over a set of vectors:
 * <tt>op(op(v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) ->
 * z[0]...z[nz-1],(*reduct_obj)</tt>.
 *
 * \param op [in] Reduction/transformation operator to apply over each
 * sub-vector and assemble the intermediate targets into <tt>reduct_obj</tt>
 * (if <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
 *
 * \param vecs [in] Array (length <tt>num_vecs</tt>) of a set of pointers to
 * non-mutable vectors to include in the operation.  The order of these
 * vectors is significant to <tt>op</tt>.  If <tt>vecs.size()==0</tt>, then
 * <tt>op</tt> is called with no non-mutable vectors.
 *
 * \param targ_vecs [in] Array (length <tt>num_targ_vecs</tt>) of a set of
 * pointers to mutable vectors to include in the operation.  The order of
 * these vectors is significant to <tt>op</tt>.  If
 * <tt>targ_vecs.size()==0</tt>, then <tt>op</tt> is called with no mutable
 * vectors.
 *
 * \param reduct_obj [in/out] Target object of the reduction operation.  This
 * object must have been created by the <tt>op.reduct_obj_create_raw()</tt>
 * function first.  The reduction operation will be added to
 * <tt>*reduct_obj</tt> if <tt>*reduct_obj</tt> has already been through a
 * reduction.  By allowing the info in <tt>*reduct_obj</tt> to be added to the
 * reduction over all of these vectors, the reduction operation can be
 * accumulated over a set of abstract vectors which can be useful for
 * implementing composite vectors for instance.  If
 * <tt>op.get_reduct_type_num_entries(...)</tt> returns <tt>num_values ==
 * 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
 * <tt>reduct_obj</tt> must be set to <tt>null</tt> and no reduction will be
 * performed.
 *
 * \param global_offset [in] (default = 0) The offset applied to the included
 * vector elements.
 *
 * <b>Preconditions:</b><ul>
 *
 * <li> [<tt>vecs.size() > 0</tt>]
 * <tt>vecs[k]->space()->isCompatible(*this->space()) == true</tt>, for <tt>k
 * = 0...vecs.size()-1</tt> (throw
 * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
 *
 * <li> [<tt>targ_vecs.size() > 0</tt>]
 * <tt>targ_vecs[k]->space()->isCompatible(*this->space()) == true</tt>, for
 * <tt>k = 0...targ_vecs.size()-1</tt> (throw
 * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
 *
 * <li> [<tt>targ_vecs.size() > 0</tt>] The vectors pointed to by
 * <tt>targ_vecs[k]</tt>, for <tt>k = 0...targ_vecs.size()-1</tt> must not
 * alias each other or any of the vectors <tt>vecs[k]</tt>, for <tt>k =
 * 0...vecs.size()-1</tt>.  <b>You have be warned!!!!</b>
 *
 * <li> <tt>global_offset >= 0</tt> (throw <tt>std::invalid_argument</tt>)
 *
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 *
 * <li> The vectors in <tt>targ_vecs[]</tt> may be modified as determined by
 * the definition of <tt>op</tt>.
 *
 * <li> [<tt>reduct_obj!=null</tt>] The reduction object <tt>reduct_obj</tt>
 * contains the combined reduction from the input state of <tt>reduct_obj</tt>
 * and the reductions that where accumulated during this this function
 * invocation.
 *
 * </ul>
 *
 * \relates VectorBase
 */
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset = 0
  )
{
  if (vecs.size())
    vecs[0]->applyOp(op, vecs, targ_vecs, reduct_obj, global_offset);
  else if (targ_vecs.size())
    targ_vecs[0]->applyOp(op, vecs, targ_vecs, reduct_obj, global_offset);
}

#ifndef THYRA_HIDE_DEPRECATED_CODE
/** \brief Deprecated.
 *
 * \relates VectorBase
 */
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar> &op,
  const int num_vecs,
  const VectorBase<Scalar>*const vecs_in[],
  const int num_targ_vecs,
  VectorBase<Scalar>*const targ_vecs_inout[],
  RTOpPack::ReductTarget *reduct_obj,
  const Ordinal global_offset = 0
  )
{
  using Teuchos::ptr;
  Array<Ptr<const VectorBase<Scalar> > > vecs(num_vecs);
  for ( int k = 0; k < num_vecs; ++k )
    vecs[k] = ptr(vecs_in[k]);
  Array<Ptr<VectorBase<Scalar> > > targ_vecs(num_targ_vecs);
  for ( int k = 0; k < num_targ_vecs; ++k )
    targ_vecs[k] = ptr(targ_vecs_inout[k]);
  applyOp<Scalar>(op, vecs(), targ_vecs(), ptr(reduct_obj), global_offset);
}
#endif // THYRA_HIDE_DEPRECATED_CODE

} // end namespace Thyra


#endif  // THYRA_VECTOR_BASE_DECL_HPP
