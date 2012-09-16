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

#include "Thyra_VectorBase.hpp"

#ifndef THYRA_EXPLICIT_VECTOR_VIEW_HPP
#define THYRA_EXPLICIT_VECTOR_VIEW_HPP


namespace Thyra {


/** \brief Create an explicit non-mutable (const) view of a <tt>VectorBase</tt> object.
 *
 * This utility class makes it easy to explicitly access a contiguous subset
 * of elements in any <tt>const VectorBase</tt> object.
 *
 * <b>Warning!</b> Creating an explicit view of an arbitrary
 * <tt>%VectorBase</tt> object may be a very expensive operation (such as with
 * distributed-memory vectors) and should only be done in special cases (such
 * as when the vector is an in-core vector).  There several specialized use
 * cases where creating these types of explicit views are necessary but in
 * most cases this should not be done.
 *
 * If one wants to modify the elements in a <tt>VectorBase</tt> object then
 * one should use the utility class <tt>DetachedVectorView</tt>.
 *
 * The default constructor, copy constructor and assignment operators are not
 * allowed.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class ConstDetachedVectorView {
public:

  /** \brief Construct an explicit non-mutable (const) view of a subset of elements.
   *
   * @param v [in] The vector that a view will be taken.  This object must be
   * maintained until <tt>*this</tt> is destroyed.
   *
   * @param rng [in] the range of element indices that the explicit view will
   * be taken.
   *
   * @param forceUnitStride [in] If <tt>true</tt> then the view will have unit
   * stride.
   *
   * Preconditions:<ul>
   *
   * <li>[<tt>rng.full_range()==false</tt>] <tt>rng.ubound() <
   * v.space()->dim()</tt>
   *
   * </ul>
   *
   * Postconditions:<ul>
   *
   * <li><tt>this->sv()</tt> returns the created view
   *
   * <li><tt>this->globalOffset()==rng.lbound()</tt>
   *
   * <li><tt>this->subDim()==rng.size()</tt>
   *
   * <li><tt>this->values()</tt> returns a pointer to a <tt>Scalar</tt> array
   *
   * <li><tt>this->stride()</tt> returns the stride between the elements
   * pointed it in <tt>this->values()</tt>
   *
   * <li>[<tt>forceUnitStride==true</tt>] <tt>this->stride()==1</tt>
   *
   * </ul>
   */
  ConstDetachedVectorView(
    const Teuchos::RCP<const VectorBase<Scalar> > &v
    ,const Range1D &rng = Range1D(), const bool forceUnitStride = false
    )
    {
      this->initialize(v,rng,forceUnitStride);
    }

  /** \brief Construct an explicit non-mutable (const) view of a subset of elements.
   *
   * @param v [in] The vector that a view will be taken.  This object must be
   * maintained until <tt>*this</tt> is destroyed.
   *
   * @param rng [in] the range of element indices that the explicit view will
   * be taken.
   *
   * @param forceUnitStride [in] If <tt>true</tt> then the view will have unit
   * stride.
   *
   * Preconditions:<ul>
   *
   * <li>[<tt>rng.full_range()==false</tt>] <tt>rng.ubound() <
   * v.space()->dim()</tt>
   *
   * </ul>
   *
   * Postconditions:<ul>
   *
   * <li><tt>this->sv()</tt> returns the created view
   *
   * <li><tt>this->globalOffset()==rng.lbound()</tt>
   *
   * <li><tt>this->subDim()==rng.size()</tt>
   *
   * <li><tt>this->values()</tt> returns a pointer to a <tt>Scalar</tt> array
   *
   * <li><tt>this->stride()</tt> returns the stride between the elements
   * pointed it in <tt>this->values()</tt>
   *
   * <li>[<tt>forceUnitStride==true</tt>] <tt>this->stride()==1</tt>
   *
   * </ul>
   */
  ConstDetachedVectorView( const VectorBase<Scalar>& v,
    const Range1D &rng = Range1D(), const bool forceUnitStride = false )
    {
      this->initialize(Teuchos::rcp(&v,false),rng,forceUnitStride);
    }

  /** \brief Free the explicit view on the <tt>VectorBase</tt> object
   * <tt>v</tt> passed to <tt>ConstDetachedVectorView()</tt>.
   */
  ~ConstDetachedVectorView()
    {
      if( sv_s_.stride() != sv_.stride() )
        delete [] const_cast<Scalar*>(sv_.values().get());
      v_->releaseDetachedView(&sv_s_);
    }

  /** \brief Returns the explicit view as an
   * <tt>RTOpPack::ConstSubVectorView<Scalar></tt> object.
   */
  const RTOpPack::ConstSubVectorView<Scalar>& sv() const { return sv_; }

  /** \brief Returns the global offset for the explicit view. */
  Teuchos_Ordinal globalOffset() const { return sv_.globalOffset(); }

  /** \brief Returns the dimension of the explicit view. */
  Teuchos_Ordinal subDim() const { return sv_.subDim(); }

  /** \brief Return a pointer to a <tt>Scalar</tt> array containing the
   * explicit view.
   */
  const Scalar* values() const { return sv_.values().get(); }

  /** \brief Return the stride between elements in the array returned from
   * <tt>this->values()</tt>.
   */
  ptrdiff_t stride() const { return sv_.stride(); }

  /** \brief Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim()-1)</tt>.
   */
  const Scalar& operator[](Teuchos_Ordinal i) const { return sv_[i]; }

  /** \brief Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim()-1)</tt>.
   */
  const Scalar& operator()(Teuchos_Ordinal i) const { return sv_(i); }

private:

  Teuchos::RCP<const VectorBase<Scalar> > v_;
  RTOpPack::ConstSubVectorView<Scalar>  sv_s_;
  RTOpPack::ConstSubVectorView<Scalar>  sv_;

  void initialize(
    const Teuchos::RCP<const VectorBase<Scalar> > &v,
    const Range1D &rng, const bool forceUnitStride
    )
    {
      v_ = v;
      v_->acquireDetachedView(rng,&sv_s_);
      if( forceUnitStride && sv_s_.stride() != 1 ) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "I don't think non-unit stride has ever been tested!");
        //const ArrayRCP<Scalar> values = Teuchos::arcp(sv_s_.subDim());
        //Teuchos_Ordinal i; const Scalar *sv_v;
        //for( sv_v = sv_s_.values().get(), i=0; i < sv_s_.subDim(); ++i, sv_v += sv_s_.stride() )
        //  values[i] = *sv_v;
        //sv_.initialize(sv_s_.globalOffset(),sv_s_.subDim(),values,1);
      }
      else {
        sv_ = sv_s_;
      }
    }
  // Not defined and not to be called
  ConstDetachedVectorView();
  ConstDetachedVectorView(const ConstDetachedVectorView<Scalar>&);
  ConstDetachedVectorView<Scalar>& operator==(const ConstDetachedVectorView<Scalar>&);
};

 
/** \brief Create an explicit mutable (non-const) view of a <tt>VectorBase</tt> object.
 *
 * This utility class makes it easy to explicitly access a contiguous subset
 * of elements in any <tt>VectorBase</tt> object and change the <tt>VectorBase
 * object</tt>.
 *
 * <b>Warning!</b> Creating an explicit view of an arbitrary
 * <tt>%VectorBase</tt> object may be a very expensive operation (such as with
 * distributed-memory vectors) and should only be done in special cases (such
 * as when the vector is an in-core vector).  There several specialized use
 * cases where creating these types of explicit views are necessary but in
 * most cases this should not be done.
 *
 * If one wants to only read the elements in a <tt>VectorBase</tt> object then
 * one should use the utility class <tt>ConstDetachedVectorView</tt>.
 *
 * The default constructor, copy constructor and assignment operators are not
 * allowed.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DetachedVectorView {
public:

  /** \brief Construct an explicit mutable (non-const) view of a subset of elements.
   *
   * @param v [in] The vector that a view will be taken.  This object must be
   * maintained until <tt>*this</tt> is destroyed.  The elements in <tt>v</tt>
   * are not guaranteed to be updated until <tt>*this</tt> is destroyed.
   *
   * @param rng [in] the range of element indices that the explicit view will
   * be taken.
   *
   * @param forceUnitStride [in] If <tt>true</tt> then the view will have unit
   * stride.
   *
   * Preconditions:<ul>
   *
   * <li>[<tt>rng.full_range()==false</tt>] <tt>rng.ubound() <
   * v.space()->dim()</tt>
   *
   * </ul>
   *
   * Postconditions:<ul>
   *
   * <li><tt>this->sv()</tt> returns the created view
   *
   * <li><tt>this->globalOffset()==rng.lbound()</tt>
   *
   * <li><tt>this->subDim()==rng.size()</tt>
   *
   * <li><tt>this->values()</tt> returns a pointer to a <tt>Scalar</tt> array
   *
   * <li><tt>this->stride()</tt> returns the stride between the elements
   * pointed it in <tt>this->values()</tt>
   *
   * <li>[<tt>forceUnitStride==true</tt>] <tt>this->stride()==1</tt>
   *
   * </ul>
   */
  DetachedVectorView(
    const Teuchos::RCP<VectorBase<Scalar> > &v
    ,const Range1D &rng = Range1D(), const bool forceUnitStride = false
    )
    {
      this->initialize(v,rng,forceUnitStride);
    }

  /** \brief Construct an explicit mutable (non-const) view of a subset of elements.
   *
   * @param v [in] The vector that a view will be taken.  This object must be
   * maintained until <tt>*this</tt> is destroyed.  The elements in <tt>v</tt>
   * are not guaranteed to be updated until <tt>*this</tt> is destroyed.
   *
   * @param rng [in] the range of element indices that the explicit view will
   * be taken.
   *
   * @param forceUnitStride [in] If <tt>true</tt> then the view will have unit
   * stride.
   *
   * Preconditions:<ul>
   * <li>[<tt>rng.full_range()==false</tt>] <tt>rng.ubound() < v.space()->dim()</tt>
   * </ul>
   *
   * Postconditions:<ul>
   *
   * <li><tt>this->sv()</tt> returns the created view
   *
   * <li><tt>this->globalOffset()==rng.lbound()</tt>
   *
   * <li><tt>this->subDim()==rng.size()</tt>
   *
   * <li><tt>this->values()</tt> returns a pointer to a <tt>Scalar</tt> array
   *
   * <li><tt>this->stride()</tt> returns the stride between the elements
   * pointed it in <tt>this->values()</tt>
   *
   * <li>[<tt>forceUnitStride==true</tt>] <tt>this->stride()==1</tt>
   *
   * </ul>
   */
  DetachedVectorView( VectorBase<Scalar>& v, const Range1D &rng = Range1D(), const bool forceUnitStride = false )
    {
      this->initialize(Teuchos::rcp(&v,false),rng,forceUnitStride);
    }

  /** \brief Commits back the the explicit view on the <tt>VectorBase</tt>
   * object <tt>v</tt> passed to <tt>DetachedVectorView()</tt>.
   */
  ~DetachedVectorView()
    {
      if( sv_s_.stride() != sv_.stride() ) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "I don't think non-unit stride has ever been tested!");
        //Teuchos_Ordinal i; Scalar *sv_v; const Scalar *values;
        //for (
        //  sv_v = sv_s_.values().get(), values = sv_.values().get(), i=0;
        //  i < sv_s_.subDim();
        //  ++i, sv_v += sv_s_.stride()
        //  )
        //{
        //  *sv_v = *values++;
        //}
        //delete [] sv_.values().get();
      }
      v_->commitDetachedView(&sv_s_);
    }

  /** \brief Returns the explicit view as an
   * <tt>RTOpPack::ConstSubVectorView<Scalar></tt> object.
   */
  const RTOpPack::SubVectorView<Scalar>& sv() const { return sv_; }

  /** \brief Returns the global offset for the explicit view. */
  Teuchos_Ordinal globalOffset() const { return sv_.globalOffset(); }

  /** \brief Returns the dimension of the explicit view. */
  Teuchos_Ordinal subDim() const { return sv_.subDim(); }

  /** \brief Return a pointer to a <tt>Scalar</tt> array containing the
   * explicit view.
   */
  Scalar* values() const { return sv_.values().get(); }

  /** \brief Return the stride between elements in the array returned from
   * <tt>this->values()</tt>.
   */
  ptrdiff_t stride() const { return sv_.stride(); }

  /** \brief Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim()-1)</tt>.
   */
  Scalar& operator[](Teuchos_Ordinal i) const { return sv_[i]; }

  /** \brief Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim()-1)</tt>. */
  Scalar& operator()(Teuchos_Ordinal i) const { return sv_(i); }

private:

  Teuchos::RCP<VectorBase<Scalar> > v_;
  RTOpPack::SubVectorView<Scalar>  sv_s_;
  RTOpPack::SubVectorView<Scalar>  sv_;

  void initialize(
    const Teuchos::RCP<VectorBase<Scalar> > &v
    ,const Range1D &rng, const bool forceUnitStride
    )
    {
      v_ = v;
      v_->acquireDetachedView(rng,&sv_s_);
      if( forceUnitStride && sv_s_.stride() != 1 ) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "I don't think non-unit stride has ever been tested!");
        //Scalar *values = new Scalar[sv_s_.subDim()];
        //Teuchos_Ordinal i; const Scalar *sv_v;
        //for( sv_v = sv_s_.values().get(), i=0; i < sv_s_.subDim(); ++i, sv_v += sv_s_.stride() )
        //  values[i] = *sv_v;
        //sv_.initialize(sv_s_.globalOffset(),sv_s_.subDim(),values,1);
      }
      else {
        sv_ = sv_s_;
      }
    }

  // Not defined and not to be called
  DetachedVectorView();
  DetachedVectorView(const DetachedVectorView<Scalar>&);
  DetachedVectorView<Scalar>& operator==(const DetachedVectorView<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_EXPLICIT_VECTOR_VIEW_HPP
