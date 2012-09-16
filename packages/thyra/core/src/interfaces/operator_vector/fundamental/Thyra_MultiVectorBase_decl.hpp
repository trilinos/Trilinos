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

#ifndef THYRA_MULTI_VECTOR_BASE_DECL_HPP
#define THYRA_MULTI_VECTOR_BASE_DECL_HPP

#include "Thyra_LinearOpBase_decl.hpp"
#include "RTOpPack_RTOpT.hpp"


namespace Thyra {


/** \brief Interface for a collection of column vectors called a multi-vector.
 * 
 * \section Thyra_MVB_outline_sec Outline
 * 
 * <ul>
 * <li>\ref Thyra_MVB_intro_sec
 * <li>\ref Thyra_MVB_views_sec
 *   <ul>
 *   <li>\ref Thyra_MVB_col_access_sec
 *   <li>\ref Thyra_MVB_subviews_sec
 *   <li>\ref Thyra_MVB_view_behavior_sec
 *   </ul>
 * <li>\ref Thyra_MVB_as_LO_sec
 *   <ul>
 *   <li>\ref Thyra_MVB_block_update_sec
 *   <li>\ref Thyra_MVB_block_inner_prod_sec
 *   </ul>
 * <li>\ref Thyra_MVB_RTOp_sec
 * <li>\ref Thyra_MVB_rtop_collection_sec
 * <li>\ref Thyra_MVB_expl_access_sec
 * <li>\ref Thyra_MVB_expl_access_utils_sec
 * <li>\ref Thyra_MVB_dev_notes_sec
 * </ul>
 * 
 * \section Thyra_MVB_intro_sec Introduction
 * 
 * The primary purpose for this interface is to allow for the convenient
 * aggregation of column vectors as a single matrix-like object.  Such an
 * orderly arrangement of column vectors into a single aggregate object allows
 * for better optimized linear algebra operations such as matrix-matrix
 * multiplication and the solution of linear systems for multiple right-hand
 * sides.  Every computing environment (serial, parallel, out-of-core etc.) 
 * should be able to define at least one reasonably efficient implementation
 * of this interface.
 *
 * \section Thyra_MVB_views_sec Changeable and non-changeable views
 *
 * This interface allows a client to create <tt>VectorBase</tt> and
 * <tt>MultiVectorBase</tt> views of single and multiple columns of
 * <tt>*this</tt> multi-vector, respectively.  These two views are described
 * separately in the next two subsections.  The behavior of the vector and
 * multi-vector views is very similar and the common behavior is described in
 * the subsection \ref Thyra_MVB_view_behavior_sec.
 * 
 * \subsection Thyra_MVB_col_access_sec Accessing the individual columns as vector views
 * 
 * The individual columns of a multi-vector can be access using the non-const
 * and const versions of the <tt>col()</tt> function.  For example, the
 * individual columns of one multi-vector can be copied to another as follows

 \code

  template<class Scalar>
  void copyOneColumnAtATime(
    const Thyra::MultiVectorBase<Scalar>   &X
    ,Thyra::MultiVectorBase<Scalar>        *Y
    )
  {
    for( int j = =; j < X.domain()->dim(); ++j )
      assign( &*Y->col(j), *X.col(j) );
  } 

 \endcode

 * In the above code fragment, the expression <tt>X.col(j)</tt> returns a
 * smart-pointer to a non-changeable column of <tt>X</tt> while the expression
 * <tt>Y->col(j)</tt> returns a smart-pointer to a changeable column of
 * <tt>Y</tt> which is being modified.
 *
 * <b>Note:</b> A modification to <tt>Y</tt> is not guaranteed to be committed
 * back to <tt>Y</tt> until the smart pointer returned from <tt>Y.col(j)</tt>
 * is deleted.
 * 
 * \subsection Thyra_MVB_subviews_sec Accessing collections of columns as multi-vector views
 * 
 * Another important aspect of this interface is the ability to allow clients
 * to access non-changeable and changeable <tt>MultiVectorBase</tt> views of
 * columns of a parent <tt>MultiVectorBase</tt> object.  These sub-views are
 * created with one of the overloaded <tt>subView()</tt> functions of which
 * there are two general forms.
 * 
 * The first form provides views of contiguous columns through the functions
 * <tt>subView(const Range1D&)</tt> and <tt>subView(const Range1D&)const</tt>.
 * For example, the following function shows how to copy the first three columns
 * of one multi-vector to the last three columns of another multi-vector.

 \code

 template<class Scalar>
 void copyThreeColumns(
   const Thyra::MultiVectorBase<Scalar> &X
   ,Thyra::MultiVectorBase<Scalar>      *Y
   )
 {
   const int m = Y->domain()->dim();
   assign( &*Y->subView(Range1D(m-3,m-1)), *X.subView(Range1D(0,2)) );
 }

 \endcode
 
 * <b>Note:</b> In the above example <tt>*Y</tt> can be the same multi-vector
 * as <tt>X</tt>.
 *
 * <b>Note:</b> In the above example <tt>*Y</tt> is not guaranteed to be
 * updated until the view returned from <tt>Y->subView(Range1D(m-3,m-1)</tt> is
 * destroyed (which occurs at the end of the statement in which it occurs in
 * this case).
 *
 * <!-- Warning! Do not reformat the below paragraph or the \ref links will break! --> 
 * The second form provides views of non-contiguous columns through
 * the functions
 * <tt>\ref Thyra_MVB_subView_noncontiguous_nonconst "subView(const int numCols, const int cols[])"</tt>
 * and
 * <tt>\ref Thyra_MVB_subView_noncontiguous_const "subView(const int numCols, const int cols[]) const"</tt>.
 * For example, the following function copies columns 1,
 * 3, and 5 from one multi-vector to columns 2, 4, and 6 of another
 * multi-vector.

 \code

 template<class Scalar>
 void copyThreeStaggeredColumns(
   const Thyra::MultiVectorBase<Scalar> &X
   ,Thyra::MultiVectorBase<Scalar>      *Y
   )
 {
   using Teuchos::tuple;
   assign( Y->subView(tuple<int>(2,4,6)()).ptr(), *X.subView(tuple<int>(1,3,5)()) );
 }

 \endcode

 * <b>Note:</b> In the above example <tt>*Y</tt> can be the same multi-vector
 * as <tt>X</tt>.
 *
 * <b>Note:</b> In the above example <tt>*Y</tt> is not guaranteed to be
 * updated until the view returned from
 * <tt>Y->subView(tuple<int>(2,4,6)())</tt> is destroyed (which occurs at
 * the end of the statement in which it occurs in this case).
 * 
 * In general, the first contiguous form of views will be more efficient that
 * the second non-contiguous form.  Therefore, user's should try to structure
 * their ANAs to use the contiguous form of multi-vector views if possible and
 * only result to the non-contiguous form of views when absolutely needed.
 *
 * \subsection Thyra_MVB_view_behavior_sec Common behavior of vector and multi-vector views
 *
 * When a view is created it may become a largely separate object from the
 * parent multi-vector and the exact relationship between the two in undefined
 * by this interface.  This is true whether we are talking about individual
 * column vector views or contiguous or non-contiguous multiple-column
 * multi-vector views which are described above.  These views and the parent
 * multivector follow the state behavior outlined \ref
 * Thyra_Op_Vec_Behavior_Of_Views_grp "here".
 *
 * If <tt>X_view</tt> is some view of a parent multi-vector <tt>X</tt> the
 * following restrictions apply:
 *
 * <ul>
 *
 * <li><b>Undefined behavior:</b> Changing the parent
 *
 * The client should not attempt to change the parent multi-vector <tt>X</tt>
 * while any view is active.  The behavior of doing so is undefined.  For
 * example, the value returned from the following function and the final state
 * of the parent multi-vector <tt>*X</tt> are undefined:

 \code

  template<class Scalar>
  Scalar undefinedBehaviorFromChangingParent(
    Thyra::MultiVectorBase<Scalar>   *X
    )
  {
    // Create the view
    RCP< Thyra::MultiVectorBase<Scalar> >
      X_view = X->subView(Teuchos::Range1D(0,0));
    // Change the parent while the view is still active
    Teuchos::assign( X, Teuchos::ScalarTraits<Scalar>::one() );
    // Above, changing the parent multi-vector may or may not change the subview
    return Teuchos::norm_1(*X_view); // The value returned is undefined
    // When the RCP X_view goes out of scope here, the state of the 
    // parent multi-vector *X is undefined!
  } 

 \endcode

 * <li><b>Undefined behavior:</b> Changing the view and accessing the parent while
 * view is still active
 *
 * The client should not attempt to change the view <tt>X_view</tt> and then
 * access the parent multi-vector <tt>X</tt> while the view is still active.
 * The behavior of doing so is undefined.  For example, the value returned
 * from the following function is undefined:

 \code

  template<class Scalar>
  Scalar undefinedBehaviorFromChaningViewAndAccessingParent(
    Thyra::MultiVectorBase<Scalar>   *X
    )
  {
    // Create the view
    RCP< Thyra::MultiVectorBase<Scalar> >
      X_view = X->subView(Teuchos::Range1D(0,0));
    // Change the view
    Teuchos::assign( *&X_view, Teuchos::ScalarTraits<Scalar>::one() );
    // Above, changing the view may or may not immediately update the parent multi-vector
    return Teuchos::norm_1(*X); // The value returned from the parent is undefined at this point
    // When the RCP X_view goes out of scope here, the parent multi-vector
    // *X is guaranteed to be updated 
  } 

 \endcode

 * <li><b>Undefined behavior:</b> Creating overlapping changeable views
 *
 * The client should not attempt to create overlapping changeable views.  If
 * any of these changeable views is modified, the the behavior of the parent
 * multi-vector is undefined.  For example, the state of the parent
 * multi-vector <tt>*X</tt> is undefined after the following function returns:

 \code

  template<class Scalar>
  Scalar undefinedBehaviorOfOverlappingViews(
    Thyra::MultiVectorBase<Scalar>   *X
    )
  {
    // Create two overlapping views
    RCP< Thyra::MultiVectorBase<Scalar> >
      X_view1 = X->subView(Teuchos::Range1D(0,0)),
      X_view2 = X->subView(Teuchos::Range1D(0,0));
    // Change one of the views but not the other
    Teuchos::assign( *&X_view2, Teuchos::ScalarTraits<Scalar>::one() );
    // Once the RCPs X_view1 and X_view2 go out of scope here,
    // the state of the parent multi-vector *X is undefined!  In some cases,
    // the intial view in X_view1 will be relected in X and in other
    // cases the view in X_view2 will be written to the parent X.
  } 

 \endcode
 
 * Note that overlapping non-changeable views of a multi-vector are just fine
 * since they do not change the state of the parent multi-vector.
 *
 * </ul>
 *
 * In general, to stay out of trouble with multi-vector views:
 *
 * <ul>
 *
 * <li>Never change or access the parent multi-vector while a changeable view
 * is active.
 *
 * <li>Never create simultaneous changeable overlapping views.
 *
 * </ul>
 *
 * Note, however, that creating simultaneous non-overlapping non-changeable or
 * changeable views is just fine as long as the parent multi-vector is not
 * modified while the views are active.  For example, the final state of
 * <tt>*X</tt> is well defined after the following function finishes
 * executing:

 \code

  template<class Scalar>
  Scalar wellDefinedBehaviorOfNonOverlappingViews(
    Thyra::MultiVectorBase<Scalar>   *X
    )
  {
    // Create two non-overlapping views
    RCP< Thyra::MultiVectorBase<Scalar> >
      X_view1 = X->subView(Teuchos::Range1D(0,0)),
      X_view2 = X->subView(Teuchos::Range1D(1,1));
    // Change the two views
    Teuchos::assign( *&X_view1, Teuchos::ScalarTraits<Scalar>::zero() );
    Teuchos::assign( *&X_view2, Teuchos::ScalarTraits<Scalar>::one() );
    // When the RCPs X_view1 and X_view2 go out of scope here,
    // the state of the parent multi-vector *X will be guaranteed to be
    // updated to the values changed in these views.
  } 

 \endcode

 * \section Thyra_MVB_as_LO_sec MultiVectorBase as a linear operator
 * 
 * The <tt>%MultiVectorBase</tt> interface is derived from the
 * <tt>LinearOpBase</tt> interface and therefore every
 * <tt>%MultiVectorBase</tt> object can be used as a linear operator which has
 * some interesting implications.  Since a linear operator can apply itself to
 * vectors and multi-vectors and a multi-vector is a linear operator, this
 * means that a multi-vector can apply itself to other vectors and
 * multi-vectors.  There are several different use cases that this
 * functionality is useful.  Two of the more important use cases are block
 * updates and block inner products.
 *
 * \subsection Thyra_MVB_block_update_sec Multi-vector block updates
 *
 * Let <tt>V</tt> and <tt>Y</tt> be multi-vector objects with the same vector
 * space with a very large number of rows <tt>m</tt> and a moderate number of
 * columns <tt>n</tt>.  Now, consider the block update of the form
 
 \verbatim

   Y = Y + V * B
 \endverbatim

 * where the multi-vector <tt>B</tt> is of dimension <tt>n x b</tt>.
 *
 * The following function shows how this block update might be performed.

 \code

 template<class Scalar>
 void myBlockUpdate(
   const Thyra::MultiVectorBase<Scalar> &V
   ,const int                           b
   ,Thyra::MultiVectorBase<Scalar>      *Y
   )
 {
   typedef Teuchos::ScalarTraits<Scalar> ST;
   // Create the multi-vector B used for the update
   RCP<Thyra::MultiVectorBase<Scalar> >
     B = createMembers(V.domain(),b);
   // Fill up B for the update
   ...
   // Do the update Y = V*B + Y
   V.apply(Thyra::NONCONJ_ELE,*B,Y,ST::one(),ST::one());
 }

 \endcode

 * In a block update, as demonstrated above, level-3 BLAS can be used to
 * provide a very high level of performance.  Note that in an SPMD program,
 * that <tt>B</tt> would be a locally replicated multi-vector and <tt>V</tt>
 * and <tt>Y</tt> would be distributed-memory multi-vectors.  In an SPMD
 * environment, there would be no global communication in a block update.
 *
 * \subsection Thyra_MVB_block_inner_prod_sec Multi-vector block inner products
 *
 * An important operation supported by the
 * <tt>LinearOpBase::applyTranspose()</tt> function is the block inner product
 * which takes the form
 
 \verbatim

   B = adjoint(V)*X
 \endverbatim

 * where <tt>V</tt> and <tt>X</tt> are tall, thin multi-vectors and <tt>B</tt>
 * is a small multi-vector.  In an SPMD environment, <tt>V</tt> and <tt>X</tt>
 * would be distributed-memory objects while <tt>B</tt> would be locally
 * replicated in each process.  The following function shows how block inner
 * product would be performed:

 \code

 template<class Scalar>
 RCP<Thyra::MultiVectorBase<Scalar> >
 void myBlockInnerProd(
   const Thyra::MultiVectorBase<Scalar>                    &V
   ,const Thyra::MultiVectorBase<Scalar>                   &X
   )
 {
   // Create the multi-vector B used for the result
   RCP<Thyra::MultiVectorBase<Scalar> >
     B = createMembers(V.domain(),X.domain()->dim());
   // Do the inner product B = adjoint(V)*X
   V.applyTranspose(Thyra::CONJ_ELE,X,&*B);
   // Return the result
   return B;
 }

 \endcode

 * In an SPMD program, the above block inner product will use level-3 BLAS to
 * multiply the local elements of <tt>V</tt> and <tt>X</tt> and will then do a
 * single global reduction to assemble the product <tt>B</tt> in all of the
 * processes.
 *
 * \section Thyra_MVB_RTOp_sec Support for reduction/transformation operations
 * 
 * Another powerful feature of this interface is the ability to apply
 * reduction/transformation operators over a sub-set of rows and columns in a
 * set of multi-vector objects using the <tt>applyOp()</tt> functions.  The
 * behavior is identical to the client extracting each column in a set of
 * multi-vectors and calling <tt>VectorBase::applyOp()</tt> individually on
 * these columns.  However, the advantage of using the multi-vector apply
 * functions is that there may be greater opportunities for increased
 * performance in a number of respects.  Also, the intermediate reduction
 * objects over a set of columns can be reduced by a secondary reduction
 * object.
 * 
 * \section Thyra_MVB_rtop_collection_sec Collection of pre-written RTOps and wrapper functions
 *
 * There already exists RTOp-based implementations of several standard vector
 * operations and some convenience functions that wrap these operators and
 * call <tt>applyOp()</tt>.  See the Operator/Vector Support Software
 * collection for these.
 * 
 * \section Thyra_MVB_expl_access_sec Explicit multi-vector coefficient access
 * 
 * This interface also allows a client to extract a sub-set of elements in an
 * explicit form as non-changeable <tt>RTOpPack::ConstSubMultiVectorView</tt> objects or
 * changeable <tt>RTOpPack::SubMultiVectorView</tt> objects using the
 * <tt>acquireDetachedView()</tt> functions.  In general, this is a very bad
 * thing to do and should be avoided.  However, there are some situations
 * where this is needed, just as is the case for vectors (see \ref
 * Thyra_VB_expl_access_sec).  The default implementation of these explicit
 * access functions use sophisticated reduction/transformation operators with
 * the <tt>applyOp()</tt> function in order to extract and set the needed
 * elements.  Therefore, all <tt>%MultiVectorBase</tt> subclasses
 * automatically support these operations (even if it is a bad idea to use
 * them).
 * 
 * \section Thyra_MVB_expl_access_utils_sec Explicit multi-vector coefficient access utilities
 * 
 * Client code in general should not directly call the above described
 * explicit sub-multi-vector access functions but should instead use the
 * utility classes <tt>ConstDetachedMultiVectorView</tt> and
 * <tt>DetachedMultiVectorView</tt> since these are easier to use and safer in
 * the event that an exception is thrown.  These classes are documented in the
 * Operator/Vector Support Software collection.
 * 
 * \section Thyra_MVB_dev_notes_sec Notes for subclass developers
 * 
 * This is a fairly bare-bones interface class without much in the way of
 * default function implementations.  The subclass
 * <tt>MultiVectorDefaultBase</tt> (contained in the Operator/Vector Support
 * Software collection) uses a default multi-vector implementation to provide
 * overrides of many of the functions and should be the first choice for
 * subclasses implementations to derive their implementations from rather than
 * starting from scratch.
 * 
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class MultiVectorBase : virtual public LinearOpBase<Scalar>
{
public:

#ifdef THYRA_INJECT_USING_DECLARATIONS
  using LinearOpBase<Scalar>::apply;
#endif

  /** @name Provide access to the columns as VectorBase objects */
  //@{

  /** \brief Calls colImpl().
   *
   * Temporary NVI function.
   */
  RCP<const VectorBase<Scalar> > col(Ordinal j) const
    { return colImpl(j); }

  /** \brief Calls nonconstColImpl().
   *
   * Temporary NVI function.
   */
  RCP<VectorBase<Scalar> > col(Ordinal j)
    { return nonconstColImpl(j); }

  //@}

  /** @name Multi-vector sub-views */
  //@{

  /** \brief Calls contigSubViewImpl().
   *
   * Temporary NVI function.
   */
  RCP<const MultiVectorBase<Scalar> >
  subView( const Range1D& colRng ) const
    {
      return contigSubViewImpl(colRng);
    }

  /** \brief Calls nonconstContigSubViewImpl().
   *
   * Temporary NVI function.
   */
  RCP<MultiVectorBase<Scalar> >
  subView( const Range1D& colRng )
    { return nonconstContigSubViewImpl(colRng); }

  /** \brief nonContigSubViewImpl().
   *
   * Temporary NVI function.
   */
  RCP<const MultiVectorBase<Scalar> >
  subView( const ArrayView<const int> &cols ) const
    { return nonContigSubViewImpl(cols); }

  /** \brief nonconstNonContigSubViewImpl().
   *
   * Temporary NVI function.
   */
  RCP<MultiVectorBase<Scalar> >
  subView( const ArrayView<const int> &cols )
    { return nonconstNonContigSubViewImpl(cols); }
  
  //@}

  /** @name Collective reduction/transformation operator apply functions */
  //@{

  /** \brief Calls mvMultiReductApplyOpImpl().
   *
   * Temporary NVI function.
   */
  void applyOp(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
    const Ordinal primary_global_offset
    ) const
    {
      mvMultiReductApplyOpImpl(primary_op, multi_vecs, targ_multi_vecs,
        reduct_objs, primary_global_offset);
    }

  /** \brief mvSingleReductApplyOpImpl().
   *
   * Temporary NVI function.
   */
  void applyOp(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const RTOpPack::RTOpT<Scalar> &secondary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal primary_global_offset
    ) const
    {
      mvSingleReductApplyOpImpl(primary_op, secondary_op, multi_vecs, targ_multi_vecs,
        reduct_obj, primary_global_offset);
    }
  
  //@}

  /** @name Explicit sub-multi-vector access */
  //@{

  /** \brief Calls acquireDetachedMultiVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void acquireDetachedView(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
    ) const
    { acquireDetachedMultiVectorViewImpl( rowRng, colRng, sub_mv ); }

  /** \brief Calls releaseDetachedMultiVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void releaseDetachedView(
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const
    { releaseDetachedMultiVectorViewImpl(sub_mv); }

  /** \brief Calls acquireNonconstDetachedMultiVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void acquireDetachedView(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar> *sub_mv
    )
    { acquireNonconstDetachedMultiVectorViewImpl(rowRng,colRng,sub_mv); }

  /** \brief Calls commitNonconstDetachedMultiVectorViewImpl().
   *
   * Temporary NVI function.
   */
  void commitDetachedView(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    )
    { commitNonconstDetachedMultiVectorViewImpl(sub_mv); }

  //@}

  /** @name Cloning */
  //@{

  /** \brief Clone the multi-vector object (if supported).
   *
   * The default implementation uses the vector space to create a
   * new multi-vector object and then uses a transformation operator
   * to assign the vector elements.  A subclass should only override
   * this function if it can do something more sophisticated
   * (i.e. lazy evaluation) but in general, this is not needed.
   */
  virtual RCP<MultiVectorBase<Scalar> > clone_mv() const = 0;

  //@}

  /** @name Overridden functions from LinearOpBase */
  //@{

  /// This function is simply overridden to return <tt>this->clone_mv()</tt>.
  RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

protected:

  /** @name Protected virtual functions to be overridden by subclasses */
  //@{

  /** \brief Return a non-changeable view of a constituent column vector.
   *
   * \param j [in] zero-based index of the column to return a view for
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>0 <= j && j < this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>this->range()->isCompatible(*return->space()) == true</tt>
   * </ul>
   *
   * See \ref Thyra_MVB_col_access_sec and \ref Thyra_MVB_view_behavior_sec
   * for the behavior of this view.
   *
   * The default implementation of this function (which is the only
   * implementation needed by most subclasses) is based on the
   * non-const version <tt>col()</tt>.
   */
  virtual RCP<const VectorBase<Scalar> > colImpl(Ordinal j) const;

  /** \brief Return a changeable view of a constituent column vector.
   *
   * \param j [in] zero-based index of the column to return a view for
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>0 <= j && j < this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>this->range()->isCompatible(*return->space()) == true</tt>
   * </ul>
   *
   * See \ref Thyra_MVB_col_access_sec and \ref Thyra_MVB_view_behavior_sec
   * for the behavior of this view.
   *
   * <b>Note:</b> <tt>*this</tt> is not guaranteed to be modified until the
   * smart pointer returned by this function is destroyed.
   */
  virtual RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j) = 0;

  /** \brief Return a non-changeable sub-view of a contiguous set of columns
   * of the this multi-vector.
   *
   * \anchor Thyra_MVB_subView_contiguous_const
   *
   * \param colRng [in] zero-based range of columns to create a view of.  Note
   * that it is valid for <tt>colRng.full_range()==true</tt> in which case the
   * view of the entire multi-vector is taken.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>colRng.ubound() < this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == Teuchos::full_range(colRng,0,this->domain()->dim()-1).size()</tt>
   * <li> <tt>*return->col(k)</tt> represents the same column vector as <tt>this->col(colRng.lbound()+k)</tt>,
   *      for <tt>k=0...Teuchos::full_range(colRng,0,this->domain()->dim()).ubound()-1</tt>
   * </ul>
   *
   * See \ref Thyra_MVB_subviews_sec and \ref Thyra_MVB_view_behavior_sec for
   * the behavior of this view.
   */
  virtual RCP<const MultiVectorBase<Scalar> >
  contigSubViewImpl( const Range1D& colRng ) const = 0;

  /** \brief Return a changeable sub-view of a contiguous set of columns of
   * the this multi-vector.
   *
   * \anchor Thyra_MVB_subView_contiguous_nonconst
   *
   * \param colRng [in] zero-based range of columns to create a view of.  Note
   * that it is valid for <tt>colRng.full_range()==true</tt> in which case the
   * view of the entire multi-vector is taken.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>colRng.ubound() < this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == Teuchos::full_range(colRng,0,this->domain()->dim()-1).size()</tt>
   * <li> <tt>*return->col(k)</tt> represents the same column vector as <tt>this->col(colRng.lbound()+k)</tt>,
   *      for <tt>k=0...Teuchos::full_range(colRng,0,this->domain()->dim()).ubound()-1</tt>
   * </ul>
   *
   * See \ref Thyra_MVB_subviews_sec and \ref Thyra_MVB_view_behavior_sec for
   * the behavior of this view.
   */
  virtual RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl( const Range1D& colRng ) = 0;

  /** \brief Return a non-changeable sub-view of a non-contiguous set of columns of this multi-vector.
   *
   * \anchor Thyra_MVB_subView_noncontiguous_const
   *
   * \param cols [in] Array of the zero-based column indexes to use in the
   * returned view.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>cols.size() <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>0 <= cols[k] < this->domain()->dim()</tt>, for <tt>k=0...cols.size()-1</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>col[k1] != col[k2]</tt>, for all <tt>k1 != k2</tt> in the range <tt>[0,cols.size()-1]</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == cols.size()</tt>
   * <li> <tt>*return->col(k)</tt> represents the same column vector as <tt>this->col(cols[k])</tt>,
   *      for <tt>k=0...cols.size()-1</tt>
   * </ul>
   *
   * See \ref Thyra_MVB_subviews_sec and \ref Thyra_MVB_view_behavior_sec for
   * the behavior of this view.
   */
  virtual RCP<const MultiVectorBase<Scalar> >
  nonContigSubViewImpl( const ArrayView<const int> &cols ) const = 0;

  /** \brief Return a changeable sub-view of a non-contiguous set of columns of this multi-vector.
   *
   * \anchor Thyra_MVB_subView_noncontiguous_nonconst
   *
   * \param cols [in] Array of the zero-based column indexes to use in the
   * returned view.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>cols.size() <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>0 <= cols[k] < this->domain()->dim()</tt>, for <tt>k=0...cols.size()-1</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>col[k1] != col[k2]</tt>, for all <tt>k1 != k2</tt> in the range <tt>[0,cols.size()-1]</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == cols.size()</tt>
   * <li> <tt>*return->col(k)</tt> represents the same column vector as <tt>this->col(cols[k])</tt>,
   *      for <tt>k=0...cols.size()-1</tt>
   * </ul>
   *
   * See \ref Thyra_MVB_subviews_sec and \ref Thyra_MVB_view_behavior_sec for
   * the behavior of this view.
   */
  virtual RCP<MultiVectorBase<Scalar> >
  nonconstNonContigSubViewImpl( const ArrayView<const int> &cols ) = 0;

  /** \brief Apply a reduction/transformation operator column by column and
   * return an array of the reduction objects.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> See the preconditions for <tt>Thyra::applyOp()</tt>
   * </ul>
   *
   * See the documentation for the function <tt>Thyra::applyOp()</tt>
   * for a description of the arguments.
   *
   * This function is not to be called directly by the client but instead
   * through the nonmember function <tt>Thyra::applyOp()</tt>.
   *
   * It is expected that <tt>this</tt> will be one of the multi-vector
   * objects in <tt>multi_vecs[]</tt> or <tt>targ_multi_vecs[]</tt>.
   *
   * The default implementation calls <tt>VectorBase::applyOp()</tt> on
   * each column <tt>this->col(j)</tt> for <tt>j = 0
   * ... this->range()->dim()-1</tt>.
   */
  virtual void mvMultiReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
    const Ordinal primary_global_offset
    ) const = 0;

  /** \brief Apply a reduction/transformation operator column by column and
   * reduce the intermediate reduction objects into a single reduction object.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> See the preconditions for <tt>Thyra::applyOp()</tt>
   * </ul>
   *
   * See the documentation for the function <tt>Thyra::applyOp()</tt>
   * for a description of the arguments.
   *
   * This function is not to be called directly by the client but instead
   * through the nonmember function <tt>Thyra::applyOp()</tt>.
   *
   * It is expected that <tt>this</tt> will be one of the multi-vector
   * objects in <tt>multi_vecs[]</tt> or <tt>targ_multi_vecs[]</tt>.
   *
   * The default implementation calls <tt>applyOp()</tt> where an
   * array of reduction objects is taken.
   */
  virtual void mvSingleReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const RTOpPack::RTOpT<Scalar> &secondary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal primary_global_offset
    ) const = 0;

  /** \brief Get a non-changeable explicit view of a sub-multi-vector.
   *
   * \param rowRng [in] The range of the rows to extract the sub-multi-vector
   * view.
   *
   * \param colRng [in] The range of the columns to extract the
   * sub-multi-vector view.
   *
   * \param sub_mv [in/out] View of the sub-multi_vector.  Prior to the first
   * call to this function, <tt>sub_mv->set_uninitialized()</tt> must be
   * called.  Technically <tt>*sub_mv</tt> owns the memory but this memory can
   * be freed only by calling <tt>this->releaseDetachedView(sub_mv)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->range().get()!=NULL && this->domain().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!rowRng.full_range()</tt>] <tt>rowRng.ubound() < this->range()->dim()</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>colRng.ubound() < this->domain()->dim()</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>*sub_mv</tt> contains an explicit non-changeable view to the elements
   *      in the row and column ranges <tt>Teuchos::full_range(rowRng,0,this->range()->dim()-1)</tt>
   *      and <tt>Teuchos::full_range(colRng,0,this->domain()->dim()-1)</tt> respectively.
   * </ul>
   *
   * <b>Note:</b> This view is to be used immediately and then released with a
   * call to <tt>releaseDetachedView()</tt>.
   *
   * Note that calling this operation might require some dynamic memory
   * allocations and temporary memory.  Therefore, it is critical that
   * <tt>this->releaseDetachedView(sub_mv)</tt> be called by client in order to
   * clean up memory and avoid memory leaks after the sub-multi-vector view is
   * finished being used.
   *
   * <b>Heads Up!</b> Note that client code in general should not directly
   * call this function but should instead use the utility class
   * <tt>ConstDetachedMultiVectorView</tt> which will also take care of calling
   * <tt>releaseDetachedView()</tt>.
   *
   * If <tt>this->acquireDetachedView(...,sub_mv)</tt> was previously
   * called on <tt>sub_mv</tt> then it may be possible to reuse this
   * memory if it is sufficiently sized.  The user is encouraged to
   * make multiple calls to
   * <tt>this->acquireDetachedView(...,sub_mv)</tt> before
   * <tt>this->releaseDetachedView(sub_mv)</tt> to finally clean up all of
   * the memory.  Of course, the same <tt>sub_mv</tt> object must be
   * passed to the same multi-vector object for this to work correctly.
   *
   * This function has a default implementation based on the vector operation
   * <tt>VectorBase::acquireDetachedView()</tt> called on the non-changeable vector
   * objects returned from <tt>col()</tt>.  Note that the footprint of the
   * reduction object (both internal and external state) will be
   * O(<tt>rowRng.size()*colRng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, however, this is a terrible implementation and
   * must be overridden if <tt>rowRng.size()</tt> is large at all.  Although,
   * this function should not even be used in cases where the multi-vector is
   * very large.  If a subclass does override this function, it must also
   * override <tt>releaseDetachedView()</tt> which has a default implementation
   * which is a companion to this function's default implementation.
   */
  virtual void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
    ) const = 0;

  /** \brief Free a non-changeable explicit view of a sub-multi-vector.
   *
   * \param sub_mv * [in/out] The memory referred to by
   * <tt>sub_mv->values()</tt> * will be released if it was allocated and
   * <tt>*sub_mv</tt> will be zeroed out using
   * <tt>sub_mv->set_uninitialized()</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->range().get()!=NULL && this->domain().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_mv</tt> must have been passed through a call to 
   *      <tt>this->acquireDetachedView(...,sub_mv)</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> See <tt>RTOpPack::ConstSubMultiVectorView::set_uninitialized()</tt> for <tt>sub_mv</tt>
   * </ul>
   *
   * The sub-multi-vector view must have been allocated by
   * <tt>this->acquireDetachedView()</tt> first.
   *
   * This function has a default implementation which is a companion
   * to the default implementation for <tt>acquireDetachedView()</tt>.  If
   * <tt>acquireDetachedView()</tt> is overridden by a subclass then this
   * function must be overridden also!
   */
  virtual void releaseDetachedMultiVectorViewImpl(
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const = 0;

  /** \brief Get a changeable explicit view of a sub-multi-vector.
   *
   * \param rowRng [in] The range of the rows to extract the sub-multi-vector
   * view.
   *
   * \param colRng [in] The range of the columns to extract the
   * sub-multi-vector view.
   *
   * \param sub_mv [in/out] Changeable view of the sub-multi-vector.  Prior to
   * the first call <tt>sub_mv->set_uninitialized()</tt> must have been called
   * for the correct behavior.  Technically <tt>*sub_mv</tt> owns the memory
   * but this memory must be committed and freed only by calling
   * <tt>this->commitDetachedView(sub_mv)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->range().get()!=NULL && this->domain().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!rowRng.full_range()</tt>] <tt>rowRng.ubound() < this->range()->dim()</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>colRng.ubound() < this->domain()->dim()</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
    *
   * <b>Postconditions:</b><ul>
   * <li> <tt>*sub_mv</tt> contains an explicit changeable view to the elements
   *      in the row and column ranges <tt>full_range(rowRng,0,this->range()->dim()-1)</tt>
   *      and <tt>full_range(colRng,0,this->domain()->dim()-1)</tt> respectively.
   * </ul>
   *
   * <b>Note:</b> This view is to be used immediately and then committed back
   * with a call to <tt>commitDetachedView()</tt>.
   *
   * Note that calling this operation might require some internal allocations
   * and temporary memory.  Therefore, it is critical that
   * <tt>this->commitDetachedView(sub_mv)</tt> is called to commit the
   * changed entries and clean up memory and avoid memory leaks after the
   * sub-multi-vector is modified.
   *
   * <b>Heads Up!</b> Note that client code in general should not directly
   * call this function but should instead use the utility class
   * <tt>DetachedMultiVectorView</tt> which will also take care of
   * calling <tt>commitDetachedView</tt>.
   *
   * If <tt>this->acquireDetachedView(...,sub_mv)</tt> was previously
   * called on <tt>sub_mv</tt> then it may be possible to reuse this
   * memory if it is sufficiently sized.  The user is encouraged to
   * make multiple calls to
   * <tt>this->acquireDetachedView(...,sub_mv)</tt> before
   * <tt>this->commitDetachedView(sub_mv)</tt> to finally clean up
   * all of the memory.  Of course the same <tt>sub_mv</tt> object
   * must be passed to the same multi-vector object for this to work
   * correctly.
   *
   * Changes to the underlying sub-multi-vector are not guaranteed to
   * become permanent until <tt>this->acquireDetachedView(...,sub_mv)</tt>
   * is called again, or <tt>this->commitDetachedView(sub_mv)</tt> is
   * called.
   *
   * This function has a default implementation based on the vector
   * operation <tt>VectorBase::acquireDetachedView()</tt> called on the changeable
   * vector objects returned from <tt>col()</tt>.  Note that the
   * footprint of the reduction object (both internal and external
   * state) will be O(<tt>rowRng.size()*colRng.size()</tt>).  For
   * serial applications this is fairly reasonable and will not be a
   * major performance penalty.  For parallel applications, however,
   * this is a terrible implementation and must be overridden if
   * <tt>rowRng.size()</tt> is large at all.  Although, this function
   * should not even be used in case where the multi-vector is very
   * large.  If a subclass does override this function, it must also
   * override <tt>commitDetachedView()</tt> which has a default
   * implementation which is a companion to this function's default
   * implementation.
   */
  virtual void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar> *sub_mv
    ) = 0;

  /** \brief Commit changes for a changeable explicit view of a sub-multi-vector.
   *
   * \param sub_mv * [in/out] The data in <tt>sub_mv->values()</tt> will be
   * written back to internal storage and the memory referred to by
   * <tt>sub_mv->values()</tt> will be released if it was allocated * and
   * <tt>*sub_mv</tt> will be zeroed out using *
   * <tt>sub_mv->set_uninitialized()</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->range().get()!=NULL && this->domain().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_mv</tt> must have been passed through a call to 
   *      <tt>this->acquireDetachedView(...,sub_mv)</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> See <tt>RTOpPack::SubMultiVectorView::set_uninitialized()</tt> for <tt>sub_mv</tt>
   * <li> <tt>*this</tt> will be updated according the the changes made to <tt>sub_mv</tt>
   * </ul>
   *
   * The sub-multi-vector view must have been allocated by
   * <tt>this->acquireDetachedView()</tt> first.
   *
   * This function has a default implementation which is a companion
   * to the default implementation for <tt>acquireDetachedView()</tt>.  If
   * <tt>acquireDetachedView()</tt> is overridden by a subclass then this
   * function must be overridden also!
   */
  virtual void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    ) = 0;

  //@}

public:

#ifndef THYRA_HIDE_DEPRECATED_CODE
  /** @name Deprecated public functions */
  //@{

  /** \brief Deprecated. */
  THYRA_DEPRECATED RCP<const MultiVectorBase<Scalar> >
  subView( const int numCols, const int cols[] ) const
    { return subView( Teuchos::arrayView<const int>(cols,numCols) ); }

  /** \brief Deprecated. */
  THYRA_DEPRECATED RCP<MultiVectorBase<Scalar> >
  subView( const int numCols, const int cols[] )
    { return subView( Teuchos::arrayView<const int>(cols,numCols) ); }

  /** \brief Deprecated. */
  THYRA_DEPRECATED void applyOp(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const int num_multi_vecs,
    const MultiVectorBase<Scalar>*const multi_vecs_in[],
    const int num_targ_multi_vecs,
    MultiVectorBase<Scalar>*const targ_multi_vecs_inout[],
    RTOpPack::ReductTarget*const reduct_objs_inout[],
    const Ordinal primary_global_offset
    ) const;

  /** \brief Deprecated. */
  THYRA_DEPRECATED void applyOp(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const RTOpPack::RTOpT<Scalar> &secondary_op,
    const int num_multi_vecs,
    const MultiVectorBase<Scalar>*const multi_vecs_in[],
    const int num_targ_multi_vecs,
    MultiVectorBase<Scalar>*const targ_multi_vecs_inout[],
    RTOpPack::ReductTarget* reduct_obj,
    const Ordinal primary_global_offset
    ) const;

  //@}
#endif // THYRA_HIDE_DEPRECATED_CODE
private:
  
  // Not defined and not to be called
  MultiVectorBase<Scalar>&
  operator=(const MultiVectorBase<Scalar>&);

};


/** \brief Apply a reduction/transformation operator column by column and
 * return an array of the reduction objects.
 *
 * ToDo: Finish documentation!
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal primary_global_offset = 0
  )
{
  if(multi_vecs.size())
    multi_vecs[0]->applyOp(primary_op, multi_vecs, targ_multi_vecs,
      reduct_objs, primary_global_offset);
  else if(targ_multi_vecs.size())
    targ_multi_vecs[0]->applyOp(primary_op, multi_vecs, targ_multi_vecs,
      reduct_objs, primary_global_offset);
}


/** \brief Apply a reduction/transformation operator column by column and
 * reduce the intermediate reduction objects into one reduction object.
 *
 * ToDo: Finish documentation!
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const RTOpPack::RTOpT<Scalar> &secondary_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal primary_global_offset = 0
  )
{
  if(multi_vecs.size())
    multi_vecs[0]->applyOp(primary_op, secondary_op, multi_vecs, targ_multi_vecs,
      reduct_obj, primary_global_offset);
  else if(targ_multi_vecs.size())
    targ_multi_vecs[0]->applyOp(primary_op, secondary_op, multi_vecs, targ_multi_vecs,
      reduct_obj, primary_global_offset);
}

#ifndef THYRA_HIDE_DEPRECATED_CODE
//
// Deprecated non-members
//


/** \brief Deprecated.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
THYRA_DEPRECATED inline
void applyOp(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const int num_multi_vecs,
  const MultiVectorBase<Scalar>*const multi_vecs[],
  const int num_targ_multi_vecs,
  MultiVectorBase<Scalar>*const targ_multi_vecs[],
  RTOpPack::ReductTarget*const reduct_objs[],
  const Ordinal primary_global_offset = 0
  )
{
  if(num_multi_vecs)
    multi_vecs[0]->applyOp(
      primary_op,
      num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs,
      reduct_objs, primary_global_offset);
  else if(num_targ_multi_vecs)
    targ_multi_vecs[0]->applyOp(
      primary_op,
      num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs,
      reduct_objs, primary_global_offset);
}


/** \brief Deprecated.
 *
 * ToDo: Finish documentation!
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
THYRA_DEPRECATED inline
void applyOp(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const RTOpPack::RTOpT<Scalar> &secondary_op,
  const int num_multi_vecs,
  const MultiVectorBase<Scalar>*const multi_vecs[],
  const int num_targ_multi_vecs,
  MultiVectorBase<Scalar>*const targ_multi_vecs[],
  RTOpPack::ReductTarget *reduct_obj,
  const Ordinal primary_global_offset = 0
  )
{
  if(num_multi_vecs)
    multi_vecs[0]->applyOp(
      primary_op, secondary_op,
      num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs,
      reduct_obj, primary_global_offset);
  else if(num_targ_multi_vecs)
    targ_multi_vecs[0]->applyOp(
      primary_op, secondary_op,
      num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs,
      reduct_obj, primary_global_offset);
}

#endif // THYRA_HIDE_DEPRECATED_CODE

} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_BASE_DECL_HPP
