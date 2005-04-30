// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_DECL_HPP
#define THYRA_MULTI_VECTOR_DECL_HPP

#include "Thyra_LinearOpBaseDecl.hpp"
//#include "Thyra_LinearOpBase.hpp"
#include "RTOpPack_RTOpT.hpp"

namespace Thyra {

/** \brief Interface for a collection of column vectors called a multi-vector.
 *
 * The primary purpose for this interface is to allow for convienent
 * aggregations of column vectors.  Such an orderly arrangement of
 * column vectors into a single aggregate object allows for better
 * optimized linear algebra operations such as matrix-matrix
 * multiplication and the solution of linear systems for multiple
 * right hand sides.  Every computing environment (serial, parallel,
 * out-of-core etc.) should be able to define at least one reasonbly
 * efficient implementation of <tt>%MultiVectorBase</tt>.
 *
 * The <tt>%MultiVectorBase</tt> interface is derived from the
 * <tt>LinearOpBase</tt> interface and therefore a <tt>%MultiVectorBase</tt>
 * can be considered as a linear operator which has some interesting
 * implications.
 *
 * Note that since, this interface is derived from <tt>LinearOpBase</tt>
 * that it must support the functions <tt>domain()</tt> and
 * <tt>range()</tt>.
 *
 * Another very powerful feature of this interface is the ability to
 * apply reduction/transformation operators over a sub-set of rows and
 * columns in a set of multi-vector objects.  The behavior is
 * identical to the client extracted the rows or columns in a set of
 * multi-vectors and called <tt>VectorBase::applyOp()</tt> itself.
 * However, the advantage of using the multi-vector functions is that
 * there may be greater opportunities for increased performance in a
 * number of respects.  Also, the intermediate reduction objects over
 * a set of rows or columns can be reduced by a secondary reduction
 * object.
 *
 * This interface also allows a client to extract a sub-set of
 * elements in an explicit form as non-mutable
 * <tt>RTOpPack::SubMultiVectorT<Scalar></tt> or mutable
 * <tt>RTOpPack::MutableSubMultiVectorT<Scalar></tt> objects using the
 * operations <tt>getSubMultiVector()</tt>.  In general, this is a
 * very bad thing to do and should be avoided at all costs.  However,
 * there are some situations where this is needed and therefore it is
 * supported.  The default implementation of these operations use
 * reduction/transformation operators with <tt>applyOp()</tt> in order
 * to extract and set the needed elements and therefore all
 * <tt>%MultiVectorBase</tt> subclasses automatically support these
 * operations (even if it is a bad idea to use them).  <b>Heads
 * Up!</b> Note that client code in general should not directly call
 * the above explicit sub-vector access functions but should use the
 * utility classes <tt>ExplicitMultiVectorView</tt> and
 * <tt>ExplicitMutableMultiVectorView</tt> instead since these are
 * easier an safer in the event that an exception is thrown.
 *
 * <b>Notes for subclass developers</b>
 *
 * Only one function override is required (in addition to the functions
 * that must be overridden from <tt>LinearOpBase</tt>, except
 * <tt>apply()</tt>, see below) to create a concrete
 * <tt>MultiVectorBase</tt> subclass and that is the non-constant version
 * of <tt>col()</tt>.  All of the other functions have default
 * implementations that are quite good in many use cases.
 *
 * Note that through the magic of the <tt>applyOp()</tt> functions that
 * this interface is able to implement the pure virtual
 * <tt>apply()</tt> function from the <tt>LinearOpBase</tt> interface.  This
 * implementation is not optimal, but will be sufficient in many
 * different contexts.
 *
 * The <tt>applyOp()</tt> functions should only be overridden if the
 * subclass can do something more efficient than simply applying the
 * reduction/transformation operators one column at a time.
 *
 * The non-const versions of <tt>subView()</tt> should be overridden
 * if the performance of these functions are important.  The default
 * implementation will not achieve near-optimal performance.
 *
 * Functions that subclasses should almost never need (or want) to
 * override are the const version of <tt>col()</tt> or the the const
 * versions of <tt>subView()</tt>.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class MultiVectorBase : virtual public LinearOpBase<Scalar> {
public:

  /** \brief . */
  using LinearOpBase<Scalar>::describe;
  /** \brief . */
  using LinearOpBase<Scalar>::apply;
  
  /** @name Provide access to the columns as VectorBase objects */
  //@{

  /** \brief Get a non-mutable view of a constituent column vector.
   *
   * \param  j  [in] One-based index the view column to return.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>1 <= j && j <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>this->range()->isCompatible(*return->space()) == true</tt>
   * </ul>
   *
   * The default implementation of this function (which is the only
   * implementaion needed by most subclasses) is based on the
   * non-const version <tt>col()</tt>.
   */
  virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > col(Index j) const;

  /** \brief Get a mutable view of a column vector.
   *
   * \param  j  [in] One-based Index the view column to return.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>1 <= j && j <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>this->range()->isCompatible(*return->space()) == true</tt>
   * </ul>
   *
   * Note that <tt>*this</tt> is not guaranteed to be modified until
   * the smart pointer returned by this function is destroyed.
   */
  virtual Teuchos::RefCountPtr<VectorBase<Scalar> > col(Index j) = 0;

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
  virtual Teuchos::RefCountPtr<MultiVectorBase<Scalar> > clone_mv() const;

  //@}

  /** @name Sub-view functions */
  //@{

  /** \brief Return a const sub-view of a contiguous set of columns of the this multi-vector.
   *
   * @param  colRng  [in] Range of columns to create a view of.  Note that it is valid for
   *                  <tt>colRng.full_range()==true</tt> in which case the view of the entire
   *                  multi-vector is taken.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>colRng.ubound() <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == RangePack::full_range(colRng,1,this->domain()->dim()).size()</tt>
   * <li> <tt>*return->col(1+k)</tt> represents the same column vector as <tt>this->col(colRng.lbound()+k)</tt>,
   *      for <tt>k=0...RangePack::full_range(colRng,1,this->domain()->dim()).ubound()</tt>
   * </ul>
   *
   * The default implementation (which is the only implementation
   * needed by most subclasses) is to return object from the non-const
   * verstion <tt>subView()</tt>.
   */
  virtual Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const Range1D& colRng ) const;
  
  /** \brief Return a non-const sub-view of a contiguous set of columns of the this multi-vector.
   *
   * @param  colRng  [in] Range of columns to create a view of.  Note that it is valid for
   *                  <tt>colRng.full_range()==true</tt> in which case the view of the entire
   *                  multi-vector is taken.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>colRng.ubound() <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == RangePack::full_range(colRng,1,this->domain()->dim()).size()</tt>
   * <li> <tt>*return->col(1+k)</tt> represents the same column vector as <tt>this->col(colRng.lbound()+k)</tt>,
   *      for <tt>k=0...RangePack::full_range(colRng,1,this->domain()->dim()).ubound()</tt>
   * </ul>
   *
   * The default implementation of this function uses <tt>MultiVectorCols</tt> but this is not
   * a good default implementation in general.
   *
   * Note that <tt>*this</tt> is not guaranteed to be modified until
   * the smart pointer returned by this function goes out of scope.
   */
  virtual Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const Range1D& colRng );

  /** \brief Return a const sub-view of a non-contiguous set of columns of this multi-vector.
   *
   * @param  numCols  [in] The number of columns to extract a view for.
   * @param  cols     [in] Array (length <tt>numCols</tt>) of the columns indices to use in the
   *                  returned view.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>numCols <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>1 <= cols[k] <= this->domain()->dim()</tt>, for <tt>k=0...numCols-1</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>col[k1] != col[k2]</tt>, for all <tt>k1 != k2</tt> in the range <tt>[0,numCols-1]</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == numCols</tt>
   * <li> <tt>*return->col(k+1)</tt> represents the same column vector as <tt>this->col(cols[k])</tt>,
   *      for <tt>k=0...numCols</tt>
   * </ul>
   *
   * The default implementation (which is the only implementation
   * needed by most subclasses) is to return object from the non-const
   * verstion <tt>subView()</tt>.
   */
  virtual Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] ) const;

  /** \brief Return a non-const sub-view of a non-contiguous set of columns of this multi-vector.
   *
   * @param  numCols  [in] The number of columns to extract a view for.
   * @param  cols     [in] Array (length <tt>numCols</tt>) of the columns indices to use in the
   *                  returned view.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>1 <= numCols <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>1 <= cols[k] <= this->domain()->dim()</tt>, for <tt>k=0...numCols-1</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>col[k1] != col[k2]</tt>, for all <tt>k1 != k2</tt> in the range <tt>[0,numCols-1]</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
   * <li> <tt>return->domain()->dim() == numCols</tt>
   * <li> <tt>*return->col(k+1)</tt> represents the same column vector as <tt>this->col(cols[k])</tt>,
   *      for <tt>k=0...numCols</tt>
   * </ul>
   *
   * The default implementation of this function uses
   * <tt>MultiVectorCols</tt> but this is not a good default
   * implementation in general.
   *
   * Note that <tt>*this</tt> is not guaranteed to be modified until
   * the smart pointer returned by this function goes out of scope.
   */
  virtual Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] );
  
  //@}

  /** @name Collective applyOp() functions */
  //@{

  /** \brief Apply a reduction/transformation operator column by column and
   * return an array of the reduction objects.
   *
   * Preconditions:<ul>
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
   * each colunn <tt>this->col(j)</tt> for <tt>j = 1
   * ... this->range()->dim()</tt>.
   */
  virtual void applyOp(
    const RTOpPack::RTOpT<Scalar>   &primary_op
    ,const int                      num_multi_vecs
    ,const MultiVectorBase<Scalar>* multi_vecs[]
    ,const int                      num_targ_multi_vecs
    ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
    ,RTOpPack::ReductTarget*        reduct_objs[]
    ,const Index                    primary_first_ele
    ,const Index                    primary_sub_dim
    ,const Index                    primary_global_offset
    ,const Index                    secondary_first_ele
    ,const Index                    secondary_sub_dim
    ) const;

  /** \breif Apply a reduction/transformation operator column by column and reduce the intermediate
   * reduction objects into one reduction object.
   *
   * Preconditions:<ul>
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
  virtual void applyOp(
    const RTOpPack::RTOpT<Scalar>   &primary_op
    ,const RTOpPack::RTOpT<Scalar>  &secondary_op
    ,const int                      num_multi_vecs
    ,const MultiVectorBase<Scalar>* multi_vecs[]
    ,const int                      num_targ_multi_vecs
    ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
    ,RTOpPack::ReductTarget         *reduct_obj
    ,const Index                    primary_first_ele
    ,const Index                    primary_sub_dim
    ,const Index                    primary_global_offset
    ,const Index                    secondary_first_ele
    ,const Index                    secondary_sub_dim
    ) const;

  //@}

  /** @name Explicit sub-multi-vector access */
  //@{

  /** \brief Get a non-mutable explicit view of a sub-multi-vector.
   *
   * @param  rowRng   [in] The range of the rows to extract the sub-multi-vector view.
   * @param  colRng   [in] The range of the columns to extract the sub-multi-vector view.
   * @param  sub_mv   [in/out] View of the sub-multi_vector.  Prior to the
   *                  first call <tt>sub_mv->set_uninitialized()</tt> must be called.
   *                  Technically <tt>*sub_mv</tt> owns the memory but this memory can be freed
   *                  only by calling <tt>this->freeSubMultiVector(sub_mv)</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!rowRng.full_range()</tt>] <tt>(rowRng.ubound() <= this->range()->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>(colRng.ubound() <= this->domain()->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*sub_mv</tt> contains an explicit non-mutable view to the elements
   *      in the row and column ranges <tt>full_range(rowRng,1,this->range()->dim())</tt>
   *      and <tt>full_range(colRng,1,this->domain()->dim())</tt> respectively.
   * </ul>
   *
   * This is only a transient view of a sub-multi-vector that is to
   * be immediately used and then released with a call to
   * <tt>freeSubMultiVector()</tt>.
   *
   * Note that calling this operation might require some dynamic
   * memory allocations and temporary memory.  Therefore, it is
   * critical that <tt>this->freeSubMultiVector(sub_mv)</tt> is called to
   * clean up memory and avoid memory leaks after the
   * sub-multi-vector is used.
   *
   * <b>Heads Up!</b> Note that client code in general should not
   * directly call this function but should use the utility class
   * <tt>ExplicitMultiVectorView</tt> which will also take care of
   * calling <tt>freeSubMultiVector()</tt>.
   *
   * If <tt>this->getSubMultiVector(...,sub_mv)</tt> was previously
   * called on <tt>sub_mv</tt> then it may be possible to reuse this
   * memory if it is sufficiently sized.  The user is encouraged to
   * make multiple calls to
   * <tt>this->getSubMultiVector(...,sub_mv)</tt> before
   * <tt>this->freeSubMultiVector(sub_mv)</tt> to finally clean up all of
   * the memory.  Of course, the same <tt>sub_mv</tt> object must be
   * passed to the same multi-vector object for this to work correctly.
   *
   * This function has a default implementation based on the vector
   * operation <tt>VectorBase::getSubVector()</tt> called on the
   * non-mutable vector objects returned from <tt>col()</tt>.  Note
   * that the footprint of the reduction object (both internal and
   * external state) will be
   * O(<tt>rowRng.size()*colRng.size()</tt>).  For serial
   * applications this is faily reasonable and will not be a major
   * performance penalty.  For parallel applications, however, this
   * is a terrible implementation and must be overridden if
   * <tt>rowRng.size()</tt> is large at all.  Although, this function
   * should not even be used in case where the multi-vector is very
   * large.  If a subclass does override this function, it must also
   * override <tt>freeSubMultiVector()</tt> which has a default
   * implementation which is a companion to this function's default
   * implementation.
   */
  virtual void getSubMultiVector(
    const Range1D                       &rowRng
    ,const Range1D                      &colRng
    ,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
    ) const;

  /** \brief Free an explicit view of a sub-multi-vector.
   *
   * @param  sub_mv
   *				[in/out] The memory refered to by <tt>sub_mv->values()</tt>
   *				will be released if it was allocated and <tt>*sub_mv</tt>
   *              will be zeroed out using <tt>sub_mv->set_uninitialized()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_mv</tt> must have been passed through a call to 
   *      <tt>this->getSubMultiVector(...,sub_mv)</tt>
   * </ul>
    *
   * Postconditions:<ul>
   * <li> See <tt>RTOpPack::SubMultiVectorT::set_uninitialized()</tt> for <tt>sub_mv</tt>
   * </ul>
   *
   * The sub-multi-vector view must have been allocated by
   * <tt>this->getSubMultiVector()</tt> first.
   *
   * This function has a default implementation which is a companion
   * to the default implementation for <tt>getSubMultiVector()</tt>.  If
   * <tt>getSubMultiVector()</tt> is overridden by a subclass then this
   * function must be overridden also!
   */
  virtual void freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* sub_mv ) const;

  /** \brief Get a mutable explicit view of a sub-multi-vector.
   *
   * @param  rowRng   [in] The range of the rows to extract the sub-multi-vector view.
   * @param  colRng   [in] The range of the columns to extract the sub-multi-vector view.
   * @param  sub_mv   [in/out] Mutable view of the sub-multi-vector.  Prior to the
   *                  first call <tt>sub_mv->set_uninitialized()</tt> must
   *                  have been called for the correct behavior.  Technically
   *                  <tt>*sub_mv</tt> owns the memory but this memory must be commited
   *                  and freed only by calling <tt>this->commitSubMultiVector(sub_mv)</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>!rowRng.full_range()</tt>] <tt>(rowRng.ubound() <= this->range()->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * <li> [<tt>!colRng.full_range()</tt>] <tt>(colRng.ubound() <= this->domain()->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
    *
   * Postconditions:<ul>
   * <li> <tt>*sub_mv</tt> contains an explicit mutable view to the elements
   *      in the row and column ranges <tt>full_range(rowRng,1,this->range()->dim())</tt>
   *      and <tt>full_range(colRng,1,this->domain()->dim())</tt> respectively.
   * </ul>
   *
   * This is only a transient view of a sub-multi-vector that is to be
   * immediately used and then committed back with a call to
   * <tt>commitSubMultiVector()</tt>.
   *
   * Note that calling this operation might require some internal
   * allocations and temporary memory.  Therefore, it is critical
   * that <tt>this->commitSubMultiVector(sub_mv)</tt> is called to
   * commit the changed entires and clean up memory and avoid memory
   * leaks after the sub-multi-vector is modified.
   *
   * <b>Heads Up!</b> Note that client code in general should not
   * directly call this function but should use the utility class
   * <tt>ExplicitMutableMultiVectorView</tt> which will also take care
   * of calling <tt>commitSubMultiVector</tt>.
   *
   * If <tt>this->getSubMultiVector(...,sub_mv)</tt> was previously
   * called on <tt>sub_mv</tt> then it may be possible to reuse this
   * memory if it is sufficiently sized.  The user is encouraged to
   * make multiple calls to
   * <tt>this->getSubMultiVector(...,sub_mv)</tt> before
   * <tt>this->commitSubMultiVector(sub_mv)</tt> to finally clean up
   * all of the memory.  Of course the same <tt>sub_mv</tt> object
   * must be passed to the same multi-vector object for this to work
   * correctly.
   *
   * Changes to the underlying sub-multi-vector are not guarrenteed to
   * become permanent until <tt>this->getSubMultiVector(...,sub_mv)</tt>
   * is called agian, or <tt>this->commitSubMultiVector(sub_mv)</tt> is
   * called.
   *
   * This function has a default implementation based on the vector
   * operation <tt>VectorBase::getSubVector()</tt> called on the mutable
   * vector objects returned from <tt>col()</tt>.  Note that the
   * footprint of the reduction object (both internal and external
   * state) will be O(<tt>rowRng.size()*colRng.size()</tt>).  For
   * serial applications this is faily reasonable and will not be a
   * major performance penalty.  For parallel applications, however,
   * this is a terrible implementation and must be overridden if
   * <tt>rowRng.size()</tt> is large at all.  Although, this function
   * should not even be used in case where the multi-vector is very
   * large.  If a subclass does override this function, it must also
   * override <tt>commitSubMultiVector()</tt> which has a default
   * implementation which is a companion to this function's default
   * implementation.
   */
  virtual void getSubMultiVector(
    const Range1D                                &rowRng
    ,const Range1D                               &colRng
    ,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
    );

  /** \brief Commit changes for a mutable explicit view of a sub-multi-vector.
   *
   * @param sub_mv
   *				[in/out] The data in <tt>sub_mv->values()</tt> will be written
   *              back internal storage and the memory refered to by
   *              <tt>sub_mv->values()</tt> will be released if it was allocated
   *				and <tt>*sub_mv</tt> will be zeroed out using
   *				<tt>sub_mv->set_uninitialized()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>sub_mv</tt> must have been passed through a call to 
   *      <tt>this->getSubMultiVector(...,sub_mv)</tt>
   * </ul>
    *
   * Postconditions:<ul>
   * <li> See <tt>RTOpPack::MutableSubMultiVectorT::set_uninitialized()</tt> for <tt>sub_mv</tt>
   * <li> <tt>*this</tt> will be updated according the the changes made to <tt>sub_mv</tt>
   * </ul>
   *
   * The sub-multi-vector view must have been allocated by
   * <tt>this->getSubMultiVector()</tt> first.
   *
   * This function has a default implementation which is a companion
   * to the default implementation for <tt>getSubMultiVector()</tt>.  If
   * <tt>getSubMultiVector()</tt> is overridden by a subclass then this
   * function must be overridden also!
   */
  virtual void commitSubMultiVector( RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv );

  //@}

  /** @name Overridden functions from LinearOpBase */
  //@{

  /** \breif This function is implemented in terms of the mulit-vector <tt>applyOp()</tt> function.
   *
   * The implementation takes care of two types of operations.  One
   * (<tt>M_trans==TRANS</tt>) is the block dot product of two
   * vectors to form scalar (stored as the vector <tt>y</tt> which
   * as one component).  The other (<tt>M_trans==NOTRANS</tt>) is
   * essentially an axpy operation where <tt>x</tt> is a vector with
   * one element.  Both of these operations are performed using
   * reduction/transformation operators.
   *
   * This implementation is near optimal but the default
   * implementation of the multi-vector verstion of <tt>apply()</tt>
   * as implemented in the base class <tt>LinearOpBase</tt> will not be
   * a near optimal implementation in its current form do to
   * multiple, sequential reductions but it could be made to be with
   * a little work.
   */
  void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const;

  /// This function is simply overridden to return <tt>this->clone_mv()</tt>.
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;

  //@}

private:
  
#ifdef DOXYGEN_COMPILE
  VectorBase<Scalar>  *columns; // doxygen only
#endif	

}; // end class MultiVectorBase<Scalar>

/** \brief Apply a reduction/transformation operator column by column and
 * return an array of the reduction objects.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar>   &primary_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget*        reduct_objs[]
#ifdef __sun
  ,const Index                    primary_first_ele
  ,const Index                    primary_sub_dim
  ,const Index                    primary_global_offset
  ,const Index                    secondary_first_ele
  ,const Index                    secondary_sub_dim
#else
  ,const Index                    primary_first_ele      = 1
  ,const Index                    primary_sub_dim        = 0
  ,const Index                    primary_global_offset  = 0
  ,const Index                    secondary_first_ele    = 1
  ,const Index                    secondary_sub_dim      = 0
#endif
  )
{
  if(num_multi_vecs)
    multi_vecs[0]->applyOp(
      primary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_objs,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
  else if(num_targ_multi_vecs)
    targ_multi_vecs[0]->applyOp(
      primary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_objs,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
}

#ifdef __sun
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar>   &primary_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget*        reduct_objs[]
  )
{
  applyOp(
          primary_op
          ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
          ,reduct_objs,1,0,0,1,0
          );
}
#endif

/** \brief Apply a reduction/transformation operator column by column and reduce the intermediate
 * reduction objects into one reduction object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar>   &primary_op
  ,const RTOpPack::RTOpT<Scalar>  &secondary_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
#ifdef __sun
  ,const Index                    primary_first_ele
  ,const Index                    primary_sub_dim
  ,const Index                    primary_global_offset
  ,const Index                    secondary_first_ele
  ,const Index                    secondary_sub_dim
#else
  ,const Index                    primary_first_ele      = 1
  ,const Index                    primary_sub_dim        = 0
  ,const Index                    primary_global_offset  = 0
  ,const Index                    secondary_first_ele    = 1
  ,const Index                    secondary_sub_dim      = 0
#endif
  )
{
  if(num_multi_vecs)
    multi_vecs[0]->applyOp(
      primary_op,secondary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_obj,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
  else if(num_targ_multi_vecs)
    targ_multi_vecs[0]->applyOp(
      primary_op,secondary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_obj,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
}

#ifdef __sun

template<class Scalar>
inline
void applyOp(
  const RTOpPack::RTOpT<Scalar>   &primary_op
  ,const RTOpPack::RTOpT<Scalar>  &secondary_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  )
{
  applyOp(
      primary_op,secondary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_obj,1,0,0,1,0
      );
}

#endif


} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_DECL_HPP
