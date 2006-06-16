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

#ifndef THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_MultiVectorBaseDecl.hpp"
#include "Thyra_LinearOpDefaultBaseDecl.hpp"

namespace Thyra {

/** \brief Node subclass that uses a default <tt>MultiVectorBase</tt>
 * implementation</tt> to provide default implementations for as many other
 * functions in <tt>MultiVectorBase</tt> interface the as is reasonable.
 *
 * <b>Notes to subclass developers</b>
 *
 * Only three function overrides are required in order to create a concrete
 * <tt>MultiVectorBase</tt> subclass: <tt>range()</tt>, <tt>domain()</tt> and
 * the non-const version of <tt>col()</tt>.  All of the other functions have
 * default implementations.  However, a good implementation will provide
 * optimized overrides of at least the functions <tt>apply()</tt> and
 * <tt>applyTranspose()</tt>.  The non-const versions of <tt>subView()</tt>
 * should be overridden if subviews are important.  The default implementation
 * will not achieve near-optimal performance in many cases.
 *
 * \ingroup Thyra_Op_Vec_Adapters_grp
 */
template<class Scalar>
class MultiVectorDefaultBase
  : virtual public MultiVectorBase<Scalar>
  , virtual protected LinearOpDefaultBase<Scalar>
{
public:

  /** \brief . */
  using LinearOpDefaultBase<Scalar>::describe;
  /** \brief . */
  using MultiVectorBase<Scalar>::applyOp;
  /** \brief . */
  using MultiVectorBase<Scalar>::col;
  /** \brief . */
  using MultiVectorBase<Scalar>::subView;

  /** \name Overridden public member functions from MultiVectorBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const Range1D& colRng ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const Range1D& colRng );
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] );
  /** \brief .
   *
   * This implementation calls <tt>VectorBase::applyOp()</tt> on each column
   * <tt>this->col(j)</tt> for <tt>j = 0 ... this->range()->dim()-1</tt>.
   */
  virtual void applyOp(
    const RTOpPack::RTOpT<Scalar>         &primary_op
    ,const int                            num_multi_vecs
    ,const MultiVectorBase<Scalar>*const  multi_vecs[]
    ,const int                            num_targ_multi_vecs
    ,MultiVectorBase<Scalar>*const        targ_multi_vecs[]
    ,RTOpPack::ReductTarget*const         reduct_objs[]
    ,const Index                          primary_first_ele_offset
    ,const Index                          primary_sub_dim
    ,const Index                          primary_global_offset
    ,const Index                          secondary_first_ele_offset
    ,const Index                          secondary_sub_dim
    ) const;
  /** \brief .
   *
   * This implementation calls <tt>applyOp()</tt> where an array of reduction
   * objects is taken.
   */
  virtual void applyOp(
    const RTOpPack::RTOpT<Scalar>         &primary_op
    ,const RTOpPack::RTOpT<Scalar>        &secondary_op
    ,const int                            num_multi_vecs
    ,const MultiVectorBase<Scalar>*const  multi_vecs[]
    ,const int                            num_targ_multi_vecs
    ,MultiVectorBase<Scalar>*const        targ_multi_vecs[]
    ,RTOpPack::ReductTarget               *reduct_obj
    ,const Index                          primary_first_ele_offset
    ,const Index                          primary_sub_dim
    ,const Index                          primary_global_offset
    ,const Index                          secondary_first_ele_offset
    ,const Index                          secondary_sub_dim
    ) const;
  /** \brief .
   *
   * This implementation is based on the vector operation
   * <tt>VectorBase::acquireDetachedView()</tt> called on the non-changeable vector
   * objects returned from <tt>col()</tt>.  Note that the footprint of the
   * reduction object (both internal and external state) will be
   * O(<tt>rowRng.size()*colRng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, however, this is a terrible implementation and
   * must be overridden if <tt>rowRng.size()</tt> is large at all.  Although,
   * this function should not even be used in cases where the multi-vector is
   * very large.  If a subclass does override this function, it must also
   * override <tt>releaseDetachedView()</tt> which has an implementation
   * which is a companion to this function's implementation.
   */
  virtual void acquireDetachedView(
    const Range1D                       &rowRng
    ,const Range1D                      &colRng
    ,RTOpPack::ConstSubMultiVectorView<Scalar>  *sub_mv
    ) const;
  /** \brief .
   *
   * This implementation is a companion to the implementation for
   * <tt>acquireDetachedView()</tt>.  If <tt>acquireDetachedView()</tt> is
   * overridden by a subclass then this function must be overridden also!
   */
  virtual void releaseDetachedView( RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv ) const;
  /** \brief .
   *
   * This implementation is based on the vector operation
   * <tt>VectorBase::acquireDetachedView()</tt> called on the changeable vector
   * objects returned from <tt>col()</tt>.  Note that the footprint of the
   * reduction object (both internal and external state) will be
   * O(<tt>rowRng.size()*colRng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, however, this is a terrible implementation and
   * must be overridden if <tt>rowRng.size()</tt> is large at all.  Although,
   * this function should not even be used in case where the multi-vector is
   * very large.  If a subclass does override this function, it must also
   * override <tt>commitDetachedView()</tt> which has an implementation
   * which is a companion to this function's implementation.
   */
  virtual void acquireDetachedView(
    const Range1D                                &rowRng
    ,const Range1D                               &colRng
    ,RTOpPack::SubMultiVectorView<Scalar>    *sub_mv
    );
  /** \brief .
   *
   * This implementation is a companion to the default implementation for
   * <tt>acquireDetachedView()</tt>.  If <tt>acquireDetachedView()</tt> is
   * overridden by a subclass then this function must be overridden also!
   */
  virtual void commitDetachedView( RTOpPack::SubMultiVectorView<Scalar>* sub_mv );
  /** \brief .
   *
   * This implementation uses the vector space to create a new multi-vector
   * object and then uses a transformation operator to assign the vector
   * elements.  A subclass should only override this function if it can do
   * something more sophisticated (i.e. lazy evaluation) but in general, this
   * is not needed.
   */
  virtual Teuchos::RefCountPtr<MultiVectorBase<Scalar> > clone_mv() const;
  //@}

};

} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
