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

#ifndef THYRA_VECTOR_MULTI_VECTOR_DECL_HPP
#define THYRA_VECTOR_MULTI_VECTOR_DECL_HPP

#include "Thyra_VectorDefaultBaseDecl.hpp"

namespace Thyra {

/** \brief Generic adapter subclass that takes any <tt>MultiVectorBase</tt>
 * and turns it into a <tt>VectorBase</tt> object where columns vectors are
 * stacked on top of one another to make one big vector.
 *
 * There are two primary purposes for this concrete subclass.  The first
 * purpose of this subclass is to provide an implementation for
 * <tt>VectorBase</tt> given that a concrete implementation for a
 * <tt>MultiVectorBase</tt> is already provided.  A linear algebra library
 * implementation should have to do almost nothing to get a
 * <tt>VectorBase</tt> implementation if a <tt>MultiVectorBase</tt> is already
 * supported.  The second purpose of this subclass is to take any
 * <tt>MultiVectorBase</tt> object with multiple columns and make it look like
 * one big vector.
 *
 * To use this subclass to provide an implementation for <tt>VectorBase</tt>,
 * implement the override of the <tt>VectorSpaceBase::createMember()</tt>
 * function as:
 *
 \code

  template<class Scalar>
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > SomeVectorSpace::createMember()
  {
    return Teuchos::rcp(new DefaultVectorMultiVector<Scalar>(this->createMembers(1)));
  }
 \endcode
 *
 * where <tt>SomeVectorSpace::createMembers(int)</tt> is overridden to
 * create the multi-vector object.
 *
 * ToDo: the functionality to support the second and more general use
 * case is not finished yet but can be put together when needed.
 */
template<class Scalar>
class DefaultVectorMultiVector : virtual public VectorDefaultBase<Scalar> {
public:

  /** \brief . */
  //using VectorBase<Scalar>::col;

  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct uninitialized (see the post-conditions for <tt>uninitialize()</tt>).
  DefaultVectorMultiVector();

  /// Calls <tt>initialize()</tt>.
  DefaultVectorMultiVector(
    const Teuchos::RefCountPtr<MultiVectorBase<Scalar> > &mv
    );

  /** \brief Initialize given a MultiVectorBase object.
   *
   * Preconditions:<ul>
   * <li><tt>mv.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li><tt>mv->domain().get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->mv().get() == mv.get()</tt>
   * <tt>[<tt>mv->domain()->dim()==1</tt>] <tt>this->space().get() == mv->range().get()</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<MultiVectorBase<Scalar> > &mv
    );

  /** \brief Set to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->mv().get() == NULL</tt>
   * <tt><tt>this->space().get() == NULL</tt>
   * </ul>
   */
  void uninitialize(
    Teuchos::RefCountPtr<MultiVectorBase<Scalar> > *mv = NULL
    );

  /// Return smart pointer to non-const reference to underlying <tt>MultiVectorBase</tt> object.
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > mv();

  /// Return smart pointer to const reference to underlying <tt>MultiVectorBase</tt> object.
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > mv() const;

  //@}

  /** @name Overridden from LinearOpBase (forwarded to this->mv()) */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  /** \brief . */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  //@}

  /** @name Overridden from MultiVectorBase (forwarded to this->mv()) */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr<VectorBase<Scalar> > col(Index j);
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > clone_mv() const;
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const Range1D& col_rng ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const Range1D& col_rng );
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] );
  /** \brief . */
  void applyOp(
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
  /** \brief . */
  void applyOp(
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
  /** \brief . */
  void acquireDetachedView(
    const Range1D                       &rowRng
    ,const Range1D                      &colRng
    ,RTOpPack::ConstSubMultiVectorView<Scalar>  *sub_mv
    ) const;
  /** \brief . */
  void releaseDetachedView( RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv ) const;
  /** \brief . */
  void acquireDetachedView(
    const Range1D                                &rowRng
    ,const Range1D                               &colRng
    ,RTOpPack::SubMultiVectorView<Scalar>    *sub_mv
    );
  /** \brief . */
  void commitDetachedView( RTOpPack::SubMultiVectorView<Scalar>* sub_mv );
  //@}

  /** @name Overridden from VectorBase (defined in terms of this->mv()) */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > space() const;
  /** \brief . */
  void applyOp(
    const RTOpPack::RTOpT<Scalar>   &op
    ,const int                      num_vecs
    ,const VectorBase<Scalar>*      vecs[]
    ,const int                      num_targ_vecs
    ,VectorBase<Scalar>*            targ_vecs[]
    ,RTOpPack::ReductTarget         *reduct_obj
    ,const Index                    first_ele
    ,const Index                    sub_dim
    ,const Index                    global_offset
    ) const;
  /** \brief . */
  void acquireDetachedView( const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const;
  /** \brief . */
  void releaseDetachedView( RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const;
  /** \brief . */
  void acquireDetachedView( const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec );
  /** \brief . */
  void commitDetachedView( RTOpPack::SubVectorView<Scalar>* sub_vec );
  //@}

private:
  
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > mv_;

}; // end class DefaultVectorMultiVector

// ///////////////////////////////////////////////////
// Inline members

template <class Scalar>
inline
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::mv()
{
  return mv_;
}

template <class Scalar>
inline
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::mv() const
{
  return mv_;
}

} // end namespace Thyra

#endif // THYRA_VECTOR_MULTI_VECTOR_DECL_HPP
