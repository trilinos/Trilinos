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

#ifndef THYRA_DEFAULT_COLUMNWISE_MULTI_VECTOR_DECL_HPP
#define THYRA_DEFAULT_COLUMNWISE_MULTI_VECTOR_DECL_HPP

#include "Thyra_MultiVectorDefaultBase_decl.hpp"
#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_VectorBase.hpp"


namespace Thyra {


/** \brief Default subclass for <tt>MultiVectorBase</tt> implemented using columns
 * of separate abstract vectors.
 *
 * This is a very bad implementation of a multi-vector but this will
 * work in situations where you need a multi-vector but some
 * underlying linear algebra library does not directly support them.
 *
 * This subclass can be used to represent a <tt>%MultiVectorBase</tt>
 * wrapper around a single <tt>VectorBase</tt> object so that a single
 * vector can be passed to a method that expects a <tt>%MultiVectorBase</tt>
 * object.
 */
template<class Scalar>
class DefaultColumnwiseMultiVector : virtual public MultiVectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/Initializers */
  //@{

  /** \brief Construct to initialized.
   *
   * Postconditions:<ul>
   * <tt> <tt>this->range().get() == NULL</tt>
   * <tt> <tt>this->domain().get() == NULL</tt>
   * </ul>
   */
  DefaultColumnwiseMultiVector();

  /** \brief Calls <tt>initialize()</tt>. */
  DefaultColumnwiseMultiVector(
    const RCP<VectorBase<Scalar> > &col_vec
    );

  /** \brief Calls <tt>initialize()</tt>. */
  DefaultColumnwiseMultiVector(
    const RCP<const VectorSpaceBase<Scalar> > &range,
    const RCP<const VectorSpaceBase<Scalar> > &domain,
    const ArrayView<const RCP<VectorBase<Scalar> > > &col_vecs = Teuchos::null
    );
  
  /** \brief Initialize given a single vector object.
   *
   * \param col_vec [in] A single column vector.  It is not allowed for
   * <tt>col_vecs==NULL</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>col_vec.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>col_vec->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <tt> <tt>this->range().get() == col_vec.space().get()</tt>
   * <tt> <tt>this->domain()->dim() == 1</tt>
   * <li> <tt>this->col(0).get() == col_vec.get()</tt>
   * </ul>
   */
  void initialize(
    const RCP<VectorBase<Scalar> > &col_vec
    );

  /** \brief Initialize given the spaces for the columns and rows and possibly
   * the column vectors.
   *
   * \param range [in] The space that the columns must lie in.  The underlying
   * vector space must not be changed while <tt>this</tt> is in use.
   *
   * \param domain [in] The space that the rows must lie in.  The underlying
   * vector space must not be changed while <tt>this</tt> is in use.  What
   * this argument really specifies is what vector type will be compatible
   * with the vectors that the client may try to use to interact with the rows
   * of this multivector.
   *
   * \param col_vecs [in] Array (size <tt>domain->dim()</tt>) of column
   * vectors to use for the columns of <tt>this</tt>.  It is allowed for
   * <tt>col_vecs==NULL</tt> in which case <tt>range->createMember()</tt> will
   * be used to create the columns of <tt>this</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>range.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>domain.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>range->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>domain->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>col_vecs != NULL</tt>]
   *      <tt>col_vecs[j].get() != NULL &&
   *        col_vecs[j]->space()->is_compatible(*range) == true</tt>,
   *      for <tt>j=0..domain->dim()-1</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <tt> <tt>this->range().get() == range.get()</tt>
   * <tt> <tt>this->domain().get() == domain.get()</tt>
   * <li> [<tt>col_vecs != NULL</tt>] 
   *        <tt>this->col(j).get() == col_vecs[j].get()</tt>,
   *        for <tt>j=0..domain->dim()-1</tt>
   * </ul>
   */
  void initialize(
    const RCP<const VectorSpaceBase<Scalar> > &range,
    const RCP<const VectorSpaceBase<Scalar> > &domain,
    const ArrayView<const RCP<VectorBase<Scalar> > > &col_vecs = Teuchos::null
    );

  /** \brief Set uninitialized. */
  void uninitialize();

  //@}

  /** @name Overridden from LinearOpBAse */
  //@{
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > domain() const;
  //@}

  /** @name Overridden from MultiVectorBase */
  //@{
  /** \brief . */
  RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j);
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl( const Range1D& col_rng );
  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief For complex <tt>Scalar</tt> types returns <tt>true</tt> for
   * <tt>NOTRANS</tt> and <tt>CONJTRANS</tt> and for real types returns true
   * for all values of <tt>M_trans</tt>.
   */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief This function is implemented in terms of the multi-vector
   * <tt>applyOp()</tt> function.
   *
   * The implementation takes care of two types of operations.  One
   * (<tt>M_trans==TRANS</tt>) is the block dot product of two vectors to form
   * scalar (stored as the vector <tt>y</tt> which as one component).  The
   * other (<tt>M_trans==NOTRANS</tt>) is essentially an axpy operation where
   * <tt>x</tt> is a vector with one element.  Both of these operations are
   * performed using reduction/transformation operators.
   *
   * This implementation is near optimal but the default implementation of the
   * multi-vector version of <tt>apply()</tt> as implemented in the base class
   * <tt>LinearOpBase</tt> will not be a near optimal implementation in its
   * current form do to multiple, sequential reductions but it could be made
   * to be with a little work.
   */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:
  
  RCP<const VectorSpaceBase<Scalar> > range_;
  RCP<const VectorSpaceBase<Scalar> > domain_;
  Array<RCP<VectorBase<Scalar> > > col_vecs_;
  
};


} // end namespace Thyra


#endif // THYRA_DEFAULT_COLUMNWISE_MULTI_VECTOR_DECL_HPP
