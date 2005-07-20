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

#ifndef THYRA_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_VectorBaseDecl.hpp"
#include "Thyra_MultiVectorDefaultBaseDecl.hpp"
#include "Thyra_SingleRhsLinearOpBaseDecl.hpp"

namespace Thyra {

/** \brief Convienent node subclass for concrete <tt>VectorBase</tt>
 * subclasses.
 *
 * This node subclass provides as many default implementations as possible for
 * virtual functions.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
class VectorDefaultBase
  : virtual public VectorBase<Scalar>
  , virtual protected MultiVectorDefaultBase<Scalar>
  , virtual protected SingleRhsLinearOpBase<Scalar>
{
public:

  /** \brief . */
  using SingleRhsLinearOpBase<Scalar>::apply;
  /** \brief . */
  using MultiVectorDefaultBase<Scalar>::describe;
  /** \brief . */
  using MultiVectorDefaultBase<Scalar>::applyOp;
  /** \brief . */
  using MultiVectorDefaultBase<Scalar>::col;

  /** @name Overridden from LinearOpBase (should never need to be overridden in subclasses) */
  //@{
  /// Returns <tt>this->space()</tt>
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  /// Returns a <tt>SerialVectorSpace</tt> object with dimension 1.
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  //@}

  /** @name Overridden from MultiVectorBase (should never need to be overridden in subclasses) */
  //@{
  /// Returns <tt>Teuchos::rcp(this,false)</tt>
  Teuchos::RefCountPtr<VectorBase<Scalar> > col(Index j);
  /// Returns <tt>this->clone_v()</tt>
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > clone_mv() const;
  /// Returns <tt>Teuchos::rcp(this,false)</tt>
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const Range1D& col_rng ) const;
  /// Returns <tt>Teuchos::rcp(this,false)</tt>
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const Range1D& col_rng );
  /// Returns <tt>Teuchos::rcp(this,false)</tt>
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] ) const;
  /// Returns <tt>Teuchos::rcp(this,false)</tt>
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] );
  /// Implemented in terms of <tt>this->getSubVector()</tt>
  void getSubMultiVector(
    const Range1D                       &rowRng
    ,const Range1D                      &colRng
    ,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
    ) const;
  /// Implemented in terms of <tt>this->freeSubVector()</tt>
  void freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* sub_mv ) const;
  /// Implemented in terms of <tt>this->getSubVector()</tt>
  void getSubMultiVector(
    const Range1D                                &rowRng
    ,const Range1D                               &colRng
    ,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
    );
  /// Implemented in terms of <tt>this->commitSubVector()</tt>
  void commitSubMultiVector( RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv );
  //@}

protected:

  /** @name Overridden from SingleRhsLinearOpBase (should never need to be overridden in subclasses) */
  //@{

  /** \brief For complex <tt>Scalar</tt> types returns <tt>true</tt> for
   * <tt>NOTRANS</tt> and <tt>CONJTRANS</tt> and for real types returns true
   * for all values of <tt>M_trans</tt>.
   */
  bool opSupported(ETransp M_trans) const;

  /** \brief. Applies vector and its adjoint (transpose) as a linear operator. */
  void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private data members

  Teuchos::RefCountPtr<VectorSpaceBase<Scalar> >  domain_; // Only initialized if *this is used as a MultiVectorBase

  // /////////////////////////////////////
  // Private member functions

  void validateColRng( const Range1D &rowRng ) const;
  void validateColIndexes(  const int numCols, const int cols[] ) const;

};

} // end namespace Thyra

#endif  // THYRA_VECTOR_DEFAULT_BASE_DECL_HPP
