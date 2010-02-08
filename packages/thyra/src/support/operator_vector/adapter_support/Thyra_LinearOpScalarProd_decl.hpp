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

#ifndef THYRA_LINEAR_OP_SCALAR_PROD_DECL_HPP
#define THYRA_LINEAR_OP_SCALAR_PROD_DECL_HPP


#include "Thyra_ScalarProdBase_decl.hpp"


namespace Thyra {

/** \brief Concrete implementation of a scalar product using a symmetric
 * positive-definite linear operator.
 *
 * This subclass will work with any <tt>VectorBase</tt> or
 * <tt>MultiVectorBase</tt> implementation who's vector spaces are compatible
 * with the underlying symmetric positive-definite linear operator object.
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class LinearOpScalarProd : public ScalarProdBase<Scalar> {
public:

  /** @name Constructors, initializers, accessors */
  //@{

  /** \brief . */
  LinearOpScalarProd();

  /** \brief . */
  LinearOpScalarProd( const RCP<const LinearOpBase<Scalar> > &op );

  /** \brief . */
  void initialize( const RCP<const LinearOpBase<Scalar> > &op );

  /** \brief . */
  const RCP<const LinearOpBase<Scalar> >& op() const;

  /** \brief . */
  void uninitialize(
    const Ptr<RCP<const LinearOpBase<Scalar> > > &op = Teuchos::null );

  //@}

protected:
  
  /** @name Overridden from ScalarProdBase */
  //@{

  /** \brief Returns <tt>false</tt>. */
  virtual bool isEuclideanImpl() const;

  /** \brief . */
  void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out
    ) const;

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > getLinearOpImpl() const;

  //@}

private:

  RCP<const LinearOpBase<Scalar> >  op_;

};


// //////////////////////////////////
// Inline members


template<class Scalar>
inline
const RCP<const LinearOpBase<Scalar> >& LinearOpScalarProd<Scalar>::op() const
{
  return op_;
}


} // end namespace Thyra


#endif  // THYRA_LINEAR_OP_SCALAR_PROD_DECL_HPP
