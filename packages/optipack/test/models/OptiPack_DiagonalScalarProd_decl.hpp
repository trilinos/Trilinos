/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef OPTIPACK_DIAGONAL_SCALAR_PROD_DECL_HPP
#define OPTIPACK_DIAGONAL_SCALAR_PROD_DECL_HPP


#include "OptiPack_Types.hpp"
#include "Thyra_ScalarProdBase.hpp"


namespace OptiPack {


/** \brief Concrete implementation of a scalar product using a diagonal
 * vector.
 *
 * This test class really shows how to create an application-defined scalar
 * product.
 */
template<class Scalar>
class DiagonalScalarProd : public Thyra::ScalarProdBase<Scalar> {
public:
  
  /** @name Consturctors/Initializers/Accessors */
  //@{

  /** \brief . */
  DiagonalScalarProd();

  /** \brief . */
  void initialize( const RCP<const Thyra::VectorBase<Scalar> > &s_diag );

  //@}

protected:
  
  /** @name Overridden protected virtual functions from ScalarProdBase */
  //@{

  /** \brief Returns <tt>false</tt>. */
  virtual bool isEuclideanImpl() const;
  
  /** \brief . */
  virtual void scalarProdsImpl(
    const Thyra::MultiVectorBase<Scalar>& X, const Thyra::MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out ) const;

  /** \brief . */
  RCP<const Thyra::LinearOpBase<Scalar> > getLinearOpImpl() const;

  //@}

private:

  RCP<const Thyra::VectorBase<Scalar> > s_diag_;

};


/** \brief Nonmember constructor. */
template<class Scalar>
RCP<DiagonalScalarProd<Scalar> >
diagonalScalarProd(const RCP<const Thyra::VectorBase<Scalar> > &s_diag)
{
  const RCP<DiagonalScalarProd<Scalar> > scalarProd =
    Teuchos::rcp(new DiagonalScalarProd<Scalar>());
  scalarProd->initialize(s_diag);
  return scalarProd;
}



} // end namespace OptiPack


#endif  // OPTIPACK_DIAGONAL_SCALAR_PROD_DECL_HPP
