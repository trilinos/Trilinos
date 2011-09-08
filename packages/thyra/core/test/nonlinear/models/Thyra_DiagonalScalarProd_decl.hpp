/*
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef THYRA_DIAGONAL_SCALAR_PROD_DECL_HPP
#define THYRA_DIAGONAL_SCALAR_PROD_DECL_HPP


#include "Thyra_ScalarProdBase.hpp"


namespace Thyra {


/** \brief Concrete implementation of a scalar product using a diagonal
 * vector.
 *
 * This test class really shows how to create an application-defined scalar
 * product.
 */
template<class Scalar>
class DiagonalScalarProd : public ScalarProdBase<Scalar> {
public:
  
  /** @name Consturctors/Initializers/Accessors */
  //@{

  /** \brief . */
  DiagonalScalarProd();

  /** \brief . */
  void initialize( const RCP<const VectorBase<Scalar> > &s_diag );

  //@}

protected:
  
  /** @name Overridden protected virtual functions from ScalarProdBase */
  //@{

  /** \brief Returns <tt>false</tt>. */
  virtual bool isEuclideanImpl() const;
  
  /** \brief . */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out ) const;

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > getLinearOpImpl() const;

  //@}

private:

  RCP<const VectorBase<Scalar> > s_diag_;

};


/** \brief Nonmember constructor.
 *
 * \relates DiagonalScalarProd
 */
template<class Scalar>
RCP<DiagonalScalarProd<Scalar> >
diagonalScalarProd(const RCP<const VectorBase<Scalar> > &s_diag)
{
  const RCP<DiagonalScalarProd<Scalar> > scalarProd =
    Teuchos::rcp(new DiagonalScalarProd<Scalar>());
  scalarProd->initialize(s_diag);
  return scalarProd;
}



} // end namespace Thyra


#endif  // THYRA_DIAGONAL_SCALAR_PROD_DECL_HPP
