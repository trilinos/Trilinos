// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DIAGONAL_SCALAR_PROD_DEF_HPP
#define THYRA_DIAGONAL_SCALAR_PROD_DEF_HPP


#include "Thyra_DiagonalScalarProd_decl.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


// Consturctors/Initializers/Accessors


template<class Scalar>
DiagonalScalarProd<Scalar>::DiagonalScalarProd()
{}


template<class Scalar>
void DiagonalScalarProd<Scalar>::initialize(
  const RCP<const VectorBase<Scalar> > &s_diag )
{
  s_diag_ = s_diag.assert_not_null();
}


// Overridden from ScalarProdBase


template<class Scalar>
bool DiagonalScalarProd<Scalar>::isEuclideanImpl() const
{
  return false;
}


template<class Scalar>
void DiagonalScalarProd<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out ) const
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  const Ordinal m = X.domain()->dim();

#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_VEC_SPACES( "DiagonalScalarProd<Scalar>::scalarProds(X,Y,sclarProds)",
    *s_diag_->space(), *Y.range() );
  THYRA_ASSERT_VEC_SPACES( "DiagonalScalarProd<Scalar>::scalarProds(X,Y,sclarProds)",
    *X.range(), *Y.range() );
  THYRA_ASSERT_VEC_SPACES( "DiagonalScalarProd<Scalar>::scalarProds(X,Y,sclarProds)",
    *X.domain(), *Y.domain() );
  TEUCHOS_ASSERT_EQUALITY( as<Ordinal>(scalarProds_out.size()), m );
#endif

  const ConstDetachedSpmdVectorView<Scalar> s_diag(s_diag_);

  const RCP<const Teuchos::Comm<Ordinal> > comm = s_diag.spmdSpace()->getComm();
  
  for (Ordinal j = 0; j < m; ++j) {

    const ConstDetachedSpmdVectorView<Scalar> x(X.col(j));
    const ConstDetachedSpmdVectorView<Scalar> y(Y.col(j));
    
    Scalar scalarProd_j = ST::zero();

    for (Ordinal i = 0; i < x.subDim(); ++i) {
      scalarProd_j += ST::conjugate(x[i]) * s_diag[i] * y[i];
    }

    if (!is_null(comm)) {
      Scalar g_scalarProd_j = 0.0;
      Teuchos::reduceAll<Ordinal,Scalar>(
        *comm, Teuchos::REDUCE_SUM,
        scalarProd_j,
        Teuchos::outArg(g_scalarProd_j)
        );
      scalarProds_out[j] = g_scalarProd_j;
    }
    else {
      scalarProds_out[j] = scalarProd_j;
    }

  }

}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DiagonalScalarProd<Scalar>::getLinearOpImpl() const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "ToDo: Implement when needed!")
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}


} // end namespace Thyra


#endif  // THYRA_DIAGONAL_SCALAR_PROD_DEF_HPP
