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

#ifndef THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP
#define THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP

#include "Thyra_TpetraEuclideanScalarProd_decl.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Thyra_TpetraVector.hpp"


namespace Thyra {


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X,
  const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar>& scalarProds_out
  ) const
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  Teuchos::RCP<const TMV> X_tpetra = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(X));
  Teuchos::RCP<const TMV> Y_tpetra = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(Y));

  if (nonnull(X_tpetra) && nonnull(Y_tpetra)) {
    // Which one do we want transposed?
    // Tpetra transposes the argument of dot.
    // Below is the order from TpetraVectorSpace::scalarProdsImpl,
    // so this would transpose Y. However, Thyra::dots (which calls
    // the RTOp) transposes the first argument, so scalarProdsImpl
    // in EuclideanScalarProd transposes X...
    X_tpetra->dot(*Y_tpetra, scalarProds_out);
  } else {
    // If one of the casts succeeded, sync that MV to host space
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    if (nonnull(X_tpetra))
      Teuchos::rcp_const_cast<TMV>(X_tpetra)->template sync<Kokkos::HostSpace>();
    if (nonnull(Y_tpetra))
      Teuchos::rcp_const_cast<TMV>(Y_tpetra)->template sync<Kokkos::HostSpace>();
#else
    if (nonnull(X_tpetra))
      Teuchos::rcp_const_cast<TMV>(X_tpetra)->sync_host ();
    if (nonnull(Y_tpetra))
      Teuchos::rcp_const_cast<TMV>(Y_tpetra)->sync_host ();
#endif

    EuclideanScalarProd<Scalar>::scalarProdsImpl(X, Y, scalarProds_out);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

  RCP<const TMV> tmv = rcp_dynamic_cast<const TMV>(mv);
  if (nonnull(tmv)) {
    return tmv->getConstTpetraMultiVector();
  }

  RCP<const TV> tv = rcp_dynamic_cast<const TV>(mv);
  if (nonnull(tv)) {
    return tv->getConstTpetraVector();
  }

  return Teuchos::null;
}


} // end namespace Thyra


#endif  // THYRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP
