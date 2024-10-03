// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef IFPACK2_BORDEREDOPERATOR_DEF_HPP
#define IFPACK2_BORDEREDOPERATOR_DEF_HPP

#include "Ifpack2_BorderedOperator_decl.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_TestForException.hpp"

namespace Ifpack2 {

template< class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
BorderedOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node >::
BorderedOperator (const Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node > >& A) : 
  A_ (A)
{ 
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, 
    Teuchos::typeName (*this) << "::BorderedOperator constructor: "
    "The input Operator A is null.");
}

template< class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
BorderedOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node >::getDomainMap() const
{
  return A_->getDomainMap();
}

template< class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
BorderedOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node >::getRangeMap() const
{
  return A_->getRangeMap();
}

template< class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
bool
BorderedOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node >::hasTransposeApply() const 
{
  return A_->hasTransposeApply();
}

template< class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
void 
BorderedOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node >::
apply (const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >& X,
       Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >& Y,
       Teuchos::ETransp mode, 
       Scalar coefAx, 
       Scalar coefY ) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::BorderedOperator::apply(): X.getNumVectors() = " 
    << X.getNumVectors() << " != Y.getNumVectors() = " 
    << Y.getNumVectors() << ".");
  A_->apply (X, Y, mode, coefAx, coefY );
}

} // namespace Ifpack2

#define IFPACK2_BORDEREDOPERATOR_INSTANT(S,LO,GO,N) \
  template class Ifpack2::BorderedOperator< S, LO, GO, N >;

#endif /* IFPACK2_BorderedOperator_DEF_HPP */
