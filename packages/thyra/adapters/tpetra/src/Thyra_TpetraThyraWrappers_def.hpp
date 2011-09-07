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

#ifndef THYRA_TPETRA_THYRA_WRAPPERS_HPP
#define THYRA_TPETRA_THYRA_WRAPPERS_HPP


#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Thyra_TpetraLinearOp.hpp"


namespace Thyra {


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
getOrCreateTpetraVectorSpace(
  const RCP<const VectorSpaceBase<Scalar> > space,
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &tpetraMap
  )
{
  using Teuchos::rcp_dynamic_cast;
  typedef TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraVectorSpace_t;
  RCP<const TpetraVectorSpace_t> tpetraSpace;
  if (nonnull(space)) {
    tpetraSpace = rcp_dynamic_cast<const TpetraVectorSpace_t>(space, true);
  }
  else {
    tpetraSpace = tpetraVectorSpace<Scalar>(tpetraMap);
  }
  return tpetraSpace;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const ScalarProdVectorSpaceBase<Scalar> >
getOrCreateLocallyReplicatedTpetraVectorSpace(
  const RCP<const VectorSpaceBase<Scalar> > space,
  const RCP<const Teuchos::Comm<int> > &tpetraComm,
  const RCP<Node> &tpetraNode,
  const int numCols
  )
{
  using Teuchos::rcp_dynamic_cast;
  typedef TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraVectorSpace_t;
  RCP<const TpetraVectorSpace_t> tpetraSpace;
  if (nonnull(space)) {
    tpetraSpace = rcp_dynamic_cast<const TpetraVectorSpace_t>(space, true);
  }
  else {
    tpetraSpace = tpetraVectorSpace<Scalar>(
      Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal>(
        numCols, tpetraComm, tpetraNode 
        )
      );
  }
  return tpetraSpace;
}


} // namespace Thyra


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Thyra::createVectorSpace(
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &tpetraMap
  )
{
  return tpetraVectorSpace<Scalar>(tpetraMap);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
Thyra::createVector(
  const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector_in,
  const RCP<const VectorSpaceBase<Scalar> > space_in
  )
{
  return tpetraVector(
    getOrCreateTpetraVectorSpace(space_in, tpetraVector_in->getMap()),
    tpetraVector_in
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Thyra::createConstVector(
  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector_in,
  const RCP<const VectorSpaceBase<Scalar> > space
  )
{
  return constTpetraVector(
    getOrCreateTpetraVectorSpace(space, tpetraVector_in->getMap()),
    tpetraVector_in
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
Thyra::createMultiVector(
  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector_in,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace
  )
{
  return tpetraMultiVector(
    getOrCreateTpetraVectorSpace(rangeSpace, tpetraMultiVector_in->getMap()),
    getOrCreateLocallyReplicatedTpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
      domainSpace, tpetraMultiVector_in->getMap()->getComm(),
      tpetraMultiVector_in->getMap()->getNode(),
      tpetraMultiVector_in->getNumVectors()
      ),
    tpetraMultiVector_in
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
Thyra::createConstMultiVector(
  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector_in,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace
  )
{
  return constTpetraMultiVector(
    getOrCreateTpetraVectorSpace(rangeSpace, tpetraMultiVector_in->getMap()),
    getOrCreateLocallyReplicatedTpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
      domainSpace, tpetraMultiVector_in->getMap()->getComm(),
      tpetraMultiVector_in->getMap()->getNode(),
      tpetraMultiVector_in->getNumVectors()
      ),
    tpetraMultiVector_in
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::createLinearOp(
  const RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator_in,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace
  )
{
  return tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
    getOrCreateTpetraVectorSpace(rangeSpace, tpetraOperator_in->getRangeMap()),
    getOrCreateTpetraVectorSpace(domainSpace, tpetraOperator_in->getDomainMap()),
    tpetraOperator_in
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::createConstLinearOp(
  const RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator_in,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace
  )
{
  return constTpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
    getOrCreateTpetraVectorSpace(rangeSpace, tpetraOperator_in->getRangeMap()),
    getOrCreateTpetraVectorSpace(domainSpace, tpetraOperator_in->getDomainMap()),
    tpetraOperator_in
    );
}


namespace Thyra {


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetraVector(const RCP<VectorBase<Scalar> > &v)
{
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraVector_t;
  return Teuchos::rcp_dynamic_cast<TpetraVector_t>(v, true)->getTpetraVector();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraVector(const RCP<const VectorBase<Scalar> > &v)
{
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraVector_t;
  return Teuchos::rcp_dynamic_cast<const TpetraVector_t>(v, true)->getConstTpetraVector();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetraMultiVector(const RCP<MultiVectorBase<Scalar> > &mv)
{

#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(mv));
#endif

  using Teuchos::rcp_dynamic_cast;
  
  typedef Thyra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    ThyraTpetraMultiVector_t;
  const RCP<ThyraTpetraMultiVector_t> tmv =
    rcp_dynamic_cast<ThyraTpetraMultiVector_t>(mv);
  if (nonnull(tmv)) {
    return tmv->getTpetraMultiVector();
  }
  
  typedef Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    ThyraTpetraVector_t;
  const RCP<ThyraTpetraVector_t> tv =
    rcp_dynamic_cast<ThyraTpetraVector_t>(mv);
  if (nonnull(tv)) {
    return tv->getTpetraVector();
  }

  TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error, the input mv = " << mv->description() << " does not support the"
    " Thyra::TpetraMultiVector or the Thyra::TpetraVector interfaces!");

  return Teuchos::null;

}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> > &mv)
{

#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(mv));
#endif

  using Teuchos::rcp_dynamic_cast;
  
  typedef Thyra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    ThyraTpetraMultiVector_t;
  const RCP<const ThyraTpetraMultiVector_t> tmv =
    rcp_dynamic_cast<const ThyraTpetraMultiVector_t>(mv);
  if (nonnull(tmv)) {
    return tmv->getConstTpetraMultiVector();
  }
  
  typedef Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    ThyraTpetraVector_t;
  const RCP<const ThyraTpetraVector_t> tv =
    rcp_dynamic_cast<const ThyraTpetraVector_t>(mv);
  if (nonnull(tv)) {
    return tv->getConstTpetraVector();
  }

  TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error, the input mv = " << mv->description() << " does not support the"
    " Thyra::TpetraMultiVector or the Thyra::TpetraVector interfaces!");

  return Teuchos::null;

}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetraOperator(const RCP<LinearOpBase<Scalar> > &op)
{
  typedef TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinearOp_t;
  return Teuchos::rcp_dynamic_cast<TpetraLinearOp_t>(op, true)->getTpetraOperator();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraOperator(const RCP<const LinearOpBase<Scalar> > &op)
{
  typedef TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinearOp_t;
  return Teuchos::rcp_dynamic_cast<const TpetraLinearOp_t>(op, true)->getConstTpetraOperator();
}


} // namespace Thyra


#endif // THYRA_TPETRA_THYRA_WRAPPERS_HPP
