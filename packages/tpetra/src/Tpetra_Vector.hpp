// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_VECTOR_HPP
#define TPETRA_VECTOR_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_VectorDecl.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, bool zeroOut) 
    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(node,map,1,zeroOut)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source)
  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(source)
  { }

  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &values)
  //REFACTOR// : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,values,values.size(),1)
  //REFACTOR// { }

  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayRCP<Scalar> &values)
  //REFACTOR// : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,values,values.size(),1)
  //REFACTOR// { }

  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> 
  //REFACTOR// Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::RCP<MultiVectorData<Scalar> > &mvdata) 
  //REFACTOR//   : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,mvdata)
  //REFACTOR// { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~Vector() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) 
  {
    this->replaceGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) 
  {
    this->sumIntoGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceMyValue(LocalOrdinal myRow, const Scalar &value) 
  {
    this->replaceMyValue(myRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoMyValue(LocalOrdinal myRow, const Scalar &value) 
  {
    this->sumIntoMyValue(myRow,0,value);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractCopy1D(Teuchos::ArrayView<Scalar> A) const 
  {
    Teuchos_Ordinal lda = this->myLength();
    this->extractCopy1D(A,lda);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractCopy1D(Scalar *A) const 
  {
    Teuchos_Ordinal lda = this->myLength();
    this->extractCopy1D(A,lda);
  }

  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractView1D(Teuchos::ArrayView<Scalar> &A) 
  //REFACTOR// {
  //REFACTOR//   Teuchos_Ordinal dummy;
  //REFACTOR//   this->extractView1D(A,dummy);
  //REFACTOR// }

  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractConstView1D(Teuchos::ArrayView<const Scalar> &A) const 
  //REFACTOR// {
  //REFACTOR//   Teuchos_Ordinal dummy;
  //REFACTOR//   this->extractConstView1D(A,dummy);
  //REFACTOR// }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const 
  {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "a.getMap(): " << std::endl << a.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( this->myLength() != a.myLength(), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors do not have the same local length.");
#endif
    Scalar dot;
    dot = DMVA::Dot(*this->lclMV_,*a.lclMV_);
    if (this->isDistributed()) {
      Scalar lcl = dot;
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lcl,&dot);
    }
    return dot;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue() const 
  {
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    Scalar sum = DMVA::Sum(*this->lclMV_);
    if (this->isDistributed()) {
      Scalar lsum = sum;
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lsum,&sum);
    }
    return sum / Teuchos::as<Scalar>(this->globalLength());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1() const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Mag norm = DMVA::Norm1(*this->lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&norm);
    }
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2() const
  {
    using Teuchos::ScalarTraits;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    Mag norm = DMVA::Norm2Squared(*this->lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&norm);
    }
    return ScalarTraits<Mag>::squareroot(norm);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf() const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Mag norm = DMVA::NormInf(*this->lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,lnorm,&norm);
    }
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const
  {
    using Teuchos::ScalarTraits;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    TEST_FOR_EXCEPTION(weights.numVectors() != 1, std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vector of weights must contain one vector.");
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "weights.getMap(): " << std::endl << weights.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( this->myLength() != weights.myLength(), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors do not have the same local length.");
#endif
    Mag norm = DMVA::WeightedNorm(*this->lclMV_,*weights.lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&norm);
    }
    return ScalarTraits<Mag>::squareroot(norm / Teuchos::as<Mag>(this->globalLength()));
  }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Scalar& Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator[](Teuchos_Ordinal index)
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(index < 0 || index >= myLength(), std::runtime_error,
  //REFACTOR//       "Tpetra::Vector::operator[](j): index j exceeds local dimension.");
  //REFACTOR//   Node &node = lclMV_->getNode();
  //REFACTOR//   Scalar * mvdata = node.template viewBuffer<Scalar>(myLength(),lclMV_->getValues(0),0);
  //REFACTOR//   return mv; // FINISH: can't return reference and release view.... what to do?
  //REFACTOR//   node.template releaseView<Scalar>(mvdata);
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// const Scalar& Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator[](Teuchos_Ordinal index) const
  //REFACTOR// {
  //REFACTOR//   ...
  //REFACTOR//   see operator[] non-const above
  //REFACTOR//   ...
  //REFACTOR// }


} // namespace Tpetra

#endif // TPETRA_VECTOR_HPP
