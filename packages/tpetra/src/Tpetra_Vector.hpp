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

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, bool zeroOut) 
//    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(node,map,1,zeroOut)
//  {}

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source)
//  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(source)
//  {}

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &values)
//  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(node,map,values,values.size(),1)
//  {}

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> 
//  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::RCP<Kokkos::MultiVector<Scalar,Node> > &lclMV) 
//    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,lclMV)
//  {}

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~Vector() {}

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) 
//  {
//    this->replaceGlobalValue(globalRow,0,value);
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) 
//  {
//    this->sumIntoGlobalValue(globalRow,0,value);
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceMyValue(LocalOrdinal myRow, const Scalar &value) 
//  {
//    this->replaceMyValue(myRow,0,value);
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoMyValue(LocalOrdinal myRow, const Scalar &value) 
//  {
//    this->sumIntoMyValue(myRow,0,value);
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy(Teuchos::ArrayView<Scalar> A) const {
//    Teuchos_Ordinal lda = this->getMyLength();
//    this->get1dCopy(A,lda);
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy(Scalar *A) const {
//    this->get1dCopy( Teuchos::arrayView<Scalar>(A,this->getMyLength()) );
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const 
//  {
//#ifdef HAVE_TPETRA_DEBUG
//    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
//        "Tpetra::Vector::dots(): Vectors do not have compatible Maps:" << std::endl
//        << "this->getMap(): " << std::endl << this->getMap() 
//        << "a.getMap(): " << std::endl << a.getMap() << std::endl);
//#else
//    TEST_FOR_EXCEPTION( this->getMyLength() != a.getMyLength(), std::runtime_error,
//        "Tpetra::Vector::dots(): Vectors do not have the same local length.");
//#endif
//    Scalar dot;
//    dot = DMVA::Dot(*this->lclMV_,*a.lclMV_);
//    if (this->isDistributed()) {
//      Scalar lcl = dot;
//      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lcl,&dot);
//    }
//    return dot;
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue() const 
//  {
//    typedef Teuchos::ScalarTraits<Scalar> SCT;
//    Scalar sum = DMVA::Sum(*this->lclMV_);
//    if (this->isDistributed()) {
//      Scalar lsum = sum;
//      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lsum,&sum);
//    }
//    return sum / Teuchos::as<Scalar>(this->getGlobalLength());
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1() const
//  {
//    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
//    Mag norm = DMVA::Norm1(*this->lclMV_);
//    if (this->isDistributed()) {
//      Mag lnorm = norm;
//      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&norm);
//    }
//    return norm;
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2() const
//  {
//    using Teuchos::ScalarTraits;
//    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
//    Mag norm = DMVA::Norm2Squared(*this->lclMV_);
//    if (this->isDistributed()) {
//      Mag lnorm = norm;
//      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&norm);
//    }
//    return ScalarTraits<Mag>::squareroot(norm);
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf() const
//  {
//    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
//    Mag norm = DMVA::NormInf(*this->lclMV_);
//    if (this->isDistributed()) {
//      Mag lnorm = norm;
//      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,lnorm,&norm);
//    }
//    return norm;
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const
//  {
//    using Teuchos::ScalarTraits;
//    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
//    TEST_FOR_EXCEPTION(weights.getNumVectors() != 1, std::runtime_error,
//        "Tpetra::Vector::normWeighted(): Vector of weights must contain one vector.");
//#ifdef HAVE_TPETRA_DEBUG
//    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
//        "Tpetra::Vector::normWeighted(): Vectors do not have compatible Maps:" << std::endl
//        << "this->getMap(): " << std::endl << this->getMap() 
//        << "weights.getMap(): " << std::endl << weights.getMap() << std::endl);
//#else
//    TEST_FOR_EXCEPTION( this->getMyLength() != weights.getMyLength(), std::runtime_error,
//        "Tpetra::Vector::normWeighted(): Vectors do not have the same local length.");
//#endif
//    Mag norm = DMVA::WeightedNorm(*this->lclMV_,*weights.lclMV_);
//    if (this->isDistributed()) {
//      Mag lnorm = norm;
//      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&norm);
//    }
//    return ScalarTraits<Mag>::squareroot(norm / Teuchos::as<Mag>(this->getGlobalLength()));
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  std::string Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
//  {
//    std::ostringstream oss;
//    oss << Teuchos::Describable::description();
//    oss << "{length="<<this->getGlobalLength()
//        << "}";
//    return oss.str();
//  }

//  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
//  {
//    using std::endl;
//    using std::setw;
//    using Teuchos::VERB_DEFAULT;
//    using Teuchos::VERB_NONE;
//    using Teuchos::VERB_LOW;
//    using Teuchos::VERB_MEDIUM;
//    using Teuchos::VERB_HIGH;
//    using Teuchos::VERB_EXTREME;
//    Teuchos::EVerbosityLevel vl = verbLevel;
//    if (vl == VERB_DEFAULT) vl = VERB_LOW;
//    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap().getComm();
//    const int myImageID = comm->getRank();
//    const int numImages = comm->getSize();
//    int width = 1;
//    for (int dec=10; dec<this->getGlobalLength(); dec *= 10) {
//      ++width;
//    }
//    Teuchos::OSTab tab(out);
//    if (vl != VERB_NONE) {
//      // VERB_LOW and higher prints description()
//      if (myImageID == 0) out << this->description() << std::endl; 
//      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
//        if (myImageID == imageCtr) {
//          if (vl != VERB_LOW) {
//            // VERB_MEDIUM and higher prints getMyLength()
//            out << "node " << setw(width) << myImageID << ": local length=" << this->getMyLength() << endl;
//            if (vl != VERB_MEDIUM) {
//              // VERB_HIGH and higher prints isConstantStride() and stride()
//              if (vl == VERB_EXTREME && this->getMyLength() > 0) {
//                Node &node = this->lclMV_->getNode();
//                Teuchos::ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(
//                                                              this->getMyLength(), 
//                                                              this->lclMV_->getValues() );
//                // VERB_EXTREME prints values
//                for (Teuchos_Ordinal i=0; i<this->getMyLength(); ++i) {
//                  out << setw(width) << this->getMap().getGlobalIndex(i) 
//                      << ": "
//                      << myview[i] << endl;
//                }
//                myview = Teuchos::null;
//              }
//            }
//            else {
//              out << endl;
//            }
//          }
//        }
//      }
//    }
//  }

} // namespace Tpetra

#endif // TPETRA_VECTOR_HPP
