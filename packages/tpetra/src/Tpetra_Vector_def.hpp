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

#ifndef TPETRA_VECTOR_DEF_HPP
#define TPETRA_VECTOR_DEF_HPP

#include <Kokkos_NodeTrace.hpp>

#include "Tpetra_MultiVector.hpp"

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Vector_decl.hpp"
#endif

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, bool zeroOut) 
    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,1,zeroOut) {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source)
  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(source) {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Scalar> &values)
  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,values,values.size(),1) {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Vector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, Teuchos::ArrayRCP<Scalar> values)
    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,values,map->getNodeNumElements(),1) {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~Vector() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) {
    this->MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) {
    this->MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) {
    this->MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue(myRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) {
    this->MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue(myRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy(Teuchos::ArrayView<Scalar> A) const {
    size_t lda = this->getLocalLength();
    this->get1dCopy(A,lda);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const {
    using Teuchos::outArg;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap()->isCompatible(*a.getMap()), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap() 
        << "a.getMap(): " << std::endl << *a.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( this->getLocalLength() != a.getLocalLength(), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors do not have the same local length.");
#endif
    Scalar gbldot;
    gbldot = MVT::Dot(this->lclMV_,a.lclMV_);
    if (this->isDistributed()) {
      Scalar lcldot = gbldot;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,lcldot,outArg(gbldot));
    }
    return gbldot;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue() const {
    using Teuchos::outArg;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    Scalar sum = MVT::Sum(this->lclMV_);
    if (this->isDistributed()) {
      Scalar lsum = sum;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,lsum,outArg(sum));
    }
    return sum / Teuchos::as<Scalar>(this->getGlobalLength());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1() const {
    using Teuchos::outArg;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Mag norm = MVT::Norm1(this->lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,lnorm,outArg(norm));
    }
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2() const {
    using Teuchos::ScalarTraits;
    using Teuchos::outArg;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    Mag norm = MVT::Norm2Squared(this->lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,lnorm,outArg(norm));
    }
    return ScalarTraits<Mag>::squareroot(norm);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf() const {
    using Teuchos::outArg;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Mag norm = MVT::NormInf(this->lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_MAX,lnorm,outArg(norm));
    }
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const {
    using Teuchos::ScalarTraits;
    using Teuchos::outArg;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap()->isCompatible(*weights.getMap()), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap() 
        << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( this->getLocalLength() != weights.getLocalLength(), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors do not have the same local length.");
#endif
    Mag norm = MVT::WeightedNorm(this->lclMV_,weights.lclMV_);
    if (this->isDistributed()) {
      Mag lnorm = norm;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,lnorm,outArg(norm));
    }
    return ScalarTraits<Mag>::squareroot(norm / Teuchos::as<Mag>(this->getGlobalLength()));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<this->getGlobalLength()
        << "}";
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<this->getGlobalLength(); dec *= 10) {
      ++width;
    }
    Teuchos::OSTab tab(out);
    if (vl != VERB_NONE) {
      // VERB_LOW and higher prints description()
      if (myImageID == 0) out << this->description() << std::endl; 
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // VERB_MEDIUM and higher prints getLocalLength()
            out << "node " << setw(width) << myImageID << ": local length=" << this->getLocalLength() << endl;
            if (vl != VERB_MEDIUM) {
              // VERB_HIGH and higher prints isConstantStride() and stride()
              if (vl == VERB_EXTREME && this->getLocalLength() > 0) {
                Teuchos::RCP<Node> node = this->lclMV_.getNode();
                KOKKOS_NODE_TRACE("Vector::describe()")
                Teuchos::ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(
                                                               this->getLocalLength(), 
                                                               MVT::getValues(this->lclMV_) );
                // VERB_EXTREME prints values
                for (size_t i=0; i<this->getLocalLength(); ++i) {
                  out << setw(width) << this->getMap()->getGlobalElement(i) 
                      << ": "
                      << myview[i] << endl;
                }
                myview = Teuchos::null;
              }
            }
            else {
              out << endl;
            }
          }
        }
      }
    }
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_VECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class Vector< SCALAR , LO , GO , NODE >; \


#endif // TPETRA_VECTOR_DEF_HPP
