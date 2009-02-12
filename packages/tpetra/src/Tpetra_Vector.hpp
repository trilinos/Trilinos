// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef TPETRA_VECTOR_HPP
#define TPETRA_VECTOR_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_BLAS.hpp>

#include "Tpetra_VectorDecl.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, bool zeroOut) 
    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(map,1,zeroOut)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal>::Vector(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &source)
  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(source)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &values)
  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(map,values,values.size(),1)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayRCP<Scalar> &values)
  : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(map,values,values.size(),1)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal> 
  Vector<Scalar,LocalOrdinal,GlobalOrdinal>::Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::RCP<MultiVectorData<Scalar> > &mvdata) 
    : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(map,mvdata)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal>::~Vector() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) 
  {
    this->replaceGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) 
  {
    this->sumIntoGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::replaceMyValue(LocalOrdinal myRow, const Scalar &value) 
  {
    this->replaceMyValue(myRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::sumIntoMyValue(LocalOrdinal myRow, const Scalar &value) 
  {
    this->sumIntoMyValue(myRow,0,value);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::extractCopy1D(Teuchos::ArrayView<Scalar> A) const 
  {
    Teuchos_Ordinal dummy;
    this->extractCopy1D(A,dummy);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::extractCopy1D(Scalar *A) const 
  {
    Teuchos_Ordinal dummy;
    this->extractCopy1D(A,dummy);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::extractView1D(Teuchos::ArrayView<Scalar> &A) 
  {
    Teuchos_Ordinal dummy;
    this->extractView1D(A,dummy);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::extractConstView1D(Teuchos::ArrayView<const Scalar> &A) const 
  {
    Teuchos_Ordinal dummy;
    this->extractConstView1D(A,dummy);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal>::dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &a) const 
  {
    Teuchos::BLAS<int,Scalar> blas;
    // compute local dot products of *this and a
    // sum these across all nodes
#ifdef TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::dot(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "a.getMap(): " << std::endl << a.getMap() << std::endl);
#else 
    TEST_FOR_EXCEPTION( this->myLength() != a.myLength(), std::runtime_error,
        "Tpetra::Vector::dot(): Vectors do not have the same local length.");
#endif
    Scalar ldot, gdot;
    ldot = blas.DOT(this->myLength(),&(*this)[0],1,&a[0],1);
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,ldot,&gdot);
    }
    else {
      gdot = ldot;
    }
    return gdot;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal>::meanValue() const 
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    // compute local components of the means
    // sum these across all nodes
    Scalar lmean = ST::zero(),
           gmean;
    const_pointer cpos = this->MVData_->ptrs_[0], 
                  cend = this->MVData_->ptrs_[0]+this->myLength();
    for (; cpos != cend; ++cpos) {
      lmean += (*cpos);
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lmean,&gmean);
    }
    else {
      gmean = lmean;
    }
    return gmean/Teuchos::as<Scalar>(this->globalLength());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal>::norm1() const
  {
    Teuchos::BLAS<int,Scalar> blas;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Mag lnorm, gnorm;
    lnorm = blas.ASUM(this->myLength(),&(*this->MVData_->ptrs_[0]),1);
    if (this->isDistributed()) {
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&gnorm);
    }
    else {
      gnorm = lnorm;
    }
    return gnorm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal>::norm2() const
  {
    using Teuchos::ScalarTraits;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    Mag lnorm = ScalarTraits<Mag>::zero(),
        gnorm;
    const_pointer cpos = this->MVData_->ptrs_[0], 
                  cend = this->MVData_->ptrs_[0]+this->myLength();
    for (; cpos != cend; ++cpos) {
      lnorm += ScalarTraits<Scalar>::magnitude( 
          (*cpos) * ScalarTraits<Scalar>::conjugate(*cpos)
        );
    }
    if (this->isDistributed()) {
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&gnorm);
    }
    else {
      gnorm = lnorm;
    }
    return ScalarTraits<Mag>::squareroot(gnorm);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal>::normInf() const
  {
    Teuchos::BLAS<int,Scalar> blas;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Mag lnorm = Teuchos::ScalarTraits<Mag>::zero(), gnorm;
    // careful!!! IAMAX returns FORTRAN-style (i.e., one-based) index. subtract ind by one
    if (this->myLength() > 0) {
      int ind = blas.IAMAX(this->myLength(),&(*this->MVData_->ptrs_[0]),1) - 1;
      lnorm = Teuchos::ScalarTraits<Scalar>::magnitude( (*this)[ind] );
    }
    if (this->isDistributed()) {
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,lnorm,&gnorm);
    }
    else {
      gnorm = lnorm;
    }
    return gnorm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Scalar,LocalOrdinal,GlobalOrdinal>::normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &weights) const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::as;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;
#ifdef TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "weights.getMap(): " << std::endl << weights.getMap() << std::endl);
#else 
    TEST_FOR_EXCEPTION( this->myLength() != weights.myLength(), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors do not have the same local length.");
#endif
    // compute local components of the norms
    // sum these across all nodes
    Mag lnorm = ScalarTraits<Mag>::zero(),
        gnorm;
    const_pointer
      wpos = weights.MVData_->ptrs_[0],
      cpos =   this->MVData_->ptrs_[0],
      cend =   this->MVData_->ptrs_[0] + this->myLength();
    for (; cpos != cend; ++cpos, ++wpos) {
      Scalar tmp = *cpos / *wpos;
      lnorm += SCT::magnitude( tmp * SCT::conjugate(tmp) );
    }
    if (this->isDistributed()) {
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&gnorm);
    }
    else {
      gnorm = lnorm;
    }
    return ScalarTraits<Mag>::squareroot(gnorm/as<Mag>(this->globalLength()));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal> 
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::print(std::ostream &os) const
  {
    using std::endl;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap().getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        if (myImageID == 0) {
          os << "Global length: " << this->globalLength() << endl;
        }
        os << "Local length: " << this->myLength() << endl;
      }
      // Do a few global ops to give I/O a chance to complete
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal> 
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal>::printValues(std::ostream &os) const
  {
    using std::endl;
    using std::setw;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap().getComm();
    const Map<LocalOrdinal,GlobalOrdinal> &map = this->getMap();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        for (int i=0; i<this->myLength(); ++i) {
          os << setw(4) << map.getGlobalIndex(i) << "\t";
          os << setw(20) << this->MVData_->ptrs_[0][i] << " ";
          os << endl;
        }
      }
      // Do a few global ops to give I/O a chance to complete
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Scalar& Vector<Scalar,LocalOrdinal,GlobalOrdinal>::operator[](Teuchos_Ordinal index)
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return this->MVData_->ptrs_[0][index];
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Scalar& Vector<Scalar,LocalOrdinal,GlobalOrdinal>::operator[](Teuchos_Ordinal index) const
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return this->MVData_->ptrs_[0][index];
  }


} // namespace Tpetra

#endif // TPETRA_VECTOR_HPP
