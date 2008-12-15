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
#include <Teuchos_Object.hpp>
#include <Teuchos_BLAS.hpp>

#include "Tpetra_VectorDecl.hpp"

namespace Tpetra {

  template<typename Ordinal,typename Scalar>
  Vector<Ordinal,Scalar>::Vector(const Map<Ordinal> &map, bool zeroOut) 
    : MultiVector<Ordinal,Scalar>(map,1,zeroOut)
  { 
    this->setLabel("Tpetra::Vector");
  }

  template<typename Ordinal, typename Scalar>
  Vector<Ordinal,Scalar>::Vector(const Vector<Ordinal,Scalar> &source)
  : MultiVector<Ordinal,Scalar>(source)
  {
    this->setLabel("Tpetra::Vector");
  }

  template<typename Ordinal,typename Scalar>
  Vector<Ordinal,Scalar>::Vector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &values)
  : MultiVector<Ordinal,Scalar>(map,values,values.size(),1)
  {
    this->setLabel("Tpetra::Vector");
  }

  template <typename Ordinal, typename Scalar> 
  Vector<Ordinal,Scalar>::Vector(const Map<Ordinal> &map, const Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > &mvdata) 
    : MultiVector<Ordinal,Scalar>(map,mvdata)
  {
    this->setLabel("Tpetra::Vector");
  }

  template<typename Ordinal, typename Scalar>
  Vector<Ordinal,Scalar>::~Vector() {}

  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::replaceGlobalValue(Ordinal globalRow, const Scalar &value) 
  {
    this->replaceGlobalValue(globalRow,0,value);
  }

  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::sumIntoGlobalValue(Ordinal globalRow, const Scalar &value) 
  {
    this->sumIntoGlobalValue(globalRow,0,value);
  }

  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::replaceMyValue(Ordinal myRow, const Scalar &value) 
  {
    this->replaceMyValue(myRow,0,value);
  }

  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::sumIntoMyValue(Ordinal myRow, const Scalar &value) 
  {
    this->sumIntoMyValue(myRow,0,value);
  }


  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::extractCopy(Teuchos::ArrayView<Scalar> A) const 
  {
    Ordinal dummy;
    this->extractCopy(A,dummy);
  }

  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::extractView(Teuchos::ArrayView<Scalar> &A) 
  {
    Ordinal dummy;
    this->extractView(A,dummy);
  }

  template<typename Ordinal,typename Scalar>
  void Vector<Ordinal,Scalar>::extractConstView(Teuchos::ArrayView<const Scalar> &A) const 
  {
    Ordinal dummy;
    this->extractConstView(A,dummy);
  }

  template<typename Ordinal,typename Scalar>
  Scalar Vector<Ordinal,Scalar>::dot(const Vector<Ordinal,Scalar> &a) const 
  {
    Teuchos::BLAS<int,Scalar> blas;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    // compute local dot products of *this and a
    // sum these across all nodes
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors must have compatible Maps.");
    Scalar ldot, gdot;
    ldot = blas.DOT(this->MVData_->cPtrs_[0].size(),this->MVData_->cPtrs_[0].getRawPtr(),ONE,&a[0],ONE);
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,ldot,&gdot);
    }
    else {
      gdot = ldot;
    }
    return gdot;
  }

  template<typename Ordinal,typename Scalar>
  Scalar Vector<Ordinal,Scalar>::meanValue() const 
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::as;
    typedef ScalarTraits<Scalar> SCT;
    const Ordinal numImages = this->getMap().getComm()->getSize();
    // compute local components of the means
    // sum these across all nodes
    Array<Scalar> lmean = ScalarTraits<Scalar>::zero(),
                  gmean;
    typename ArrayView<const Scalar>::iterator cpos = this->MVData_->cPtrs_[0].begin(),
                                               cend = this->MVData_->cPtrs_[0].end();
    for (; cpos != cend; ++cpos) {
      lmean += *cpos;
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lmean,&gmean);
    }
    else {
      gmean = lmean;
    }
    return gmean/as<Scalar>(numImages);
  }

  template<typename Ordinal,typename Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Ordinal,Scalar>::norm1() const
  {
    Teuchos::BLAS<int,Scalar> blas;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    Mag lnorm, gnorm;
    lnorm = blas.ASUM(this->MVData_->cPtrs_[0].size(),this->MVData_->cPtrs_[0].getRawPtr(),ONE);
    if (this->isDistributed()) {
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,lnorm,&gnorm);
    }
    else {
      gnorm = lnorm;
    }
    return gnorm;
  }

  template<typename Ordinal,typename Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Ordinal,Scalar>::norm2() const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    Mag lnorm = ScalarTraits<Mag>::zero(),
        gnorm;
    typename Teuchos::ArrayView<const Scalar>::iterator cpos = this->MVData_->cPtrs_[0].begin(),
                                                        cend = this->MVData_->cPtrs_[0].end();
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

  template<typename Ordinal,typename Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Ordinal,Scalar>::normInf() const
  {
    Teuchos::BLAS<int,Scalar> blas;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    Mag lnorm, gnorm;
    // careful!!! IAMAX returns FORTRAN-style (i.e., one-based) index. subtract ind by one
    Ordinal ind = blas.IAMAX(this->MVData_->cPtrs_[0].size(),this->MVData_->cPtrs_[0].getRawPtr(),ONE) - ONE;
    lnorm = Teuchos::ScalarTraits<Scalar>::magnitude( this->MVData_->cPtrs_[0][ind] );
    if (this->isDistributed()) {
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,lnorm,&gnorm);
    }
    else {
      gnorm = lnorm;
    }
    return gnorm;
  }

  template<typename Ordinal,typename Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType Vector<Ordinal,Scalar>::normWeighted(const Vector<Ordinal,Scalar> &weights) const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::as;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;
    const Ordinal numImages = this->getMap().getComm()->getSize();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
        "Tpetra::Vector::normWeighted(): Vectors must have compatible Maps.");
    // compute local components of the norms
    // sum these across all nodes
    Mag lnorm = ScalarTraits<Mag>::zero(),
        gnorm;
    typename ArrayView<const Scalar>::iterator 
      wpos = weights.MVData_->cPtrs_[0].begin(),
      cpos =   this->MVData_->cPtrs_[0].begin(),
      cend =   this->MVData_->cPtrs_[0].end();
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
    return ScalarTraits<Mag>::squareroot(gnorm/as<Mag>(numImages));
  }

  template <typename Ordinal, typename Scalar> 
  void Vector<Ordinal,Scalar>::print(std::ostream &os) const
  {
    using std::endl;
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm = this->getMap().getComm();
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


  template <typename Ordinal, typename Scalar> 
  void Vector<Ordinal,Scalar>::printValues(std::ostream &os) const
  {
    using std::endl;
    using std::setw;
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm = this->getMap().getComm();
    const Map<Ordinal> &map = this->getMap();
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

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::abs(const Vector<Ordinal,Scalar> &a) 
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::abs(): Vectors must have compatible Maps.");
    typename ArrayView<Scalar>::iterator curpos = this->MVData_->ptrs_[0].begin(),
                                           apos =     a.MVData_->ptrs_[0].begin();
    for (; curpos != this->MVData_->ptrs_[0].end(); ++curpos, ++apos) {
      *curpos = Teuchos::ScalarTraits<Scalar>::magnitude(*apos);
    }
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::reciprocal(const Vector<Ordinal,Scalar> &a) 
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::reciprocal(): Vectors must have compatible Maps.");
    typename Teuchos::ArrayView<Scalar>::iterator curpos = this->MVData_->ptrs_[0].begin(),
                                                    apos =     a.MVData_->ptrs_[0].begin();
    for (; curpos != this->MVData_->ptrs_[0].end(); ++curpos, ++apos) {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION( ScalarTraits<Scalar>::magnitude(*apos) <= ScalarTraits<Scalar>::sfmin() ||
          *apos == ScalarTraits<Scalar>::sfmin(), std::runtime_error,
          "Tpetra::Vector::reciprocal(): element of a was zero or too small to invert: " << *apos );
#endif
      *curpos = ScalarTraits<Scalar>::one()/(*apos);
    }
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::scale(const Scalar &alpha, const Vector<Ordinal,Scalar> &a) 
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::scale(): Vectors must have compatible Maps.");
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    if (alpha == Teuchos::ScalarTraits<Scalar>::zero()) {
      putScalar(alpha); // set me = 0.0
    }
    else if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      *this = a;        // set me = a
    }
    else {
      // set me == alpha*a
      ArrayView<Scalar> &curpos = this->MVData_->ptrs_[0],
                        &apos   =     a.MVData_->ptrs_[0];
      // copy a to *this
      blas.COPY(curpos.size(),apos.getRawPtr(),ONE,curpos.getRawPtr(),ONE);
      // then scale *this in-situ
      blas.SCAL(curpos.size(),alpha,curpos.getRawPtr(),ONE);
    }
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::update(const Scalar &alpha, const Vector<Ordinal,Scalar> &a, const Scalar &beta) 
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename Teuchos::ArrayView<Scalar>::iterator avi;
    typedef typename Teuchos::ArrayView<const Scalar>::iterator avci;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    if (alpha == ST::zero()) {
      scale(beta);
      return;
    }
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()), std::runtime_error,
        "Tpetra::Vector::update(): MultiVectors must have compatible Maps.");

    if (beta == ST::zero()) { // this = alpha*a
      scale(alpha,a);
      return;
    }
    else if (beta == ST::one()) { // this = this + alpha*a
      if (alpha == ST::one()) {   // this = this + a
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = a.MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos) { *curpos = (*curpos) + (*Apos); }
      }
      else {                      // this = this + alpha*a
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = a.MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos) { *curpos = (*curpos) + alpha*(*Apos); }
      }
    }
    else {                        // this = beta*this + alpha*a
      if (alpha == ST::one()) {   // this = beta*this + a
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = a.MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos) { *curpos = beta*(*curpos) + (*Apos); }
      }
      else {                      // this = beta*this + alpha*a
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = a.MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos) { *curpos = beta*(*curpos) + alpha*(*Apos); }
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::update(const Scalar &alpha, const Vector<Ordinal,Scalar> &a, const Scalar &beta, const Vector<Ordinal,Scalar> &b, const Scalar &gamma)
  {
    typedef typename Teuchos::ArrayView<Scalar>::iterator avi;
    typedef typename Teuchos::ArrayView<const Scalar>::iterator avci;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    if (alpha == ST::zero()) {
      update(beta,b,gamma);
      return;
    }
    else if (beta == ST::zero()) {
      update(alpha,a,gamma);
      return;
    }
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(a.getMap()) || !this->getMap().isCompatible(b.getMap()),
        std::runtime_error, "Tpetra::Vector::update(): MultiVectors must have compatible Maps.");
    // determine if alpha==1 xor beta==1
    // if only one of these is 1.0, rename it alpha
    Teuchos::Ptr<const Vector<Ordinal,Scalar> > Aptr = Teuchos::ptr(&a), Bptr = Teuchos::ptr(&b);
    Teuchos::Ptr<const Scalar> lalpha = Teuchos::ptr(&alpha),
                               lbeta  = Teuchos::ptr(&beta);
    if (alpha!=ST::one() && beta==ST::one()) {
      // switch them
      Aptr = Teuchos::ptr(&b);
      Bptr = Teuchos::ptr(&a);
      lalpha = Teuchos::ptr(&beta);
      lbeta  = Teuchos::ptr(&alpha);
    }

    if (gamma == ST::zero()) { // this = lalpha*a + lbeta*b
      if (*lalpha == ST::one()) {
        if (*lbeta == ST::one()) { // this = gamma*this + a + b
          avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*Apos) + (*Bpos); }
        }
        else { // this = a + lbeta*b
          avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*Apos) + (*lbeta)*(*Bpos); }
        }
      }
      else { // this = lalpha*a + lbeta*b
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
      }
    }
    else if (gamma == ST::one()) { // this = this + lalpha*a + lbeta*b
      if ((*lalpha) == ST::one()) {
        if ((*lbeta) == ST::one()) { // this = this + a + b
          avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*Apos) + (*Bpos); }
        }
        else { // this = this + a + lbeta*b
          avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*Apos) + (*lbeta)*(*Bpos); }
        }
      }
      else { // this = this + lalpha*a + lbeta*b
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
      }
    }
    else { // this = gamma*this + lalpha*a + lbeta*b
      if ((*lalpha) == ST::one()) {
        if ((*lbeta) == ST::one()) { // this = gamma*this + a + b
          avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*Apos) + (*Bpos); }
        }
        else { // this = gamma*this + a + lbeta*b
          avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*Apos) + (*lbeta)*(*Bpos); }
        }
      }
      else { // this = gamma*this + lalpha*a + lbeta*b
        avi curpos = this->MVData_->ptrs_[0].begin(), cend = this->MVData_->ptrs_[0].end(); avci Apos = Aptr->MVData_->cPtrs_[0].begin(), Bpos = Bptr->MVData_->cPtrs_[0].begin();
        for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
      }
    }
  }


  template<typename Ordinal,typename Scalar>
  Scalar& Vector<Ordinal,Scalar>::operator[](Ordinal index)
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return this->MVData_->ptrs_[0][index];
  }


  template<typename Ordinal,typename Scalar>
  const Scalar& Vector<Ordinal,Scalar>::operator[](Ordinal index) const
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return this->MVData_->cPtrs_[0][index];
  }


} // namespace Tpetra

#endif // TPETRA_VECTOR_HPP
