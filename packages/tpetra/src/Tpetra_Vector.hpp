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

#include "Tpetra_VectorDecl.hpp"
#include "Tpetra_VectorData.hpp"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#ifdef HAVE_TPETRA_TBB
#include "Tpetra_TBB_TaskScheduler.hpp"
#include "Tpetra_TBB_Vec_Kernels.hpp"
#endif

#include <Teuchos_Object.hpp>

namespace Tpetra {

  template<typename OrdinalType, typename ScalarType>
  Vector<OrdinalType,ScalarType>::Vector(VectorSpace<OrdinalType, ScalarType> const& VectorSpace) 
    : DistObject<OrdinalType, ScalarType>(VectorSpace.elementSpace(), 
        VectorSpace.platform().createScalarComm(), "Tpetra::Vector")
      , VectorData_()
  {
    ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
    //OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero(); // *WARNING-UNUSED*
    OrdinalType const length = VectorSpace.getNumMyEntries();

    VectorData_ = Teuchos::rcp(new VectorData<OrdinalType, ScalarType>(VectorSpace, length, scalarZero));
  }

  template<typename OrdinalType, typename ScalarType>
  Vector<OrdinalType,ScalarType>::Vector(
        const ScalarType* vectorEntries, OrdinalType numEntries
      , VectorSpace<OrdinalType, ScalarType> const& VectorSpace)
    : DistObject<OrdinalType, ScalarType>(VectorSpace.elementSpace()
      , VectorSpace.platform().createScalarComm(), "Tpetra::Vector")
      , VectorData_()
  {
    ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const length = VectorSpace.getNumMyEntries();

    VectorData_ = Teuchos::rcp(new VectorData<OrdinalType, ScalarType>(VectorSpace, length, scalarZero));

    if(numEntries != length)
      throw reportError("numEntries = " + toString(numEntries) + ".  Should be = " + toString(length) + ".", -1);
    for(OrdinalType i = ordinalZero; i < length; i++)
      VectorData_->scalarArray_[i] = vectorEntries[i];
  }

  template<typename OrdinalType, typename ScalarType>
  Vector<OrdinalType,ScalarType>::Vector(Vector<OrdinalType, ScalarType> const& Source)
    : DistObject<OrdinalType, ScalarType>(Source)
      , VectorData_(Source.VectorData_)
  {}

  template<typename OrdinalType, typename ScalarType>
  Vector<OrdinalType,ScalarType>::~Vector() {}

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::submitEntries(OrdinalType numEntries, OrdinalType const* indices, ScalarType const* values) {
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    for(OrdinalType i = ordinalZero; i < numEntries; i++)
      VectorData_->scalarArray_[indices[i]] += values[i];
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::setAllToScalar(ScalarType const value) {
    OrdinalType const max = getNumMyEntries();
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    for(OrdinalType i = ordinalZero; i < max; i++)
      VectorData_->scalarArray_[i] = value;
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::setAllToRandom() {
    OrdinalType const max = getNumMyEntries();
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    for(OrdinalType i = ordinalZero; i < max; i++)
      VectorData_->scalarArray_[i] = Teuchos::ScalarTraits<ScalarType>::random();
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::extractCopy(ScalarType* userArray) const {
    OrdinalType const max = getNumMyEntries();
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    for(OrdinalType i = ordinalZero; i < max; i++)
      userArray[i] = VectorData_->scalarArray_[i];
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::extractView(ScalarType** userPointerArray) const {
    OrdinalType const max = getNumMyEntries();
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    for(OrdinalType i = ordinalZero; i < max; i++)
      userPointerArray[i] = &VectorData_->scalarArray_[i];
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::dotProduct(Vector<OrdinalType, ScalarType> const& x) const {
    if(! vectorSpace().isCompatible(x.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType length = getNumMyEntries();
    OrdinalType ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();

#ifdef HAVE_TPETRA_TBB
    ScalarType localDP = task_scheduler_is_alive() ?
      threaded_vector_dot(length, scalarPointer(), x.scalarPointer()) :
      BLAS().DOT(length, scalarPointer(), ordinalOne, x.scalarPointer(), ordinalOne);
#else 
    // call BLAS routine to calculate local dot product
    ScalarType localDP = BLAS().DOT(length, scalarPointer(), ordinalOne, x.scalarPointer(), ordinalOne);
#endif

    // use Comm call to sum all local dot products
    ScalarType globalDP;
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_SUM,localDP,&globalDP);

    // update flops counter: 2n-1
    updateFlops(length + length - ordinalOne);

    return(globalDP);
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::absoluteValue(Vector<OrdinalType, ScalarType> const& x) {
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const length = getNumMyEntries();

    for(OrdinalType i = ordinalZero; i < length; i++)
      VectorData_->scalarArray_[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(x[i]);
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::reciprocal(Vector<OrdinalType, ScalarType> const& x) {
    if(! vectorSpace().isCompatible(x.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    ScalarType const scalarOne = Teuchos::ScalarTraits<ScalarType>::one();
    OrdinalType const length = getNumMyEntries();

    for(OrdinalType i = ordinalZero; i < length; i++)
      VectorData_->scalarArray_[i] = scalarOne / x[i];

    // update flops counter: n
    updateFlops(length);
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::scale(ScalarType scalarThis) {
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);

    // update flops counter: n
    updateFlops(length);
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x) {
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // this = x
    scalarArray() = x.scalarArray();
    // this = this * scalarX
    BLAS().SCAL(length, scalarX, scalarPointer(), ordinalOne);

    // update flops counter: n
    updateFlops(length);
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarThis) {
    if(! vectorSpace().isCompatible(x.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType const length = getNumMyEntries();

    //If we have Intel TBB thread support enabled, call the threaded_vector_update
    //otherwise, use a combination of BLAS.SCAL and BLAS.AXPY.

#     ifdef HAVE_TPETRA_TBB
    if (task_scheduler_is_alive()) {
      threaded_vector_update(length, scalarThis, scalarPointer(), scalarX, x.scalarPointer());
    }
    else {
      OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
      // calculate this *= scalarThis
      BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);

      // calculate this += scalarX * x
      BLAS().AXPY(length, scalarX, x.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);
    }
#     else
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    // calculate this *= scalarThis
    BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);

    // calculate this += scalarX * x
    BLAS().AXPY(length, scalarX, x.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);
#     endif

    // update flops counter: 3n
    updateFlops(length + length + length);
  }


  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarY, 
      Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis) {
    if(!vectorSpace().isCompatible(x.vectorSpace()) ||
        !vectorSpace().isCompatible(y.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // calculate this *= scalarThis
    BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);
    // calculate this += scalarX * x
    BLAS().AXPY(length, scalarX, x.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);
    // calculate this += scalarY * y
    BLAS().AXPY(length, scalarY, y.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);

    // update flops counter: 5n
    updateFlops(length + length + length + length + length);
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::norm1() const {
    // 1-norm = sum of abs. values of vector entries
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // compute local 1-norm
    ScalarType localNorm = BLAS().ASUM(length, scalarPointer(), ordinalOne);

    // call comm's sumAll method to compute global 1-norm
    ScalarType globalNorm;
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_SUM,localNorm,&globalNorm);

    // update flops counter: n-1
    updateFlops(length - ordinalOne);

    return(globalNorm);
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::norm2() const {
    // 2-norm = square root of the sum of the squares of the abs. values of vector entries
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // add up squares of entries
    ScalarType localSum = Teuchos::ScalarTraits<ScalarType>::zero();
    for(OrdinalType i = ordinalZero; i < length; i++)
      // the general (complex) case requires the
      // conjugate of the first term, which is not
      // needed for real numbers.
      //localSum += scalarArray()[i] * scalarArray()[i];
      localSum += Teuchos::ScalarTraits<ScalarType>::conjugate(scalarArray()[i]) * scalarArray()[i];

    // calculate global sum
    ScalarType globalSum;
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_SUM,localSum,&globalSum);

    // update flops counter: 2n
    updateFlops(length + length);

    // return square root of global sum
    return(Teuchos::ScalarTraits<ScalarType>::squareroot(globalSum));
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::normInf() const {
    // inf-norm = abs. value of the max value
    return(Teuchos::ScalarTraits<ScalarType>::magnitude(maxValue()));
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::normWeighted(Vector<OrdinalType, ScalarType> const& weights) const {
    if(!vectorSpace().isCompatible(weights.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // add up this[i] * weights[i]
    ScalarType localSum = Teuchos::ScalarTraits<ScalarType>::zero();
    for(OrdinalType i = ordinalZero; i < length; i++) {
      ScalarType temp = scalarArray()[i] * weights[i];
      localSum += temp * temp;
    }

    // get global sum
    ScalarType globalSum;
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_SUM,localSum,&globalSum);

    // divide by global length, and then take square root of that
    globalSum /= static_cast<ScalarType>(getNumGlobalEntries());

    // update flops counter: 3n
    updateFlops(3 * getNumGlobalEntries());

    return(Teuchos::ScalarTraits<ScalarType>::squareroot(globalSum));
  }
  
  // finish: move these somewhere?
  template<typename OrdinalType, typename ScalarType>
  struct less_mag : public binary_function<ScalarType, ScalarType, bool> {
    bool operator()(ScalarType x, ScalarType y) { 
      return Teuchos::ScalarTraits<ScalarType>::magnitude(x) < 
        Teuchos::ScalarTraits<ScalarType>::magnitude(y); 
    }
  };

  template<typename OrdinalType, typename ScalarType>
  struct greater_mag : public binary_function<ScalarType, ScalarType, bool> {
    bool operator()(ScalarType x, ScalarType y) { 
      return Teuchos::ScalarTraits<ScalarType>::magnitude(x) > 
        Teuchos::ScalarTraits<ScalarType>::magnitude(y); 
    }
  };


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::minValue() const {
    ScalarType localMin, globalMin;
    localMin = *(max_element(scalarArray().begin(), scalarArray().end(), less_mag<OrdinalType,ScalarType>())); // use STL max_element, takes constant time

    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_SUM,localMin,&globalMin);
    return globalMin;
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::maxValue() const {
    ScalarType localMax, globalMax;
    localMax = *(max_element(scalarArray().begin(), scalarArray().end(), greater_mag<OrdinalType,ScalarType>())); // use STL max_element, takes constant time

    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_MAX,localMax,&globalMax);
    return globalMax;
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::meanValue() const {
    ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
    // use STL accumulate, takes linear time
    ScalarType localTotal = accumulate(scalarArray().begin(), scalarArray().end(), scalarZero); 

    ScalarType globalTotal;
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    Teuchos::reduceAll(vectorSpace().comm(),Teuchos::REDUCE_SUM,localTotal,&globalTotal);

    // update flops counter: n
    updateFlops(getNumGlobalEntries());

    return(globalTotal / static_cast<ScalarType>(getNumGlobalEntries()));
  }



  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::elementwiseMultiply(ScalarType scalarXY, Vector<OrdinalType, ScalarType> const& x, 
      Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis) {
    if(!vectorSpace().isCompatible(x.vectorSpace()) ||
        !vectorSpace().isCompatible(y.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // calculate x@y into temp vector
    vector<ScalarType> xytemp(length);
    transform(x.scalarArray().begin(), x.scalarArray().end(), 
        y.scalarArray().begin(), xytemp.begin(), multiplies<ScalarType>());

    // calculate this *= scalarThis
    BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);

    // calculate this = scalarXY * temp + this
    BLAS().AXPY(length, scalarXY, &xytemp[ordinalZero], 
        ordinalOne, scalarPointer(), ordinalOne);

    // update flops counter: n
    updateFlops(length);
  }



  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::elementwiseReciprocalMultiply(ScalarType scalarXY, 
      Vector<OrdinalType, ScalarType> const& x, 
      Vector<OrdinalType, ScalarType> const& y, 
      ScalarType scalarThis) {
    if(!vectorSpace().isCompatible(x.vectorSpace()) ||
        !vectorSpace().isCompatible(y.vectorSpace()))
      throw Teuchos::Object::reportError("Vector sizes do not match.", 2);

    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const length = getNumMyEntries();

    // calculate y@x into temp vector
    vector<ScalarType> xytemp(length);
    transform(y.scalarArray().begin(), y.scalarArray().end(), 
        x.scalarArray().begin(), xytemp.begin(), divides<ScalarType>());

    // calculate this *= scalarThis
    BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);

    // calculate this += scalarXY * temp
    BLAS().AXPY(length, scalarXY, &xytemp[ordinalZero], ordinalOne, 
        scalarPointer(), ordinalOne);

    // update flops counter: 2n
    updateFlops(length + length);
  }



  template<typename OrdinalType, typename ScalarType>
  ScalarType Vector<OrdinalType,ScalarType>::getSeed() const {
    return(VectorData_->seed_);
  }


  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::setSeed(ScalarType seed) {
    VectorData_->seed_ = seed;
  }



  template<typename OrdinalType, typename ScalarType>
  ScalarType& Vector<OrdinalType,ScalarType>::operator[](OrdinalType index) {
    return(VectorData_->scalarArray_[index]);
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType const& Vector<OrdinalType,ScalarType>::operator[](OrdinalType index) const {
    return(VectorData_->scalarArray_[index]);
  }



  template<typename OrdinalType, typename ScalarType>
  OrdinalType Vector<OrdinalType,ScalarType>::getNumMyEntries() const {
    return(vectorSpace().getNumMyEntries());
  }


  template<typename OrdinalType, typename ScalarType>
  OrdinalType Vector<OrdinalType,ScalarType>::getNumGlobalEntries() const {
    return(vectorSpace().getNumGlobalEntries());
  }



  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::print(ostream& os) const {
    OrdinalType const myImageID = vectorSpace().comm().getRank();
    OrdinalType const numImages = vectorSpace().comm().getSize();
    OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();

    for (OrdinalType imageCtr = ordinalZero; imageCtr < numImages; imageCtr++) {
      if (myImageID == imageCtr) {
        if (myImageID == ordinalZero) {
          os << Teuchos::Object::label() << endl;
          os << "Number of Global Entries  = " << getNumGlobalEntries() << endl;
        }
        os <<   "ImageID = " << myImageID << endl;
        os <<           "Number of Local Entries   = " << getNumMyEntries() << endl;
        os <<           "Contents: ";
        printValues(os);
      }
      vectorSpace().comm().barrier();
    }
  }

  template<typename OrdinalType, typename ScalarType>
  void Vector<OrdinalType,ScalarType>::printValues(ostream& os) const {
    for(typename std::vector<ScalarType>::const_iterator i = scalarArray().begin(); 
        i != scalarArray().end(); i++)
      os << *i << " ";
    os << endl;
  }


  template<typename OrdinalType, typename ScalarType>
  VectorSpace<OrdinalType, ScalarType> const& Vector<OrdinalType,ScalarType>::vectorSpace() const {
    return(VectorData_->VectorSpace_);
  }


  template<typename OrdinalType, typename ScalarType>
  Vector<OrdinalType, ScalarType>& Vector<OrdinalType,ScalarType>::operator = (Vector<OrdinalType, ScalarType> const& Source) {
    VectorData_ = Source.VectorData_;
    return(*this);
  }



  template<typename OrdinalType, typename ScalarType>
  ScalarType* Vector<OrdinalType,ScalarType>::scalarPointer() {
    if(VectorData_->scalarArray_.empty())
      return(0);
    else
      return(&VectorData_->scalarArray_[Teuchos::OrdinalTraits<OrdinalType>::zero()]);
  }
  template<typename OrdinalType, typename ScalarType>
  ScalarType const* Vector<OrdinalType,ScalarType>::scalarPointer() const {
    if(VectorData_->scalarArray_.empty())
      return(0);
    else
      return(&VectorData_->scalarArray_[Teuchos::OrdinalTraits<OrdinalType>::zero()]);
  }

  template<typename OrdinalType, typename ScalarType>
  Teuchos::BLAS<OrdinalType, ScalarType> const& Vector<OrdinalType,ScalarType>::BLAS() const {
    return(VectorData_->BLAS_);
  }

  template<typename OrdinalType, typename ScalarType>
  std::vector<ScalarType>& Vector<OrdinalType,ScalarType>::scalarArray() {
    return(VectorData_->scalarArray_);
  }
  template<typename OrdinalType, typename ScalarType>
  std::vector<ScalarType>const & Vector<OrdinalType,ScalarType>::scalarArray() const{
    return(VectorData_->scalarArray_);
  }

  template<typename OrdinalType, typename ScalarType>
  bool Vector<OrdinalType,ScalarType>::checkSizes(DistObject<OrdinalType, ScalarType> const& sourceObj) {
    // first check that sourceObj is actually a Vector, and not some other kind of DistObject
    /*try {
      Vector<OrdinalType, ScalarType> const& sourceVector = dynamic_cast<Vector<OrdinalType, ScalarType> const&>(sourceObj);
      }
      catch(std::bad_cast bc) {
      return(false);
      }*/

    // ???

    return(true);
  }

  template<typename OrdinalType, typename ScalarType>
  int Vector<OrdinalType,ScalarType>::copyAndPermute(DistObject<OrdinalType, ScalarType> const& sourceObj,
      OrdinalType const numSameIDs,
      OrdinalType const numPermuteIDs,
      std::vector<OrdinalType> const& permuteToLIDs,
      std::vector<OrdinalType> const& permuteFromLIDs) {
    // cast sourceObj to a Tpetra::Vector so we can actually do something with it
    Vector<OrdinalType, ScalarType> const& sourceVector = dynamic_cast<Vector<OrdinalType, ScalarType> const&>(sourceObj);

    std::vector<ScalarType> const& sourceArray = sourceVector.scalarArray();
    std::vector<ScalarType>& destArray = scalarArray();

    // the first numImportIDs GIDs are the same between source and target,
    // we can just copy them
    for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < numSameIDs; i++)
      destArray[i] = sourceArray[i];

    // next, do permutations
    for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < numPermuteIDs; i++)
      destArray[permuteToLIDs[i]] = sourceArray[permuteFromLIDs[i]];

    return(0);
  }

  template<typename OrdinalType, typename ScalarType>
  int Vector<OrdinalType,ScalarType>::packAndPrepare(DistObject<OrdinalType, ScalarType> const& sourceObj,
      OrdinalType const numExportIDs,
      std::vector<OrdinalType> const& exportLIDs,
      std::vector<ScalarType>& exports,
      OrdinalType& packetSize,
      Distributor<OrdinalType> const& distor) {
    // cast sourceObj to a Tpetra::Vector so we can actually do something with it
    Vector<OrdinalType, ScalarType> const& sourceVector = dynamic_cast<Vector<OrdinalType, ScalarType> const&>(sourceObj);

    // For a vector, we only need to send the value
    exports.clear();
    for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < numExportIDs; i++)
      exports.push_back(sourceVector[exportLIDs[i]]);

    // packetSize = 1
    packetSize = Teuchos::OrdinalTraits<OrdinalType>::one();

    return(0);
  }

  template<typename OrdinalType, typename ScalarType>
  int Vector<OrdinalType,ScalarType>::unpackAndCombine(OrdinalType const numImportIDs,
      std::vector<OrdinalType> const& importLIDs,
      std::vector<ScalarType> const& imports,
      Distributor<OrdinalType> const& distor,
      CombineMode const CM) {
    if(CM == Insert || CM == Replace) {
      // copy values from scalarExports
      for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < numImportIDs; i++)
        scalarArray().at(importLIDs[i]) = imports[i];
    }
    else if(CM == Add) {
      // sum values from scalarExports
      submitEntries(numImportIDs, &importLIDs.front(), &imports.front());
    }
    else
      throw Teuchos::Object::reportError("Unknown CombineMode", -99);

    return(0);
  }

} // namespace Tpetra

#endif // TPETRA_VECTOR_HPP
