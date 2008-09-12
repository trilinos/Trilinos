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
#include "Tpetra_VectorData.hpp"

namespace Tpetra {

  template<typename Ordinal,typename Scalar>
  Vector<Ordinal,Scalar>::Vector(const Map<Ordinal> &map) 
    : DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::Vector")
    , VectorData_()
  {
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
    VectorData_ = Teuchos::rcp(new VectorData<Ordinal,Scalar>(map, ZERO));
  }

  template<typename Ordinal,typename Scalar>
  Vector<Ordinal,Scalar>::Vector(const Teuchos::ArrayView<const Scalar> &values, const Map<Ordinal> &map)
  : DistObject<Ordinal,Scalar>(map,map.getPlatform()->createComm(), "Tpetra::Vector")
  , VectorData_()
  {
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
    const Ordinal OZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal length = map.getNumMyEntries();

    VectorData_ = Teuchos::rcp(new VectorData<Ordinal,Scalar>(map, ZERO));

    TEST_FOR_EXCEPTION(values.size() != length, std::invalid_argument,
      "Tpetra::Vector::constructor(values,map): number of values does not match size of map.");
    std::copy(values.begin(),values.end(), VectorData_->values_().begin());
  }

  template<typename Ordinal, typename Scalar>
  Vector<Ordinal,Scalar>::Vector(const Vector<Ordinal,Scalar> &source)
  : DistObject<Ordinal,Scalar>(source)
  , VectorData_(source.VectorData_)
  {}

  template<typename Ordinal, typename Scalar>
  Vector<Ordinal,Scalar>::~Vector() {}

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::submitEntries(const Teuchos::ArrayView<const Ordinal> &indices,
                                             const Teuchos::ArrayView<const Scalar>  &values)
  {
    TEST_FOR_EXCEPTION(indices.size() != values.size(), std::invalid_argument,
        "Tpetra::Vector::submitEntries(indices,values): finish.");
    const Ordinal ordinalZero = Teuchos::OrdinalTraits<Ordinal>::zero();
    typename Teuchos::ArrayView<const Ordinal>::const_iterator ind = indices.begin();
    typename Teuchos::ArrayView<const Scalar>::const_iterator  val = values.begin();
    for (; ind != indices.end(); ++ind, ++val) {
      // if TEUCHOS_DEBUG, this is bounds checked
      VectorData_->values_[*ind] += *val;
    }
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::setAllToScalar(const Scalar &value) {
    TEST_FOR_EXCEPT(true); // finish
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::setAllToRandom() {
    TEST_FOR_EXCEPT(true); // finish
  }

  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::dotProduct(const Vector<Ordinal,Scalar> &x) const 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::dotProduct(): Vectors must have compatible Maps.");
    Teuchos::BLAS<Ordinal,Scalar> blas;
    // call BLAS routine to calculate local dot product
    // use Comm call to sum all local dot products
    // update flops counter
    TEST_FOR_EXCEPT(true);
    return Teuchos::ScalarTraits<Scalar>::zero();
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::absoluteValue(const Vector<Ordinal,Scalar> &x) 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::absoluteValue(): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::reciprocal(const Vector<Ordinal,Scalar> &x) 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::reciprocal(): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::scale(const Scalar &scalarThis) 
  {
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::scale(const Scalar &scalarX, const Vector<Ordinal,Scalar> &x) 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::scale(): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::update(const Scalar &scalarX, const Vector<Ordinal,Scalar> &x, 
                                      const Scalar &scalarThis) 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::update(): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::update(const Scalar &scalarX, const Vector<Ordinal,Scalar> &x, 
                                      const Scalar &scalarY, const Vector<Ordinal,Scalar> &y, 
                                      const Scalar &scalarThis) 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::update(x,y): Vectors must have compatible Maps.");
    TEST_FOR_EXCEPTION( !getMap().isCompatible(y.getMap()), 
        std::runtime_error, "Tpetra::Vector::update(x,y): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::norm1() const {
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::norm2() const {
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::normInf() const {
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::normWeighted(const Vector<Ordinal,Scalar> &weights) const 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(weights.getMap()), 
        std::runtime_error, "Tpetra::Vector::normWeighted(): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  // finish: move these somewhere?
  /*
  template<typename Ordinal, typename Scalar>
  struct vector_less_mag : public binary_function<Scalar,Scalar,bool> {
    bool operator()(Scalar x, Scalar y) { 
      return Teuchos::ScalarTraits<Scalar>::magnitude(x) < 
        Teuchos::ScalarTraits<Scalar>::magnitude(y); 
    }
  };

  template<typename Ordinal, typename Scalar>
  struct vector_greater_mag : public binary_function<Scalar, Scalar, bool> {
    bool operator()(Scalar x, Scalar y) { 
      return Teuchos::ScalarTraits<Scalar>::magnitude(x) > 
        Teuchos::ScalarTraits<Scalar>::magnitude(y); 
    }
  };
  */


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::minValue() const {
    TEST_FOR_EXCEPT(true);
    return Teuchos::ScalarTraits<Scalar>::zero();
  }


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::maxValue() const {
    TEST_FOR_EXCEPT(true);
    return Teuchos::ScalarTraits<Scalar>::zero();
  }


  template<typename Ordinal, typename Scalar>
  Scalar Vector<Ordinal,Scalar>::meanValue() const {
    TEST_FOR_EXCEPT(true);
    // update flops
    return Teuchos::ScalarTraits<Scalar>::zero();
  }


  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::elementwiseMultiply(const Scalar &scalarXY, const Vector<Ordinal,Scalar> &x, const Vector<Ordinal,Scalar> &y, 
                                                   const Scalar &scalarThis)
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::elementwiseMultiply(x,y): Vectors must have compatible Maps.");
    TEST_FOR_EXCEPTION( !getMap().isCompatible(y.getMap()), 
        std::runtime_error, "Tpetra::Vector::elementwiseMultiply(x,y): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::elementwiseReciprocalMultiply(
      Scalar scalarXY, const Vector<Ordinal, Scalar> &x, const Vector<Ordinal, Scalar> &y, 
      const Scalar &scalarThis) 
  {
    TEST_FOR_EXCEPTION( !getMap().isCompatible(x.getMap()), 
        std::runtime_error, "Tpetra::Vector::elementwiseReciprocalMultiply(x,y): Vectors must have compatible Maps.");
    TEST_FOR_EXCEPTION( !getMap().isCompatible(y.getMap()), 
        std::runtime_error, "Tpetra::Vector::elementwiseReciprocalMultiply(x,y): Vectors must have compatible Maps.");
    // do work
    // update flops counter
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  const Scalar & Vector<Ordinal,Scalar>::getSeed() const {
    return(VectorData_->seed_);
  }


  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::setSeed(const Scalar &seed) {
    VectorData_->seed_ = seed;
  }



  template<typename Ordinal, typename Scalar>
  Scalar& Vector<Ordinal,Scalar>::operator[](Ordinal index) {
    return VectorData_->scalarArray_[index];
  }


  template<typename Ordinal, typename Scalar>
  const Scalar & Vector<Ordinal,Scalar>::operator[](Ordinal index) const {
    return VectorData_->scalarArray_[index];
  }


  template<typename Ordinal, typename Scalar>
  Ordinal Vector<Ordinal,Scalar>::getNumMyEntries() const {
    return VectorData_->map_.getNumMyEntries();
  }


  template<typename Ordinal, typename Scalar>
  Ordinal Vector<Ordinal,Scalar>::getNumGlobalEntries() const {
    return VectorData_->map_.getNumGlobalEntries();
  }


  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::print(std::ostream& os) const {
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::printValues(std::ostream& os) const {
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  const Map<Ordinal> & Vector<Ordinal,Scalar>::getMap() const {
    return VectorData_->map_;
  }


  template<typename Ordinal, typename Scalar>
  Vector<Ordinal,Scalar> & Vector<Ordinal,Scalar>::operator=(const Vector<Ordinal,Scalar> &source) 
  {
    VectorData_ = source.VectorData_;
    return *this;
  }


  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayView<Scalar> Vector<Ordinal,Scalar>::scalarPointer() 
  {
    return VectorData_->values_;
  }

  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayView<const Scalar> Vector<Ordinal,Scalar>::scalarPointer() const
  {
    return VectorData_->values_;
  }

  template<typename Ordinal, typename Scalar>
  bool Vector<Ordinal,Scalar>::checkSizes(const DistObject<Ordinal,Scalar> &sourceObj) {
    // first check that sourceObj is actually a Vector, and not some other kind of DistObject
    /*try {
      Vector<Ordinal, Scalar> const& sourceVector = dynamic_cast<Vector<Ordinal, Scalar> const&>(sourceObj);
      }
      catch(std::bad_cast bc) {
      return(false);
      }*/

    // ???
    TEST_FOR_EXCEPT(true);
    return true;
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::copyAndPermute(
      const DistObject<Ordinal,Scalar> &sourceObj,
      Ordinal numSameIDs, Ordinal numPermuteIDs,
      const Teuchos::ArrayView<const Ordinal> &permuteToLIDs, 
      const Teuchos::ArrayView<const Ordinal> &permuteFromLIDs) 
  {
    // cast sourceObj to a Tpetra::Vector so we can actually do something with it
    // const Vector<Ordinal,Scalar> &sourceVector = dynamic_cast<const Vector<Ordinal,Scalar> &>(sourceObj);

    /*
    const std::vector<Scalar> &sourceArray = sourceVector.scalarArray();
    std::vector<Scalar> &destArray = scalarArray();

    // the first numImportIDs GIDs are the same between source and target,
    // we can just copy them
    for(Ordinal i = Teuchos::OrdinalTraits<Ordinal>::zero(); i < numSameIDs; i++)
      destArray[i] = sourceArray[i];

    // next, do permutations
    for(Ordinal i = Teuchos::OrdinalTraits<Ordinal>::zero(); i < numPermuteIDs; i++)
      destArray[permuteToLIDs[i]] = sourceArray[permuteFromLIDs[i]];

     */
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::packAndPrepare(
      const DistObject<Ordinal,Scalar> &sourceObj,
      Ordinal numExportIDs,
      const Teuchos::ArrayView<const Ordinal> &exportLIDs, 
      const Teuchos::ArrayView<Scalar> &exports,
      Ordinal &packetSize,
      Distributor<Ordinal> &distor) 
  {
    /*
    // cast sourceObj to a Tpetra::Vector so we can actually do something with it
    Vector<Ordinal, Scalar> const& sourceVector = dynamic_cast<Vector<Ordinal, Scalar> const&>(sourceObj);

    // For a vector, we only need to send the value
    exports.clear();
    for(Ordinal i = Teuchos::OrdinalTraits<Ordinal>::zero(); i < numExportIDs; i++)
      exports.push_back(sourceVector[exportLIDs[i]]);

    // packetSize = 1
    packetSize = Teuchos::OrdinalTraits<Ordinal>::one();
    */
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void Vector<Ordinal,Scalar>::unpackAndCombine(
      Ordinal numImportIDs,
      const Teuchos::ArrayView<const Ordinal> &importLIDs,
      const Teuchos::ArrayView<const Scalar> &imports,
      Distributor<Ordinal> &distor,
      CombineMode CM) 
  {
    /*
    if(CM == Insert || CM == Replace) {
      // copy values from scalarExports
      for(Ordinal i = Teuchos::OrdinalTraits<Ordinal>::zero(); i < numImportIDs; i++)
        scalarArray().at(importLIDs[i]) = imports[i];
    }
    else if(CM == Add) {
      // sum values from scalarExports
      submitEntries(numImportIDs, &importLIDs.front(), &imports.front());
    }
    else
      throw Teuchos::Object::reportError("Unknown CombineMode", -99);
      */

    TEST_FOR_EXCEPT(true);
  }

} // namespace Tpetra

#endif // TPETRA_VECTOR_HPP
