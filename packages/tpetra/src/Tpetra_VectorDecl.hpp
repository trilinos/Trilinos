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

#ifndef TPETRA_VECTOR_DECL_HPP
#define TPETRA_VECTOR_DECL_HPP

#include "Tpetra_MultiVectorDecl.hpp"

namespace Tpetra {

  //! Tpetra::Vector: A class for constructing and using distributed vectors.

  template<typename Ordinal, typename Scalar>
  class Vector : public MultiVector<Ordinal,Scalar> {

  public:
  
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    Vector(const Map<Ordinal> &map, bool zeroOut=true);

    //! Vector copy constructor.
    Vector(const Vector<Ordinal,Scalar> &source);

    //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
    Vector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &values);

    //! Destructor.  
    virtual ~Vector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified location with specified value.
    void replaceGlobalValue(Ordinal globalRow, const Scalar &value);

    //! Adds specified value to existing value at the specified location.
    void sumIntoGlobalValue(Ordinal globalRow, const Scalar &value);

    //! Replace current value at the specified location with specified values.
    void replaceMyValue(Ordinal myRow, const Scalar &value);

    //! Adds specified value to existing value at the specified location.
    void sumIntoMyValue(Ordinal myRow, const Scalar &value);

    //! Initialize all values in a multi-vector with specified value.
    void putScalar(const Scalar &value);

    //! Set multi-vector values to random numbers.
    void random();

    //! Replace the underlying Map with a compatible one.
    void replaceMap(const Map<Ordinal> &map);

    //! Instruct a local (non-distributed) MultiVector to sum values across all nodes.
    void reduce();


    //@}

    //! @name Extraction methods
    //@{

    using MultiVector<Ordinal,Scalar>::extractCopy; // overloading, not hiding
    //! Return multi-vector values in user-provided two-dimensional array.
    void extractCopy(Teuchos::ArrayView<Scalar> A) const;

    using MultiVector<Ordinal,Scalar>::extractView; // overloading, not hiding
    //! Return non-const non-persisting view of values in a one-dimensional array.
    void extractView(Teuchos::ArrayView<Scalar> &A);

    using MultiVector<Ordinal,Scalar>::extractConstView; // overloading, not hiding
    //! Return const non-persisting view of values in a one-dimensional array.
    void extractConstView(Teuchos::ArrayView<const Scalar> &A) const;

    //@}

    //! @name Mathematical methods
    //@{ 

    // FINISH: expand documentation of these functions

    using MultiVector<Ordinal,Scalar>::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    Scalar dot(const Vector<Ordinal,Scalar> &a) const;

    using MultiVector<Ordinal,Scalar>::abs; // overloading, not hiding
    //! Puts element-wise absolute values of input Vector in this Vector: a = abs(this)
    void abs(const Vector<Ordinal,Scalar> &a);

    using MultiVector<Ordinal,Scalar>::reciprocal; // overloading, not hiding
    //! Puts element-wise reciprocal values of input Vector in this Vector: this(i) = 1/a(i).
    void reciprocal(const Vector<Ordinal,Scalar> &a);

    using MultiVector<Ordinal,Scalar>::scale; // overloading, not hiding
    //! Replace this Vector with scaled Vector a, this = alpha*a.
    void scale(const Scalar &alpha, const Vector<Ordinal,Scalar> &a);

    using MultiVector<Ordinal,Scalar>::update; // overloading, not hiding
    //! Update this Vector with scaled Vector a, this = beta*this + alpha*a.
    void update(const Scalar &alpha, const Vector<Ordinal,Scalar> &a, const Scalar &beta);

    //! Update this Vector with scaled Vectors a and b, this = gamma*this + alpha*a + beta*b.
    void update(const Scalar &alpha, const Vector<Ordinal,Scalar> &a, const Scalar &beta, const Vector<Ordinal,Scalar> &b, const Scalar &gamma);

    using MultiVector<Ordinal,Scalar>::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const;

    using MultiVector<Ordinal,Scalar>::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const;

    using MultiVector<Ordinal,Scalar>::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const;

    using MultiVector<Ordinal,Scalar>::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const MultiVector<Ordinal,Scalar> &weights) const;

    using MultiVector<Ordinal,Scalar>::minValue; // overloading, not hiding
    //! Compute minimum value of this Vector.
    Scalar minValue() const;

    using MultiVector<Ordinal,Scalar>::maxValue; // overloading, not hiding
    //! Compute maximum value of this Vector.
    Scalar maxValue() const;

    using MultiVector<Ordinal,Scalar>::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    Scalar meanValue() const;

    //@{ 

    //! @name Element access methods
    //@{ 

    //! [] operator, nonconst version
    Scalar& operator[](Ordinal index);

    //! [] operator, const version
    const Scalar & operator[](Ordinal index) const;

    //@}


    //! @name Attribute access methods
    //@{ 

    //! Returns number of vector entries owned by this image.
    Ordinal getNumMyEntries() const;

    //! Returns number of vector entries across all images.
    Ordinal getNumGlobalEntries() const;

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method, used by overloaded << operator.
    void print(std::ostream &os) const;

    void printValues(std::ostream &os) const;

    //@}

    //! @name Misc. 
    //@{ 

    //! Returns a const reference to the VectorSpace this Vector belongs to.
    const Map<Ordinal> & getMap() const;

    //! Assignment Operator
    Vector<Ordinal,Scalar> & operator=(const Vector<Ordinal,Scalar> &source);

    //@}

    //! @name Expert/Developer Use Only.
    //@{ 

    // Returns pointer to Scalar array inside of scalarArray
    Teuchos::ArrayView<Scalar> scalarPointer();

    Teuchos::ArrayView<const Scalar> scalarPointer() const;

    //@}

  private:

    Teuchos::RCP<VectorData<Ordinal,Scalar> > VectorData_;

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Ordinal,Scalar> & sourceObj);

    void copyAndPermute(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numSameIDs,
               Ordinal numPermuteIDs,
               const Teuchos::ArrayView<const Ordinal> & permuteToLIDs,
               const Teuchos::ArrayView<const Ordinal> & permuteFromLIDs);

    void packAndPrepare(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numExportIDs,
               const Teuchos::ArrayView<const Ordinal> & exportLIDs,
               const Teuchos::ArrayView<Scalar> & exports,
               Ordinal &packetSize,
               Distributor<Ordinal> &distor);
  
    void unpackAndCombine(Ordinal numImportIDs,
               const Teuchos::ArrayView<const Ordinal> & importLIDs,
               const Teuchos::ArrayView<const Scalar> & imports,
               Distributor<Ordinal> &distor,
               CombineMode CM);

  }; // class Vector

} // namespace Tpetra

#endif // TPETRA_VECTOR_DECL_HPP
