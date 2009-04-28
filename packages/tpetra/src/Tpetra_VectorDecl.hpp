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

#ifndef TPETRA_VECTOR_DECL_HPP
#define TPETRA_VECTOR_DECL_HPP

#include "Tpetra_MultiVectorDecl.hpp"

namespace Tpetra {

  //! Tpetra::Vector: A class for constructing and using distributed vectors.

  template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal>
  class Vector : public MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> {

    // need this so that MultiVector::operator() can call Vector's private constructor
    friend class MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>;

  public:

    /** \name Convenient typedefs */
    //@{ 

    /*! Non-const pointer-like typedef. 
        In a debug build (<tt>--enable-teuchos-debug</tt>), this is an 
        ArrayRCP<Scalar>. In a non-debug build, it is a <tt>Scalar *</tt>. In either case, 
        the syntax is the same: pointer arithmetic and indexing are supported. */
    typedef typename Teuchos::ArrayView<Scalar>::iterator pointer;

    /*! Const pointer-like typedef. 
        In a debug build (<tt>--enable-teuchos-debug</tt>), this is an 
        ArrayRCP<const Scalar>. In a non-debug build, it is a <tt>const Scalar *</tt>. In either case, 
        the syntax is the same: pointer arithmetic and indexing are supported. */
    typedef typename Teuchos::ArrayView<const Scalar>::iterator const_pointer;

    //@}
  
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, bool zeroOut=true);

    //! Vector copy constructor.
    Vector(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &source);

    //! \brief Set multi-vector values from an array using a C pointer.
    /*! \c CopyView indicates whether the data will be copied from the input array or if the Vector object will encapsulate the data throughout its existence.
        \c OwnsMem indicates whether the Vector object owns the memory and is therefore responsible for deleting it. 
     */
    Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos::DataAccess CopyView, Scalar *ArrayOfPtrs, bool OwnsMem = false);

    //! \brief Set multi-vector values from an array using Teuchos memory management classes. (copy)
    Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &A);

    //! \brief Set multi-vector values from an using Teuchos memory management classes. (view)
    Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayRCP<Scalar> &A);

    //! Destructor.  
    virtual ~Vector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified location with specified value.
    void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

    //! Adds specified value to existing value at the specified location.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

    //! Replace current value at the specified location with specified values.
    void replaceMyValue(LocalOrdinal myRow, const Scalar &value);

    //! Adds specified value to existing value at the specified location.
    void sumIntoMyValue(LocalOrdinal myRow, const Scalar &value);

    //@}

    //! @name Extraction methods
    //@{

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::extractCopy1D; // overloading, not hiding
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void extractCopy1D(Teuchos::ArrayView<Scalar> A) const;
    //! Return multi-vector values in user-provided two-dimensional array (using a C pointer).
    void extractCopy1D(Scalar *A) const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::extractView1D; // overloading, not hiding
    //! Return non-const non-persisting view of values in a one-dimensional array (using C pointers).
    inline void extractView1D(Scalar * &A);
    //! Return non-const non-persisting view of values in a one-dimensional array (using Teuchos memory management classes).
    inline void extractView1D(Teuchos::ArrayView<Scalar> &A);

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::extractConstView1D; // overloading, not hiding
    //! Return const non-persisting view of values in a one-dimensional array (using C pointers).
    inline void extractConstView1D(const Scalar * &A) const;
    //! Return const non-persisting view of values in a one-dimensional array (using Teuchos memory management classes).
    inline void extractConstView1D(Teuchos::ArrayView<const Scalar> &A) const;

    //@}

    //! @name Mathematical methods
    //@{ 

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    Scalar dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &a) const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &weights) const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    Scalar meanValue() const;

    //@} 

    //! @name Entry access methods
    //@{ 

    //! [] operator, nonconst version
    Scalar& operator[](Teuchos_Ordinal index);

    //! [] operator, const version
    const Scalar & operator[](Teuchos_Ordinal index) const;

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method, used by overloaded << operator.
    void print(std::ostream &os) const;

    void printValues(std::ostream &os) const;

    //@}

    protected:

    // Advanced Vector constuctor for creating views.
    Vector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::RCP<MultiVectorData<Scalar> > &mvdata);


  }; // class Vector

} // namespace Tpetra

#endif // TPETRA_VECTOR_DECL_HPP
