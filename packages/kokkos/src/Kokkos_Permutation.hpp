//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2003) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef KOKKOS_PERMUTATION_H
#define KOKKOS_PERMUTATION_H

#include "Kokkos_ConfigDefs.hpp"

namespace Kokkos {

//! Kokkos::Permutation: Kokkos permutation base class.

/*! The Kokkos::Permutation specifies the interface that any permutation class which is intended
    for use with the Kokkos::Vector and Kokkos::MultiVector classes must implement.

    Formally a permutation is a rearrangement of the rows of the identity matrix.  Permutations can
    applied to the left of a matrix or vector, rearranging the rows or entries, respectively.  A
    permutation can also be applied on the right of a matrix (or a vector tranpose), thereby rearranging
    the columns of the matrix.  Permutations also have the property that the inverse and 
    the transpose are identical.

    This class is most important as a way for a Kokkos::SparseOperation object to define a row, column or
    two-sided permutation of the matrix that presumably improves performance of the sparse operation.
    The Kokkos::SparseOperation object is responsible for producing the left and right permutations objects
    on demand.  In turn, these permutation objects can be used to pre-permute vectors and multivectors so
    that their ordering is compatible with what is required for the sparse operation.

*/    

  template<typename OrdinalType, typename ScalarType>
  class Permutation {
  public:

    //@{ \name Constructors/Destructor.

    //! Default constructor
    Permutation(void):
      isIdentity_(true),
      dataInitialized_(false),
      length_(0),
      indices_(0) {};
  
    //! Non-identity permutation constructor.
    Permutation (const OrdinalType * indices, OrdinalType length):
      isIdentity_(indices==0),
      dataInitialized_(indices!=0),
      length_(length),
      indices_(indices) {};

    //! Copy constructor.
    Permutation(const Permutation& source):
      isIdentity_(source.isIdentity),
      dataInitialized_(source.dataInitialized_),
      length_(source.length_),
      indices_(source.indices_) {
    };

    //! Permutation Destructor
    virtual ~Permutation(){};
    //@}


    //@{ \name Computational methods.
	
    //! Returns the result of a Kokkos::Permutation applied to a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to permute.
      \param y (Out) A Kokkos::Vector containing results.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y) const;

    //! Returns the result of a Kokkos::Permutation applied to multiple vectors in x, results in y.
    /*! 
      \param x (In) A Kokkos::MultiVector to permute.
      \param y (Out) A Kokkos::MultiVector containing results.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y) const;

    //! Returns the result of a Kokkos::Permutation inverse applied to a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to permute.
      \param y (Out) A Kokkos::Vector containing results.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int applyInverse(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y) const;

    //! Returns the result of a Kokkos::Permutation inverse applied to multiple vectors in x, results in y.
    /*! 
      \param x (In) A Kokkos::MultiVector to permute.
      \param y (Out) A Kokkos::MultiVector containing results.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int applyInverse(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y) const;

    //@}
	
    //@{ \name Permutation Attribute access methods.
	
    //! Length (dimension) of permutation operator.
    virtual OrdinalType getLength() const {return(length_);};
	
    //! Returns true if the permutation is the identity, otherwise returns false.
    virtual bool getIsIdentity() const {return(isIdentity_);};
	
    //@}

    bool isIdentity_;
    bool dataInitialized_;
    OrdinalType length_;

    OrdinalType * indices_;
  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int Permutation<OrdinalType, ScalarType>::apply(const Vector<OrdinalType, ScalarType>& x, 
						    Vector<OrdinalType, ScalarType> & y) const {

    if (x.getLength()!=y.getLength()) return(-1); // Incompatible x and y

    ScalarType * xp = x.getValues();
    ScalarType * yp = y.getValues();
    OrdinalType length = x.getLength();

    for(OrdinalType i = 0; i < length; i++)
      yp[indices_[i]] = xp[i];
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int Permutation<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
						  MultiVector<OrdinalType, ScalarType> & y) const {
    if (x.getNumRows()!=y.getNumRows()) return(-1); // Incompatible x and y length
    OrdinalType numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-2); // Not the same number of vectors in x and y

    ScalarType ** xpp = x.getValues();
    ScalarType ** ypp = y.getValues();
    OrdinalType length = x.getNumRows();

    for (OrdinalType k=0; k<numVectors; k++) {
      ScalarType * xp = xpp[k];
      ScalarType * yp = ypp[k];
      for(OrdinalType i = 0; i < length; i++)
	yp[indices_[i]] = xp[i];
    }
    return(0);
  }


  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int Permutation<OrdinalType, ScalarType>::applyInverse(const Vector<OrdinalType, ScalarType>& x, 
							 Vector<OrdinalType, ScalarType> & y) const {

    if (x.getLength()!=y.getLength()) return(-1); // Incompatible x and y

    ScalarType * xp = x.getValues();
    ScalarType * yp = y.getValues();
    OrdinalType length = x.getLength();

    for(OrdinalType i = 0; i < length; i++)
      yp[i] = xp[indices_[i]];
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int Permutation<OrdinalType, ScalarType>::applyInverse(const MultiVector<OrdinalType, ScalarType>& x, 
							 MultiVector<OrdinalType, ScalarType> & y) const {
    if (x.getNumRows()!=y.getNumRows()) return(-1); // Incompatible x and y length
    OrdinalType numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-2); // Not the same number of vectors in x and y

    ScalarType ** xpp = x.getValues();
    ScalarType ** ypp = y.getValues();
    OrdinalType length = x.getNumRows();

    for (OrdinalType k=0; k<numVectors; k++) {
      ScalarType * xp = xpp[k];
      ScalarType * yp = ypp[k];
      for(OrdinalType i = 0; i < length; i++)
	yp[i] = xp[indices_[i]];
    }
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_PERMUTATION_H */
