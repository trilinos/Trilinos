//@HEADER
// ************************************************************************
// 
//                Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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

#ifndef KOKKOS_HBMATRIX_H
#define KOKKOS_HBMATRIX_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

//! Kokkos::HbMatrix: Kokkos compressed index sparse matrix base class.

/*! The Kokkos::HbMatrix implements the Kokkos::CisMatrix interface 
    using either a Harwell-Boeing matrix or generalized form of one.

*/    

  template<class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class HbMatrix {
  public:

    //! @name Constructors/Destructor

    //@{

    //! Default HbMatrix constuctor.
    HbMatrix(Node &node = DefaultNode::getDefaultNode());

    //! Copy constructor.
    HbMatrix(const HbMatrix& source);

    //! HbMatrix Destructor
    virtual ~HbMatrix();
    //@}

    //! @name Harwell-Boeing Format Initialization Methods

    //@{

    //! Initialize structure of matrix
    /*!
      This interface supports matrices that are stored in the Harwell-Boeing format.
      \param numRows (In)  Row dimension.
      \param pntr (In)  Array of offsets into indx.  indx[pntr[i]] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \param indx (In)  Packed array of indices.  indx[pntr[i]] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
      */
    int initializeStructure(Ordinal numRows, Ordinal * pntr, Ordinal * indx);

    //! Initialize structure of matrix
    /*!
      This interface supports matrices that are stored in the Harwell-Boeing format.
      \param values (In)  Packed array of matrix values. values[pntr[i]] contains the first entry of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
      */
    int initializeValues(Scalar * values);

    //@}

    //! @name Matrix entry access methods

    //@{

    //! Returns number of entries in ith row/column, and pointer to an array of these indices.
    /*! Extract the number of entries and a pointer to the indices in the ith row/column of the matrix.  Note that
        the indices are not copied by this method.  Memory allocation is handled by the matrix object itself.

        \param i (In) The row (if isRowOriented() is true) or column that should be returned.
        \param numRowEntries (Out) The number of entries in the ith row/column.
        \param indices (Out) A pointer to the list of indices in the ith row/column.

        \return Integer error code, set to 0 if successful.
    */
    int getIndices(Ordinal i, Ordinal & numRowEntries, Ordinal *& indices) const;

    //! Returns a pointer to an array of values for the ith row/column.
    /*! Extract the values in the ith row/column of the matrix.  Note that
        the values are not copied by this method.  Memory allocation is handled by the matrix object itself.

        \param i (In) The row (if isRowOriented() is true) or column that should be returned.
        \param numEntries (Out) The number of entries in the ith row/column.
        \param indices (Out) A pointer to the list of indices in the ith row/column.

        \return Integer error code, set to 0 if successful.
    */
    int getValues(Ordinal i, Scalar *& values) const;


    //@}

    //! @name Matrix Attribute set methods

    //@{

    //! Set whether or not the compressed index matrix has no entries below the diagonal, assumed false.
    virtual int setIsUpperTriangular(bool tf) {isUpperTriangular_=tf; return(0);};

    //! Set whether or not the compressed index matrix has no entries above the diagonal, assumed false.
    virtual int setIsLowerTriangular(bool tf) {isLowerTriangular_=tf; return(0);};

    //! Set whether or not the compressed index matrix has explicit diagonal entries.
    virtual int setHasDiagonalEntries(bool tf) {hasDiagonalEntries_=tf; return(0);};

    //! Set whether or not the compressed index matrix should be treated as unit diagonal, assumed false
    /*! \warning This method will not accept a "true" argument unless setHasDiagonalEntries() has be 
        called with a "false" argument first.
    */
    virtual int setHasImplicitUnitDiagonal(bool tf) {
      if (tf && hasDiagonalEntries_) return(-1); // Cannot set unit diagonal unless there are no explicit diagonal entries
      hasImplicitUnitDiagonal_=tf; 
      return(0);};

    //@}

    //! @name Validity tests

    //@{

    //! Check if the matrix structure is valid for user-assertion of Upper/Lower Triangular and implicit unit diagonal.
    /*!
      \return Integer error code, equals zero if structure is valid, negative value means potential fatal assertion.
      Error codes: <ul>
      <li> -1 - User asserts upper triangular, but it's not.
      <li> 1 - User does not assert upper triangular, but it is.
      <li> -2 - User asserts lower triangular, but it's not.
      <li> 2 - User does not assert lower triangular, but it is.
      <li> -3 - User asserts implicit unit diagonal, but values are given.
      </ul>
    */
    virtual int checkStructure() const;

    //@}

    //! @name Matrix Attribute access methods

    //@{

    //! Returns true if the compressed index matrix has no entries below the diagonal.
    virtual bool getIsUpperTriangular() const {return(isUpperTriangular_);};

    //! Returns true if the compressed index matrix has no entries above the diagonal.
    virtual bool getIsLowerTriangular() const {return(isLowerTriangular_);};

    //! Returns true if the compressed index matrix has no diagonal entries, but should be treated as unit diagonal.
    virtual bool getHasImplicitUnitDiagonal() const {return(hasImplicitUnitDiagonal_);};

    //! Number of rows
    Ordinal getNumRows() const {return(numRows_);};

    //! Number of matrix entries
    Ordinal getNumEntries() const {return(numEntries_);};

    //@}

  protected:

    Ordinal numRows_;
    Ordinal numEntries_;

    Scalar * values_;
    Ordinal * indices_;

    Ordinal * profile_;
    mutable bool isUpperTriangular_;
    mutable bool isLowerTriangular_;
    bool hasImplicitUnitDiagonal_;
    mutable bool hasDiagonalEntries_;
    Node &node_;
  };

} // namespace Kokkos

using namespace Kokkos;
//==============================================================================
template<typename Ordinal, typename Scalar>
 HbMatrix<Ordinal, Scalar>::HbMatrix() 
  : numRows_(0),
    numEntries_(0),
    values_(0),
    allValues_(0),
    indices_(0),
    allIndices_(0),
    pntr_(0),
    profile_(0),
    isUpperTriangular_(false),
    isLowerTriangular_(false),
    hasImplicitUnitDiagonal_(false),
    hasDiagonalEntries_(true)
{
}

//==============================================================================
template<typename Ordinal, typename Scalar>
HbMatrix<Ordinal, Scalar>::HbMatrix(const HbMatrix<Ordinal, Scalar> &matrix) 
  : numRows_(matrix.numRows_),
    numEntries_(matrix.numEntries_),
    values_(matrix.values_),
    allValues_(matrix.allValues_),
    indices_(matrix.indices_),
    allIndices_(matrix.allIndices_),
    pntr_(matrix.pntr_),
    profile_(matrix.profile_),
    isUpperTriangular_(matrix.isUpperTriangular_),
    isLowerTriangular_(matrix.isLowerTriangular_),
    hasImplicitUnitDiagonal_(matrix.hasImplicitUnitDiagonal_),
    hasDiagonalEntries_(matrix.hasDiagonalEntries_)

{
}

//==============================================================================
template<typename Ordinal, typename Scalar>
HbMatrix<Ordinal, Scalar>::~HbMatrix(){}

//==============================================================================
template<typename Ordinal, typename Scalar>
int HbMatrix<Ordinal, Scalar>::initializeStructure(Ordinal numRows, Ordinal * pntr, Ordinal * indx) { 
  numRows_ = numRows;
  pntr_ = pntr;
  allIndices_ = indx;

  numEntries_ = pntr_[numRows_];

  return(0);
}

//==============================================================================
template<typename Ordinal, typename Scalar>
int HbMatrix<Ordinal, Scalar>::initializeValues(Scalar * values) { 
  allValues_ = values;
  return(0);
}

//==============================================================================
template<typename Ordinal, typename Scalar>
int HbMatrix<Ordinal, Scalar>::getIndices(Ordinal i, Ordinal & numRowEntries, Ordinal *& indices) const { 
  numRowEntries = pntr_[i+1] - pntr_[i];
  indices = allIndices_ + pntr_[i];
  if (numRowEntries==0) return(1); // Warning, no nonzeros in this row/column
  if (numRowEntries<0) return(-1); // Fatal, nonsense data
  if (indices==0) return(-1); // Fatal, null pointer, but numRowEntries is positive
  return(0);
}


//==============================================================================
template<typename Ordinal, typename Scalar>
int HbMatrix<Ordinal, Scalar>::getValues(Ordinal i, Scalar *& values) const { 
  Ordinal numRowEntries;
  numRowEntries = pntr_[i+1] - pntr_[i];
  values = allValues_ + pntr_[i];
  if (numRowEntries==0) return(1); // Warning, no nonzeros in this row/column
  if (numRowEntries<0) return(-1); // Fatal, nonsense data
  if (values==0) return(-1); // Fatal, null pointer, but numRowEntries is positive
  return(0);
}


//==============================================================================
template<typename Ordinal, typename Scalar>
int HbMatrix<Ordinal, Scalar>::checkStructure() const { 

  bool isUpper = true;
  bool isLower = true;
  bool noDiag = true;
  int numOuter = numRows_;

  for (Ordinal i=0; i< numOuter; i++) {
    Ordinal j0 = pntr_[i];
    Ordinal j1 = pntr_[i+1];
    for (Ordinal jj=j0; jj<j1; jj++) {
      Ordinal j = allIndices_[jj];
      if (i<j) isUpper = false;
      if (i>j) isLower = false;
      if (i==j) noDiag = false;
    }
  } 
  int ierr = 0;
  // Test the values of upper, lower, and noDiag and make sure they are compatible with user-asserted values

  if (isUpperTriangular_ && !isUpper) { // User asserts upper triangular, but it's not
    ierr = -1;
    isUpperTriangular_ = isUpper;
  }
  else if (!isUpperTriangular_ && isUpper) { // User does not assert upper triangular, but it is
    ierr = 1;
    isUpperTriangular_ = isUpper;
  }
  if (isLowerTriangular_ && !isLower) { // User asserts lower triangular, but it's not
    ierr = -2;
    isLowerTriangular_ = isLower;
  }
  else if (!isLowerTriangular_ && isLower) { // User does not assert lower triangular, but it is
    ierr = 2;
    isLowerTriangular_ = isLower;
  }
  if (!hasDiagonalEntries_ && !noDiag) { // User asserts no diagonal, but values are given
    ierr = -3;
    hasDiagonalEntries_ = true;
  }
  return(ierr);
}

#endif /* KOKKOS_HBMATRIX_H */
