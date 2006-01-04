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

// Author: Jim Willenbring jmwille@sandia.gov 08-12-2005

#ifndef KOKKOS_OSKIMATRIX_HPP
#define KOKKOS_OSKIMATRIX_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CisMatrix.hpp"

namespace Kokkos {

//! Kokkos::OskiMatrix: Kokkos OSKI compressed index sparse matrix class.

/*! The Kokkos::OskiMatrix implements the Kokkos::CisMatrix interface. 
    SAY MORE HERE

*/    

  template<typename OrdinalType, typename ScalarType>
  class OskiMatrix: public virtual CisMatrix<OrdinalType,ScalarType> {
  public:

    //@{ \name Constructors/Destructor.

    //! Default OskiMatrix constuctor.
    OskiMatrix(void);
  
    //! Copy constructor.
    OskiMatrix(const OskiMatrix& source);

    //! HbMatrix Destructor
    virtual ~OskiMatrix();
    //@}
  
    //@{ \name OskiMatrix Initialization Methods
	
    //! Initialize structure of the matrix
    /*!
      \param numRows (In)  Row dimension.
      \param numCols (In)  Column dimension.
      \param isRowOriented - (In) If true, the compressed index storage will be interpreted as row indices.  
      If false, then indices will be interpreted as column indices.
      \param pntr (In)  Array of offsets into indx.  indx[pntr[i]] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \param indx (In)  Packed array of indices.  indx[pntr[i]] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
			    OrdinalType * pntr, OrdinalType * indx);

    //! Initialize structure of matrix
    /*!
      \param values (In)  Packed array of matrix values. values[pntr[i]] contains the first entry of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(ScalarType * values);
	
    //@}

    //@{ \name Matrix Format Initialization Methods
 
    //! Initialize structure of matrix
    /*!
      \param numRows (In)  Row dimension.
      \param numCols (In)  Column dimension.
      \param isRowOriented (In) If true, the compressed index storage will be interpreted as row indices.
      If false, then indices will be interpreted as column indices.
      \param profile (In)  Array of index counts for indx.  pntr[i] equals the number of entries in the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \param indx (In)  An array of pointers to arrays of indices.  indx[i][0] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
			    OrdinalType * profile, OrdinalType ** indx);
 
    //! Initialize structure of matrix
    /*!
      \param values (In)  An array of pointers to arrays of matrix values. values[[i][0] contains the first entry of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(ScalarType ** values);
    //@}

    //@{ \name Matrix entry access methods.

    //! Returns number of entries in ith row/column, and pointer to an array of these indices.
    /*! Extract the number of entries and a pointer to the indices in the ith row/column of the matrix.  Note that
        the indices are not copied by this method.  Memory allocation is handled by the matrix object itself.

	\param i (In) The row (if isRowOriented() is true) or column that should be returned.
	\param numRowEntries (Out) The number of entries in the ith row/column.
	\param indices (Out) A pointer to the list of indices in the ith row/column.

	\return Integer error code, set to 0 if successful.
    */
    int getIndices(OrdinalType i, OrdinalType & numRowEntries, OrdinalType *& indices) const;

    //! Returns a pointer to an array of values for the ith row/column.
    /*! Extract the values in the ith row/column of the matrix.  Note that
        the values are not copied by this method.  Memory allocation is handled by the matrix object itself.

	\param i (In) The row (if isRowOriented() is true) or column that should be returned.
	\param numEntries (Out) The number of entries in the ith row/column.
	\param indices (Out) A pointer to the list of indices in the ith row/column.

	\return Integer error code, set to 0 if successful.
    */
    int getValues(OrdinalType i, ScalarType *& values) const;
	
	
    //@}

    //@{ \name Matrix Attribute set methods.

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

    //@{ \name Validity tests.

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

    //@{ \name Matrix Attribute access methods.
	
    //! Returns true if the compressed index matrix was formed using an OSKI CSR or CSC matrix.
    bool getIsOskiMatrix() const {return(isOskiMatrix_);};
	
    //! Returns true if the compressed index matrix should be interpreted as a row matrix.
    bool getIsRowOriented() const {return(isRowOriented_);};
	
    //! Returns true if the compressed index matrix has no entries below the diagonal.
    virtual bool getIsUpperTriangular() const {return(isUpperTriangular_);};
	
    //! Returns true if the compressed index matrix has no entries above the diagonal.
    virtual bool getIsLowerTriangular() const {return(isLowerTriangular_);};
	
    //! Returns true if the compressed index matrix has no diagonal entries, but should be treated as unit diagonal.
    virtual bool getHasImplicitUnitDiagonal() const {return(hasImplicitUnitDiagonal_);};
	
    //! Number of rows
    OrdinalType getNumRows() const {return(numRows_);};
	
    //! Number of columns
    OrdinalType getNumCols() const {return(numCols_);};
	
    //! Number of matrix entries
    OrdinalType getNumEntries() const {return(numEntries_);};
	
    //@}

  protected:

    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numEntries_;

    ScalarType ** values_;
    ScalarType * allValues_;

    OrdinalType ** indices_;
    OrdinalType * allIndices_;

    OrdinalType * pntr_;
    OrdinalType * profile_;
    bool isRowOriented_;
    mutable bool isUpperTriangular_;
    mutable bool isLowerTriangular_;
    bool hasImplicitUnitDiagonal_;
    mutable bool hasDiagonalEntries_;

    bool isOskiMatrix_;
  };

} // namespace Kokkos

using namespace Kokkos;
//==============================================================================
template<typename OrdinalType, typename ScalarType>
 OskiMatrix<OrdinalType, ScalarType>::OskiMatrix() 
  : numRows_(0),
    numCols_(0),
    numEntries_(0),
    values_(0),
    allValues_(0),
    indices_(0),
    allIndices_(0),
    pntr_(0),
    profile_(0),
    isRowOriented_(true),
    isUpperTriangular_(false),
    isLowerTriangular_(false),
    hasImplicitUnitDiagonal_(false),
    hasDiagonalEntries_(true),
    isOskiMatrix_(true)
{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
OskiMatrix<OrdinalType, ScalarType>::OskiMatrix(const OskiMatrix<OrdinalType, ScalarType> &matrix) 
  : numRows_(matrix.numRows_),
    numCols_(matrix.numCols_),
    numEntries_(matrix.numEntries_),
    values_(matrix.values_),
    allValues_(matrix.allValues_),
    indices_(matrix.indices_),
    allIndices_(matrix.allIndices_),
    pntr_(matrix.pntr_),
    profile_(matrix.profile_),
    isRowOriented_(matrix.isRowOriented_),
    isUpperTriangular_(matrix.isUpperTriangular_),
    isLowerTriangular_(matrix.isLowerTriangular_),
    hasImplicitUnitDiagonal_(matrix.hasImplicitUnitDiagonal_),
    hasDiagonalEntries_(matrix.hasDiagonalEntries_),
    isOskiMatrix_(matrix.isOskiMatrix_)

{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
OskiMatrix<OrdinalType, ScalarType>::~OskiMatrix(){}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
							   OrdinalType * pntr, OrdinalType * indx) { 
  // Oski version

  numRows_ = numRows;
  numCols_ = numCols;
  isRowOriented_ = isRowOriented;
  pntr_ = pntr;
  allIndices_ = indx;


  if (isRowOriented_) numEntries_ = pntr_[numRows_];
  else numEntries_ = pntr_[numCols_];
  
  isOskiMatrix_ = true;
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::initializeValues(ScalarType * values) { 
  // Oskversion

  if (!isOskiMatrix_) return(-1); // initialization methods must match

  allValues_ = values;
  
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
							   OrdinalType * profile, OrdinalType ** indx) { 
  // Oski version
  
  numRows_ = numRows;
  numCols_ = numCols;
  isRowOriented_ = isRowOriented;
  profile_ = profile;
  indices_ = indx;
  
  numEntries_ = 0;
  if (isRowOriented_) 
    for (int i=0; i<numRows_; i++)
      numEntries_ += profile_[i];
  else 
    for (int i=0; i<numCols_; i++)
      numEntries_ += profile_[i];
  
  isOskiMatrix_ = false;
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::initializeValues(ScalarType ** values) { 
  // Oski version
  
  if (isOskiMatrix_) return(-1); // initialization methods must match

  values_ = values;
  
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::getIndices(OrdinalType i, OrdinalType & numRowEntries, OrdinalType *& indices) const { 
  // Oski version
  
  if (isOskiMatrix_) {
    numRowEntries = pntr_[i+1] - pntr_[i];
    indices = allIndices_ + pntr_[i];
  }
  else {
    numRowEntries = profile_[i];
    indices = indices_[i];
  }

  if (numRowEntries==0) return(1); // Warning, no nonzeros in this row/column
  if (numRowEntries<0) return(-1); // Fatal, nonsense data
  if (indices==0) return(-1); // Fatal, null pointer, but numRowEntries is positive
  
  return(0);
}


//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::getValues(OrdinalType i, ScalarType *& values) const { 
  // Oski version
  
  OrdinalType numRowEntries;
  if (isOskiMatrix_) {
    numRowEntries = pntr_[i+1] - pntr_[i];
    values = allValues_ + pntr_[i];
  }
  else {
    numRowEntries = profile_[i];
    values = values_[i];
  }

  if (numRowEntries==0) return(1); // Warning, no nonzeros in this row/column
  if (numRowEntries<0) return(-1); // Fatal, nonsense data
  if (values==0) return(-1); // Fatal, null pointer, but numRowEntries is positive
  
  return(0);
}


//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::checkStructure() const { 

  bool isUpper = true;
  bool isLower = true;
  bool noDiag = true;
  int numOuter = numCols_;
  if (isRowOriented_) numOuter = numRows_;

  if (isClassicHbMatrix_) {
    for (OrdinalType i=0; i< numOuter; i++) {
      OrdinalType j0 = pntr_[i];
      OrdinalType j1 = pntr_[i+1];
      for (OrdinalType jj=j0; jj<j1; jj++) {
	OrdinalType j = allIndices_[jj];
	if (i<j) isUpper = false;
	if (i>j) isLower = false;
	if (i==j) noDiag = false;
      }
    } 
  }
  else {
    for (OrdinalType i=0; i< numOuter; i++) {
      OrdinalType * curIndices = indices_[i];
      OrdinalType j1 = pntr_[i+1];
      for (OrdinalType jj=0; jj<j1; jj++) {
	OrdinalType j = curIndices[jj];
	if (i<j) isUpper = false;
	if (i>j) isLower = false;
	if (i==j) noDiag = false;
      }
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

#endif /* KOKKOS_OSKIMATRIX_HPP */
