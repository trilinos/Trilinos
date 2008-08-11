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
// ** should be able to use an OSKI matrix as a classic HB matrix
// ** Users can potentially extract the pointers they need, might be nice
// to have an "is a" relationship eventually.

#ifndef KOKKOS_OSKIMATRIX_HPP
#define KOKKOS_OSKIMATRIX_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CisMatrix.hpp"
extern "C" {
  #include <oski/oski.h>
}

namespace Kokkos {

//! Kokkos::OskiMatrix: Kokkos OSKI compressed index sparse matrix class.

/*! The Kokkos::OskiMatrix implements the Kokkos::CisMatrix interface. 
    SAY MORE HERE

*/    

  template<typename OrdinalType, typename ScalarType>
  class OskiMatrix: public virtual CisMatrix<OrdinalType,ScalarType> {
  public:

    //! @name Constructors/Destructor

    //@{

    //! Default OskiMatrix constuctor.
    OskiMatrix(void);
  
    //! Copy constructor.
    OskiMatrix(const OskiMatrix& source);

    //! HbMatrix Destructor
    virtual ~OskiMatrix();
    //@}
  
    //! @name OskiMatrix Initialization Methods

    //@{

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
      \param type (In) Matrix type.  OSKI supports: MAT_GENERAL (default), 
      MAT_TRI_UPPER, MAT_TRI_LOWER, MAT_SYMM_UPPER, MAT_SYMM_LOWER, 
      MAT_SYMM_FULL, MAT_HERM_UPPER, MAT_HERM_LOWER, MAT_HERM_FULL.
      \param diagExplicit (In) Non-zero diagonal entries are specified 
      explicitly (default true).  False implies all diagonal entries are 
      implicitly zero.
      \param indexBase (In) Array indicies can start at zero (default) or one.
      \param indexUnsorted (In) Indices within each row (column) appear in any
      order (default: true).  False indicates indices are sorted in increasing 
      order.
      \param indexRepeated (In) Indicies may appear multiple times (default:
      true).  False indicates all indices are unique.

      \return Integer error code, set to 0 if successful.
    */
    int initializeStructure(OrdinalType numRows, OrdinalType numCols, 
	bool isRowOriented, OrdinalType * pntr, OrdinalType * indx, 
	oski_inmatprop_t type=MAT_GENERAL, bool diagExplicit=true, int indexBase=0, 
	bool indexUnsorted=true, bool indexRepeated=true);

    //! Initialize matrix values
    /*!
      \param values (In)  Packed array of matrix values. values[pntr[i]] contains the first entry of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(ScalarType * values);
	
    //@}

    //! @name Matrix Format Initialization Methods

    //@{
 
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
/*    int initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
			    OrdinalType * profile, OrdinalType ** indx);
*/ 
    //! Initialize structure of matrix
    /*!
      \param values (In)  An array of pointers to arrays of matrix values. values[[i][0] contains the first entry of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(ScalarType ** values);
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
	
    //! Returns true if the compressed index matrix was formed using an OSKI CSR or CSC matrix.
//    bool getIsOskiMatrix() const {return(isOskiMatrix_);};
	
    //! Returns true if the compressed index matrix should be interpreted as a row matrix.
    bool getIsRowOriented() const {return(isRowOriented_);};
	
    //! Returns true if the compressed index matrix has no entries below the diagonal.
    virtual bool getIsUpperTriangular() const {return(matrixType_==MAT_TRI_UPPER);};
	
    //! Returns true if the compressed index matrix has no entries above the diagonal.
    virtual bool getIsLowerTriangular() const {return(matrixType_==MAT_TRI_LOWER);};
	
    //! Returns true if the compressed index matrix has no diagonal entries, but should be treated as unit diagonal.
    virtual bool getHasImplicitUnitDiagonal() const {return(diagExplicit_==MAT_UNIT_DIAG_IMPLICIT);};
	
    //! Number of rows
    OrdinalType getNumRows() const {return(numRows_);};
	
    //! Number of columns
    OrdinalType getNumCols() const {return(numCols_);};
	
    //! Number of matrix entries
    OrdinalType getNumEntries() const {return(numEntries_);};
	
    //! Underlying OSKI Matrix
    oski_matrix_t getA_tunable() const{return(A_tunable_);};

    //@}

  protected:
    oski_matrix_t A_tunable_;
    bool dataInitialized_;

//remove the ones with a "No" later and replace with OSKI equiv.
    OrdinalType numRows_;//No
    OrdinalType numCols_;//No
    OrdinalType numEntries_;//Yes?

    //ScalarType ** values_;//Yes
    ScalarType * allValues_;//Yes

    OrdinalType ** indices_;//Yes?
    OrdinalType * allIndices_;//Yes?

    OrdinalType * pntr_;//No
    //OrdinalType * profile_;//Yes?
    bool isRowOriented_;//Yes?
    //mutable bool isUpperTriangular_;//No? **change these to look at MatrixType
    //mutable bool isLowerTriangular_;//No?
    //bool hasImplicitUnitDiagonal_;//No?
    mutable bool hasDiagonalEntries_;//No?

    oski_inmatprop_t matrixType_;
    oski_inmatprop_t diagExplicit_;
    oski_inmatprop_t indexBase_;
    oski_inmatprop_t indexUnsorted_;
    oski_inmatprop_t indexRepeated_;

   bool isOskiMatrix_;//No? - no equiv, create dataInitialized_ 
  };

} // namespace Kokkos

using namespace Kokkos;
//==============================================================================
template<typename OrdinalType, typename ScalarType>
 OskiMatrix<OrdinalType, ScalarType>::OskiMatrix() 
  : numRows_(0),
    numCols_(0),
    numEntries_(0),
    //values_(0),
    allValues_(0),
    //indices_(0),
    allIndices_(0),
    pntr_(0),
    //profile_(0),
    isRowOriented_(true),
    //isUpperTriangular_(false),
    //isLowerTriangular_(false),
    //hasImplicitUnitDiagonal_(false),
    hasDiagonalEntries_(true),
    dataInitialized_(false)
{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
OskiMatrix<OrdinalType, ScalarType>::OskiMatrix(const OskiMatrix<OrdinalType, ScalarType> &matrix) 
  : numRows_(matrix.numRows_),
    numCols_(matrix.numCols_),
    numEntries_(matrix.numEntries_),
    //values_(matrix.values_),
    allValues_(matrix.allValues_),
    indices_(matrix.indices_),
    allIndices_(matrix.allIndices_),
    pntr_(matrix.pntr_),
    //profile_(matrix.profile_),
    isRowOriented_(matrix.isRowOriented_),
    //isUpperTriangular_(matrix.isUpperTriangular_),
    //isLowerTriangular_(matrix.isLowerTriangular_),
    //hasImplicitUnitDiagonal_(matrix.hasImplicitUnitDiagonal_),
    hasDiagonalEntries_(matrix.hasDiagonalEntries_),
    dataInitialized_(matrix.dataInitialized_)
{
    if (dataInitialized_) 
      A_tunable_ = oski_CopyMat(matrix.A_tunable_);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
OskiMatrix<OrdinalType, ScalarType>::~OskiMatrix(){
    if (dataInitialized_) oski_DestroyMat(A_tunable_);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::initializeStructure(OrdinalType numRows, 
	OrdinalType numCols, bool isRowOriented, OrdinalType * pntr, 
	OrdinalType * indx, oski_inmatprop_t type, bool diagExplicit, 
	int indexBase, bool indexUnsorted, bool indexRepeated) { 
  
  if (dataInitialized_) return(-1);
  numRows_ = numRows;
  numCols_ = numCols;
  isRowOriented_ = isRowOriented;
  pntr_ = pntr;
  allIndices_ = indx;
  if (type == MAT_GENERAL || type == MAT_TRI_UPPER || 
	type == MAT_TRI_LOWER || type == MAT_SYMM_UPPER || 
	type == MAT_SYMM_LOWER || type == MAT_SYMM_FULL ||
	type == MAT_HERM_UPPER || type == MAT_HERM_LOWER || 
	type == MAT_HERM_FULL)
    matrixType_ = type;
  else
    return (-1);

  if (diagExplicit)
    diagExplicit_ = MAT_DIAG_EXPLICIT;
  else
    diagExplicit_ = MAT_UNIT_DIAG_IMPLICIT;

  if (indexBase==1)
    indexBase_ = INDEX_ONE_BASED;
  else if (indexBase==0)
    indexBase_ = INDEX_ZERO_BASED;
  else
    return (-2);

  if (indexUnsorted)
    indexUnsorted_ = INDEX_UNSORTED;
  else
    indexUnsorted_ = INDEX_SORTED;

  if (indexRepeated)
    indexRepeated_ = INDEX_REPEATED;
  else
    indexRepeated_ = INDEX_UNIQUE;

  if (isRowOriented_) numEntries_ = pntr_[numRows_];
  else numEntries_ = pntr_[numCols_];
  
  //isOskiMatrix_ = true;//eventually this should be removed
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::initializeValues(ScalarType * values) { 
  //if (!isOskiMatrix_) return(-1); // initialization methods must match
  if (dataInitialized_) return(-1);//Cannot reinitialize values using this 
					// method at this time

  allValues_ = values;
  //OSKI Matrix constructor
 
  if (isRowOriented_) 
    A_tunable_ = oski_CreateMatCSR(pntr_,allIndices_,allValues_,numRows_,numCols_,SHARE_INPUTMAT,5,matrixType_,diagExplicit_,indexBase_,indexUnsorted_,indexRepeated_);
  else
    A_tunable_ = oski_CreateMatCSC(pntr_,allIndices_,allValues_,numRows_,numCols_,SHARE_INPUTMAT,5,matrixType_,diagExplicit_,indexBase_,indexUnsorted_,indexRepeated_);

  if (A_tunable_ == NULL) return(-2);
  dataInitialized_ = true;
  return(0);
}
/* We aren't able to support general hb format at this time
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
*/
//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OskiMatrix<OrdinalType, ScalarType>::getIndices(OrdinalType i, OrdinalType & numRowEntries, OrdinalType *& indices) const { 
  // Oski version
  
//  if (isOskiMatrix_) {
    numRowEntries = pntr_[i+1] - pntr_[i];
    indices = allIndices_ + pntr_[i];
//  }
//  else {
//    numRowEntries = profile_[i];
//    indices = indices_[i];
//  }

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
//  if (isOskiMatrix_) {
    numRowEntries = pntr_[i+1] - pntr_[i];
    values = allValues_ + pntr_[i];
//  }
//  else {
//    numRowEntries = profile_[i];
//    values = values_[i];
//  }

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

//  if (isClassicHbMatrix_) {
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
//  }
//  else {
//    for (OrdinalType i=0; i< numOuter; i++) {
//      OrdinalType * curIndices = indices_[i];
//      OrdinalType j1 = pntr_[i+1];
//      for (OrdinalType jj=0; jj<j1; jj++) {
//	OrdinalType j = curIndices[jj];
//	if (i<j) isUpper = false;
//	if (i>j) isLower = false;
//	if (i==j) noDiag = false;
//      }
//    } 
//  }
  int ierr = 0;
  // Test the values of upper, lower, and noDiag and make sure they are compatible with user-asserted values

  if (matrixType_==MAT_TRI_UPPER && !isUpper) { // User asserts upper triangular, but it's not
    ierr = -1;
  }
  else if (!(matrixType_==MAT_TRI_UPPER || matrixType_==MAT_SYMM_UPPER) && isUpper) { // User does not assert upper triangular, but it is
    ierr = 1;
    //May not be able to change at this point if values are initiated
  }
  if (matrixType_==MAT_TRI_LOWER && !isLower) { // User asserts lower triangular, but it's not
    ierr = -2;
  }
  else if (!(matrixType_==MAT_TRI_LOWER || matrixType_==MAT_SYMM_LOWER) && isLower) { // User does not assert lower triangular, but it is
    ierr = 2;
    //May not be able to change at this point if values are initiated
  }
  if (!hasDiagonalEntries_ && !noDiag) { // User asserts no diagonal, but values are given
    ierr = -3;
    hasDiagonalEntries_ = true;
  }
  return(ierr);
}

#endif /* KOKKOS_OSKIMATRIX_HPP */

