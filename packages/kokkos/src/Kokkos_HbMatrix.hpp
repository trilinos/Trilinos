//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_HBMATRIX_H
#define KOKKOS_HBMATRIX_H
#include "Kokkos_CisMatrix.hpp"

namespace Kokkos {



//! Kokkos::HbMatrix: Kokkos compressed index sparse matrix base class.

/*! The Kokkos::HbMatrix implements the Kokkos::CisMatrix interface 
    using either a Harwell-Boeing matrix or
    generalized form of one.

*/    

  template<typename OrdinalType, typename ScalarType>
  class HbMatrix: public virtual CisMatrix<OrdinalType,ScalarType> {
  public:

    //@{ \name Constructors/Destructor.

    //! Default HbMatrix constuctor.
    HbMatrix(void);
  
    //! Copy constructor.
    HbMatrix(const HbMatrix& source);

    //! HbMatrix Destructor
    virtual ~HbMatrix();
    //@}
  
    //@{ \name Classical Harwell-Boeing Format Initialization Methods
	
    //! Initialize structure of matrix (Classical Harwell-Boeing Format)
    /*!
      This interface supports matrices that are stored in the classical Harwell-Boeing format.
      \param numRows (In)  Row dimension.
      \param numCols (In)  Column dimension.
      \param isRowOriented - If true, the compressed index storage will be interpreted as row indices.  
      If false, then indices will be interpreted as column indices.
      \param pntr (In)  Array of offsets into indx.  indx[pntr[i]] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \param indx (In)  Packed array of indices.  indx[pntr[i]] contains the first index of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
			    OrdinalType * pntr, OrdinalType * indx);

    //! Initialize structure of matrix (Classical Harwell-Boeing Format)
    /*!
      This interface supports matrices that are stored in the classical Harwell-Boeing format.
      \param values (In)  Packed array of matrix values. values[pntr[i]] contains the first entry of the ith row
      (if isRowOriented is true) or ith column (if isRowOriented is false).
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(ScalarType * values);
	
    //@}

    //@{ \name Generalized Harwell-Boeing Format Initialization Methods
 
    //! Initialize structure of matrix (Generalized Harwell-Boeing Format)
    /*!
      This interface supports matrices that are stored in a generalized Harwell-Boeing format.
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
 
    //! Initialize structure of matrix (Generalized Harwell-Boeing Format)
    /*!
      This interface supports matrices that are stored in the classical Harwell-Boeing format.
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

    //@{ \name Matrix Attribute access methods.
	
    //! Returns true if the compressed index matrix was formed using a classic HB matrix.
    bool getIsClassicHbMatrix() const {return(isClassicHbMatrix_);};
	
    //! Returns true if the compressed index matrix should be interpreted as a row matrix.
    bool getIsRowOriented() const {return(isRowOriented_);};
	
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

    bool isClassicHbMatrix_;
  };

} // namespace Kokkos

using namespace Kokkos;
//==============================================================================
template<typename OrdinalType, typename ScalarType>
 HbMatrix<OrdinalType, ScalarType>::HbMatrix() 
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
    isClassicHbMatrix_(true)
{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
HbMatrix<OrdinalType, ScalarType>::HbMatrix(const HbMatrix<OrdinalType, ScalarType> &matrix) 
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
    isClassicHbMatrix_(matrix.isClassicHbMatrix_)

{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
HbMatrix<OrdinalType, ScalarType>::~HbMatrix(){}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int HbMatrix<OrdinalType, ScalarType>::initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
							   OrdinalType * pntr, OrdinalType * indx) { 
  // Classic HB version

  numRows_ = numRows;
  numCols_ = numCols;
  isRowOriented_ = isRowOriented;
  pntr_ = pntr;
  allIndices_ = indx;


  if (isRowOriented_) numEntries_ = pntr_[numRows_];
  else numEntries_ = pntr_[numCols_];
  
  isClassicHbMatrix_ = true;
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int HbMatrix<OrdinalType, ScalarType>::initializeValues(ScalarType * values) { 

  // Classic HB version

  if (!isClassicHbMatrix_) return(-1); // initialization methods must match

  allValues_ = values;
  
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int HbMatrix<OrdinalType, ScalarType>::initializeStructure(OrdinalType numRows, OrdinalType numCols, bool isRowOriented,
							   OrdinalType * profile, OrdinalType ** indx) { 
  // Generalized HB version
  
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
  
  isClassicHbMatrix_ = false;
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int HbMatrix<OrdinalType, ScalarType>::initializeValues(ScalarType ** values) { 
  // Generalized HB version
  
  if (isClassicHbMatrix_) return(-1); // initialization methods must match

  values_ = values;
  
  return(0);
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int HbMatrix<OrdinalType, ScalarType>::getIndices(OrdinalType i, OrdinalType & numRowEntries, OrdinalType *& indices) const { 
  // Generalized HB version
  
  if (isClassicHbMatrix_) {
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
int HbMatrix<OrdinalType, ScalarType>::getValues(OrdinalType i, ScalarType *& values) const { 
  // Generalized HB version
  
  OrdinalType numRowEntries;
  if (isClassicHbMatrix_) {
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

#endif /* KOKKOS_HBMATRIX_H */
