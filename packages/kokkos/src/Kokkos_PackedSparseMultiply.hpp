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

#ifndef KOKKOS_PACKEDSPARSEMULTIPLY_H
#define KOKKOS_PACKEDSPARSEMULTIPLY_H

#include "Kokkos_BaseSparseMultiply.hpp" 


namespace Kokkos {

//! Kokkos::BaseSparseMultiply: A reference class for computing sparse matrix multiplication operations.

/*! The Kokkos::BaseSparseMultiply provide basic functionality for computing sparse matrix times vector, or
    sparse matrix times multivector operations.  This class is templated on the ordinal (integer) and scalar (floating
    point) types, so it can compute using any reasonable data type.

  <b>Constructing Kokkos::BaseSparseMultiply objects</b>

  Constructing Kokkos::BaseSparseMultiply objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::BaseSparseMultiply instance:  The constructor takes no arguments.
  <li> Register the structure of a Kokkos::CisMatrix object using initializeStructure(): 
       We provide this method so that derived implementations can
       take advantage of multiple problems that have the same structure.  In this situation, initializeStructure() would
       be called once and then initializeValues() would be called repeatedly, amortizing the cost of setting up the structure.
       This method may be called only once.
  <li> Register the values of a Kokkos::CisMatrix object using initializeValues(): This method is used to pass values to the
       multiply class.  It can be called repeatedly if multiple matrices have the same structure.
  </ol>

  <b> Counting Floating Point Operations </b>

  Each Kokkos::BaseSparseMultiply object keeps track of the number
  of floating point operations performed using the specified object as the \e this argument
  to the function.  The getFlops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate  performance
  numbers.  The resetFlops() function resets the floating point counter.

*/    

  template<typename OrdinalType, typename ScalarType>
  class PackedSparseMultiply: public virtual BaseSparseMultiply<OrdinalType, ScalarType> {
  public:

    //@{ \name Constructors/Destructor.
    //! PackedSparseMultiply constuctor with variable number of indices per row.
    PackedSparseMultiply();
  
    //! Copy constructor.
    PackedSparseMultiply(const PackedSparseMultiply& source);
	
    //! PackedSparseMultiply Destructor
    virtual ~PackedSparseMultiply();
    //@}
    //@{ \name Abstract Kokkos::CisMatrix Interface Initialization Methods
 
    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepStructure (In) This argument is unused by this implementation of the BaseSparseMultiply
             class since structure and values will be copied.
      \return Integer error code, set to 0 if successful.
    */
    virtual int initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepStructure = false);
 
    //! Initialize values of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepValues (In) This argument is unused by this implementation of the BaseSparseMultiply
             class since structure and values will be copied.
      \param checkStructure (In) If set to true, the structure of A will be checked against the structure of
      the matrix passed in to the initializeStructure() methods.  This parameter is false by default.

      \return Integer error code, set to 0 if successful, returns - 1 if checkStructure is true and structure is changed.
    */
    virtual int initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepValues = false,
				 bool checkStructure = false);
 
    //@}

    //@{ \name Computational methods.
	
    //! Returns the result of a Kokkos_PackedSparseMultiply multiplied by a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to multiply by.
      \param y (Out) A Kokkos::Vector containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;

    //! Returns the result of a Kokkos_PackedSparseMultiply multiplied by multiple vectors in x, results in y.
    /*! 
      \param x (In) A Kokkos::MultiVector to multiply by.
      \param y (Out) A Kokkos::MultiVector containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;
    //@}
	
    //@{ \name Operator attribute access methods.

    //! Returns false for this implementation.
    /*! This implementation will not use the user's copy of the matrix structure.
    */
    virtual bool getCanUseStructure() const {return(false);};

    //! Returns false for this implementation.
    /*! This implementation will not use the user's copy of the matrix values.
    */
    virtual bool getCanUseValues() const {return(false);};
	
    //@}
  
  protected:

    void copyEntries();
    void deleteStructureAndValues();

    struct EntryStruct {
      OrdinalType index;
      ScalarType value;
    }; 
    typedef struct EntryStruct Entry;
    
    Entry * allEntries_;
  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  PackedSparseMultiply<OrdinalType, ScalarType>::PackedSparseMultiply() 
    : BaseSparseMultiply<OrdinalType, ScalarType>(),
      allEntries_(0) {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  PackedSparseMultiply<OrdinalType, ScalarType>::PackedSparseMultiply(const PackedSparseMultiply<OrdinalType, ScalarType> &source) 
    : BaseSparseMultiply<OrdinalType, ScalarType>(source),
      allEntries_(source.allEntries_) {

    copyEntries();
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void PackedSparseMultiply<OrdinalType, ScalarType>::copyEntries() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // pntr_, profile_, allIndices_, indices_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;

    if (allEntries_!=0) {
      Entry * tmp_entries = new Entry[numEntries_];
      for (i=0; i< numEntries_; i++) {
	tmp_entries[i].index = allEntries_[i].index;
      	tmp_entries[i].value = allEntries_[i].value;
      }
      allEntries_ = tmp_entries;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void PackedSparseMultiply<OrdinalType, ScalarType>::deleteStructureAndValues() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // pntr_, profile_, allIndices_, indices_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;

    if (profile_!=0) {
      delete [] profile_;
      profile_ = 0;
    }

    if (allEntries_!=0) {
      delete [] allEntries_;
      allEntries_ = 0;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  PackedSparseMultiply<OrdinalType, ScalarType>::~PackedSparseMultiply(){

    deleteStructureAndValues();

  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
								       bool willKeepStructure) {


    if (haveStructure_) return(-1); // Can only call this one time!

    matrixForStructure_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    willKeepStructure_ = willKeepStructure;
    isRowOriented_ = A.getIsRowOriented();
    numRows_ = A.getNumRows();
    numCols_ = A.getNumCols();
    numEntries_ = A.getNumEntries();
    numRC_ = numCols_;
    if (isRowOriented_) numRC_ = numRows_;

    profile_ = new OrdinalType[numRC_];

    OrdinalType numRCEntries;
    OrdinalType * indicesRC;

      
    allEntries_ = new Entry[numEntries_]; // Allocate storage for all entries at once
    
    OrdinalType offset = 0;
    for (i=0; i< numRC_; i++) {
      int ierr = A.getIndices(i, numRCEntries, indicesRC);
      if (ierr!=0) return(ierr);
      profile_[i] = numRCEntries;
      Entry * curRC = allEntries_+offset;
      for (j=0; j<numRCEntries; j++) curRC[j].index = indicesRC[j];
      offset += numRCEntries;
    }

    costOfMatVec_ = 2.0 * ((double) numEntries_);
    haveStructure_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, 
							       bool willKeepValues, bool checkStructure) {

    if (!haveStructure_) return(-1); // Must have structure first!

    matrixForValues_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    willKeepValues_ = willKeepValues;

    ScalarType * valuesRC;

    OrdinalType offset = 0;
    for (i=0; i<numRC_; i++) {
      int ierr = A.getValues(i, valuesRC);
      if (ierr!=0) return(ierr);
      Entry * curRC = allEntries_+offset;
      OrdinalType numRCEntries = profile_[i];
      for (j=0; j<numRCEntries; j++) curRC[j].value = valuesRC[j];
      offset += numRCEntries;
    }
    haveValues_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::apply(const Vector<OrdinalType, ScalarType>& x, 
						    Vector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {

    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA) return(-2); // Unsupported at this time
    if (x.getLength()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getLength()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x				    

    OrdinalType i, j, curNumEntries;
    Entry * curEntries = allEntries_;

    OrdinalType * profile = profile_;

    ScalarType * xp = x.getValues();
    ScalarType * yp = y.getValues();

    if ((isRowOriented_ && !transA) ||
	(!isRowOriented_ && transA)) {

      
      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	ScalarType sum = 0.0;
	for(j = 0; j < curNumEntries; j++)
	  sum += curEntries[j].value * xp[curEntries[j].index];
	yp[i] = sum;
	curEntries += curNumEntries;
      }
    }
    else {
      
      for(i = 0; i < numRC_; i++)
	yp[i] = 0.0; // Initialize y for transpose multiply

      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	for(j = 0; j < curNumEntries; j++)
	  yp[curEntries[j].index] += curEntries[j].value * xp[i];
	curEntries += curNumEntries;	
      }
    }
    updateFlops(costOfMatVec_);
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
						    MultiVector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {
    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA) return(-2); // Unsupported at this time
    if (x.getNumRows()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getNumRows()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x
    OrdinalType numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-5); // Not the same number of vectors in x and y

    OrdinalType i, j, k, curNumEntries;
    Entry * curEntries = allEntries_;

    OrdinalType * profile = profile_;

    ScalarType ** xpp = x.getValues();
    ScalarType ** ypp = y.getValues();

    if ((isRowOriented_ && !transA) ||
	(!isRowOriented_ && transA)) {

      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	for (k=0; k<numVectors; k++) {
	  ScalarType * xp = xpp[k];
	  ScalarType * yp = ypp[k];
	  ScalarType sum = 0.0;
	  for(j = 0; j < curNumEntries; j++)
	    sum += curEntries[j].value * xp[curEntries[j].index];
	  yp[i] = sum;
	}
	curEntries += curNumEntries;
      }
    }
    else {
      
      for (k=0; k<numVectors; k++) {
	ScalarType * yp = ypp[k];
	for(i = 0; i < numRC_; i++)
	  yp[i] = 0.0; // Initialize y
      }
      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	for (k=0; k<numVectors; k++) {
	  ScalarType * xp = xpp[k];
	  ScalarType * yp = ypp[k];
	  for(j = 0; j < curNumEntries; j++)
	    yp[curEntries[j].index] += curEntries[j].value * xp[i];
	}
	curEntries += curNumEntries;
      }
    }
    updateFlops(costOfMatVec_ * ((double) numVectors));
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_PACKEDSPARSEMULTIPLY_H */
