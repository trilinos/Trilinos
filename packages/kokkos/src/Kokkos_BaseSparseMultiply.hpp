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

#ifndef KOKKOS_BASESPARSEMULTIPLY_H
#define KOKKOS_BASESPARSEMULTIPLY_H

#include "Kokkos_CompObject.hpp" 


namespace Kokkos {

//! Kokkos::BaseSparseMultiply: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Kokkos::BaseSparseMultiply enables the piecewise construction and use of real-valued double-precision sparse matrices
  where matrix entries are intended for row access.

  At this time, the primary function provided by Kokkos::BaseSparseMultiply is matrix times vector and matrix 
  times multi-vector multiplication.  It is also possible to extract matrix rows from a constructed matrix.

  <b>Constructing Kokkos::BaseSparseMultiply objects</b>

  Constructing Kokkos::BaseSparseMultiply objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::BaseSparseMultiply instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
  </ol>

  Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
  create new entries.

  <b> Counting Floating Point Operations </b>

  Each Kokkos::BaseSparseMultiply object keeps track of the number
  of \e serial floating point operations performed using the specified object as the \e this argument
  to the function.  The Flops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate parallel performance
  numbers.  The ResetFlops() function resets the floating point counter.


*/    

  template<typename OrdinalType, typename ScalarType>
  class BaseSparseMultiply: public CompObject {
  public:

    //@{ \name Constructors/Destructor.
    //! BaseSparseMultiply constuctor with variable number of indices per row.
    BaseSparseMultiply();
  
    //! Copy constructor.
    BaseSparseMultiply(const BaseSparseMultiply& source);
	
    //! BaseSparseMultiply Destructor
    virtual ~BaseSparseMultiply();
    //@}
    //@{ \name Abstract Kokkos::CisMatrix Interface Initialization Methods
 
    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepStructure (In) If set to true, the user is asserting that the strucuture of the matrix, as
      defined in the getIndices() method of the CisMatrix object A will be kept.  Specifically, the pointer to an 
      array of indices returned for each i in that method will continue to point to valid index data.  By default,
      this argument is set to false, implying that the calling routine is \e not required to maintain the validity 
      of this data.  If the calling routine is planning to keep this data anyway, setting this argument to true can
      reduce the overall memory requirements.
      \return Integer error code, set to 0 if successful.
    */
    virtual int initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepStructure = false);
 
    //! Initialize values of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepValues (In) If set to true, the user is asserting that the strucuture of the matrix, as
      defined in the getIndices() method of the CisMatrix object A will be kept.  Specifically, the pointer to an 
      array of indices returned for each i in that method will continue to point to valid index data.  By default,
      this argument is set to false, implying that the calling routine is \e not required to maintain the validity 
      of this data.  If the calling routine is planning to keep this data anyway, setting this argument to true can
      reduce the overall memory requirements.
      \param checkStructure (In) If set to true, the structure of A will be checked against the structure of
      the matrix passed in to the initializeStructure() methods.  This parameter is false by default.
      \return Integer error code, set to 0 if successful.
    */
    virtual int initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepValues = false,
				 bool checkStructure = false);
 
    //@}

    //@{ \name Computational methods.
	
    //! Returns the result of a Kokkos_BaseSparseMultiply multiplied by a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to multiply by.
      \param y (Out) A Kokkos::Vector containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;

    //! Returns the result of a Kokkos_BaseSparseMultiply multiplied by multiple vectors in x, results in y.
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

    //! Returns true if this implementation of Kokkos::BaseSparseMultiply can benefit from the user keeping the passed in structure.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        analysis of structure.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepStructure parameter 
	passed in to the initializeStructure() method.
    */
    virtual bool getCanUseStructure() const {return(true);};

    //! Returns true if this implementation of Kokkos::BaseSparseMultiply can benefit from the user keeping the passed in values.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        copying of values.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepValues parameter 
	passed in to the initializeValues() method.
    */
    virtual bool getCanUseValues() const {return(true);};

    virtual const CisMatrix<OrdinalType, ScalarType> & getMatrix() const {
      if (matrixForValues_==0) return(*matrixForStructure_);
      else return(*matrixForValues_);
    };
	
    //@}
  
  protected:
    void copyStructure();
    void deleteStructure();
    void copyOrdinals(OrdinalType * vec, OrdinalType len);
    void copyScalars(ScalarType * vec, OrdinalType len);
    void copyValues();
    void deleteValues();

    CisMatrix<OrdinalType, ScalarType> * matrixForStructure_;
    CisMatrix<OrdinalType, ScalarType> * matrixForValues_;

    bool willKeepStructure_;
    bool willKeepValues_;
    bool isRowOriented_;
    bool haveStructure_;
    bool haveValues_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numRC_;
    OrdinalType numEntries_;

    ScalarType ** values_;

    OrdinalType ** indices_;
    OrdinalType * profile_;

  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseMultiply<OrdinalType, ScalarType>::BaseSparseMultiply() 
    : CompObject(),
      matrixForStructure_(0),
      matrixForValues_(0),
      willKeepStructure_(false),
      willKeepValues_(false),
      isRowOriented_(true),
      haveStructure_(false),
      haveValues_(false),
      numRows_(0),
      numCols_(0),
      numRC_(0),
      numEntries_(0),
      values_(0),
      indices_(0),
      profile_(0)  {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseMultiply<OrdinalType, ScalarType>::BaseSparseMultiply(const BaseSparseMultiply<OrdinalType, ScalarType> &source) 
    : CompObject(source),
      matrixForStructure_(source.matrixForStructure_),
      matrixForValues_(source.matrixForValues_),
      willKeepStructure_(source.willKeepStructure_),
      willKeepValues_(source.willKeepValues_),
      isRowOriented_(source.isRowOriented_),
      haveStructure_(source.haveStructure_),
      haveValues_(source.haveValues_),
      numRows_(source.numRows_),
      numCols_(source.numCols_),
      numRC_(source.numRC_),
      numEntries_(source.numEntries_),
      values_(source.values_),
      indices_(source.indices_),
      profile_(source.profile_) {

    copyStructure();
    copyValues();
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyStructure() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // pntr_, profile_, allIndices_, indices_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;

    if (profile_!=0) copyOrdinals(profiles_, numRC_);
    if (indices_!=0) {
      OrdinalType ** tmp_indices =indices_;
      indices_ = new OrdinalType*[numRC_];
      for (i=0; i< numRC_; i++) {
	indices_[i] = tmp_indices[i];
	if (!willKeepStructure_) 
	  copyOrdinals(indices_[i], profiles_[i]);
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyValues() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // allValues_, values_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;

    if (values_!=0) {
      ScalarType ** tmp_values =values_;
      values_ = new ScalarType*[numRC_];
      for (i=0; i< numRC_; i++) {
	values_[i] = tmp_values[i];
	if (!willKeepValues_) 
	copyScalars(values_[i], profiles_[i]);
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::deleteStructure() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // pntr_, profile_, allIndices_, indices_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;

    if (profile_!=0) delete [] profile_;
    if (indices_!=0) {
      if (!willKeepStructure_) 
	for (i=0; i< numRC_; i++) delete [] indices_[i];
      delete [] indices_;
      }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::deleteValues() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // allValues_, values_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;

    if (values_!=0) {
      if (!willKeepValues_) 
	for (i=0; i< numRC_; i++) delete [] values_[i];
      delete [] values_;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyOrdinals(OrdinalType * vec, OrdinalType len) {

    OrdinalType i;
    OrdinalType * tmp_pntr = vec;
    vec = new OrdinalType[len];
    for (i=0; i<len; i++) vec[i] = tmp_pntr[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyScalars(ScalarType * vec, OrdinalType len) {

    OrdinalType i;
    ScalarType * tmp_vals = vec;
    vec = new OrdinalType[len];
    for (i=0; i<len; i++) vec[i] = tmp_vals[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseMultiply<OrdinalType, ScalarType>::~BaseSparseMultiply(){

    deleteStructure();
    deleteValues();

  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseMultiply<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
								       bool willKeepStructure) {


    if (haveStructure_) return(-1); // Can only call this one time!

    matrixForStructure_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    willKeepStructure_ = willKeepStructure;
    isRowOriented_ = A.getIsRowOriented();
    numRows_ = A.getNumRows();
    numCols_ = A.getNumCols();
    numEntries_ = A.getNumEntries();
    if (isRowOriented_) numRC_ = numRows_;

    profile_ = new OrdinalType[numRC_];
    
    indices_ = new OrdinalType*[numRC_];

    OrdinalType numRCEntries;
    OrdinalType * indicesRC;

    for (i=0; i<numRC_; i++) {
      int ierr = A.getIndices(i, numRCEntries, indicesRC);
      if (ierr!=0) return(ierr);
      profile_[i] = numRCEntries;
      indices_[i] = indicesRC;
    }

    // If user will not keep structure, we must copy it
    if (!willKeepStructure_)
      for (i=0; i< numRC_; i++) {
	OrdinalType numIndices = profile_[i];
	OrdinalType * new_indices = new OrdinalType[numIndices];
	OrdinalType * old_indices = indices_[i];
	for (j=0; j<numIndices; j++) new_indices[j] = old_indices[j];
      }
    haveStructure_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseMultiply<OrdinalType, ScalarType>::initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, 
							       bool willKeepValues, bool checkStructure) {

    if (!haveStructure_) return(-1); // Must have structure first!

    matrixForValues_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    willKeepValues_ = willKeepValues;

    
    values_ = new ScalarType*[numRC_];

    ScalarType * valuesRC;

    for (i=0; i<numRC_; i++) {
      int ierr = A.getValues(i, valuesRC);
      if (ierr!=0) return(ierr);
      values_[i] = valuesRC;
    }

    // If user will not keep structure, we must copy it
    if (!willKeepValues_)
      for (i=0; i< numRC_; i++) {
	OrdinalType numIndices = profile_[i];
	ScalarType * new_values = new ScalarType[numIndices];
	ScalarType * old_values = values_[i];
	for (j=0; j<numIndices; j++) new_values[j] = old_values[j];
      }
    haveValues_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseMultiply<OrdinalType, ScalarType>::apply(const Vector<OrdinalType, ScalarType>& x, 
						    Vector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {

    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA) return(-2); // Unsupported at this time
    if (x.getLength()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getLength()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x				    

    OrdinalType i, j, curNumEntries;
    OrdinalType * curIndices;
    ScalarType * curValues;

    OrdinalType * profile = profile_;
    OrdinalType ** indices = indices_;
    ScalarType ** values = values_;

    ScalarType * xp = x.getValues();
    ScalarType * yp = y.getValues();

    if ((isRowOriented_ && !transA) ||
	(!isRowOriented_ && transA)) {

      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	curIndices = *indices++;
	curValues  = *values++;
	ScalarType sum = 0.0;
	for(j = 0; j < curNumEntries; j++)
	  sum += curValues[j] * xp[curIndices[j]];
	yp[i] = sum;
      }
    }
    else {
      
      for(i = 0; i < numRC_; i++)
	yp[i] = 0.0; // Initialize y for transpose multiply

      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	curIndices = *indices++;
	curValues  = *values++;
	for(j = 0; j < curNumEntries; j++)
	  yp[curIndices[j]] += curValues[j] * xp[i];
      }
    }
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseMultiply<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
						    MultiVector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {
    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA) return(-2); // Unsupported at this time
    if (x.getNumRows()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getNumRows()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x
    OrdinalType numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-5); // Not the same number of vectors in x and y

    OrdinalType i, j, k, curNumEntries;
    OrdinalType * curIndices;
    ScalarType * curValues;

    OrdinalType * profile = profile_;
    OrdinalType ** indices = indices_;
    ScalarType ** values = values_;

    ScalarType ** xp = x.getValues();
    ScalarType ** yp = y.getValues();

    for (k=0; k<numVectors; k++)
      for(i = 0; i < numRC_; i++)
	yp[k][i] = 0.0; // Initialize y

    if ((isRowOriented_ && !transA) ||
	(!isRowOriented_ && transA)) {

      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	curIndices = *indices++;
	curValues  = *values++;
	for (k=0; k<numVectors; k++) {
	  ScalarType sum = 0.0;
	  for(j = 0; j < curNumEntries; j++)
	    sum += curValues[j] * xp[k][curIndices[j]];
	  yp[k][i] = sum;
	}
      }
    }
    else {
      
      for (k=0; k<numVectors; k++)
	for(i = 0; i < numRC_; i++)
	  yp[k][i] = 0.0; // Initialize y
      
      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	curIndices = *indices++;
	curValues  = *values++;
	for (k=0; k<numVectors; k++) {
	  for(j = 0; j < curNumEntries; j++)
	    yp[k][curIndices[j]] += curValues[j] * xp[k][i];
	}
      }
    }
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_BASESPARSEMULTIPLY_H */
