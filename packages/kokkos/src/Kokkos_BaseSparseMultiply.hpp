//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
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

#ifndef KOKKOS_BASESPARSEMULTIPLY_H
#define KOKKOS_BASESPARSEMULTIPLY_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_Permutation.hpp"
#include "Kokkos_CisMatrix.hpp" 
#include "Kokkos_SparseOperation.hpp" 
#include "Kokkos_Fortran_wrappers.hpp"

namespace Kokkos {

//! Kokkos::BaseSparseMultiply: A reference class for computing sparse matrix multiplication operations.

/*! The Kokkos::BaseSparseMultiply provide basic functionality for computing sparse matrix times vector, or
    sparse matrix times multivector operations.  This class is templated on the ordinal (integer) and 
    scalar (floating point) types, so it can compute using any reasonable data type.  It 
    implements the Kokkos::SparseOperation base class.

  <b>Constructing Kokkos::BaseSparseMultiply objects</b>

  Constructing Kokkos::BaseSparseMultiply objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::BaseSparseMultiply instance:  The constructor takes no arguments.
  <li> Register the structure of a Kokkos::CisMatrix object using initializeStructure(): 
       We provide this method so that derived implementations can
       take advantage of multiple problems that have the same structure.  In this situation, 
       initializeStructure() would
       be called once and then initializeValues() would be called repeatedly, amortizing the cost of 
       setting up the structure.
       This method may be called only once.
  <li> Register the values of a Kokkos::CisMatrix object using initializeValues(): This method is 
       used to pass values to the
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
  class BaseSparseMultiply: public virtual SparseOperation<OrdinalType, ScalarType> {
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
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface
.
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

      \return Integer error code, set to 0 if successful, returns - 1 if checkStructure is true and structure is changed.
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

    //! Returns a reference to the most recent Kokkos::CisMatrix that was passed into the \e this object.
    virtual const CisMatrix<OrdinalType, ScalarType> & getMatrix() const {
      if (matrixForValues_==0) return(*matrixForStructure_);
      else return(*matrixForValues_);};

    //! Returns a reference to the left Kokkos::Permutation object, which is the identity for this implementation.
    virtual const Permutation<OrdinalType, ScalarType> & getLeftPermutation() const {
      return(leftPermutation_);};

    //! Returns a reference to the right Kokkos::Permutation object, which is the identity for this implementation.
    virtual const Permutation<OrdinalType, ScalarType> & getRightPermutation() const {
      return(rightPermutation_);};
	
    //@}
  
  protected:
    void copyProfile();
    void copyStructure();
    void deleteStructureAndProfile();
    void copyOrdinals(OrdinalType len, OrdinalType * vecIn, OrdinalType * vecOut);
    void copyScalars(OrdinalType len, ScalarType * vecIn, ScalarType * vecOut);
    void copyValues();
    void deleteValues();

    CisMatrix<OrdinalType, ScalarType> * matrixForStructure_;
    CisMatrix<OrdinalType, ScalarType> * matrixForValues_;
    Permutation<OrdinalType, ScalarType> leftPermutation_;
    Permutation<OrdinalType, ScalarType> rightPermutation_;

    bool willKeepStructure_;
    bool willKeepValues_;
    bool isRowOriented_;
    bool haveStructure_;
    bool haveValues_;
    bool hasUnitDiagonal_;
    bool hasClassicHbStructure_;
    bool hasClassicHbValues_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numRC_;
    OrdinalType numEntries_;

    ScalarType ** values_;

    OrdinalType ** indices_;
    OrdinalType * profile_;
    ScalarType * allValues_;
    OrdinalType * allIndices_;
    double costOfMatVec_;

  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseMultiply<OrdinalType, ScalarType>::BaseSparseMultiply() 
    : SparseOperation<OrdinalType, ScalarType>(),
      matrixForStructure_(0),
      matrixForValues_(0),
      willKeepStructure_(false),
      willKeepValues_(false),
      isRowOriented_(true),
      haveStructure_(false),
      haveValues_(false),
      hasUnitDiagonal_(false),
      hasClassicHbStructure_(false),
      hasClassicHbValues_(false),
      numRows_(0),
      numCols_(0),
      numRC_(0),
      numEntries_(0),
      values_(0),
      indices_(0),
      profile_(0),
      allValues_(0),
      allIndices_(0),
      costOfMatVec_(0.0) {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseMultiply<OrdinalType, ScalarType>::BaseSparseMultiply(const BaseSparseMultiply<OrdinalType,
								  ScalarType> &source) 
    : SparseOperation<OrdinalType, ScalarType>(source),
      matrixForStructure_(source.matrixForStructure_),
      matrixForValues_(source.matrixForValues_),
      leftPermutation_(source.leftPermutation_),
      rightPermutation_(source.rightPermutation_),
      willKeepStructure_(source.willKeepStructure_),
      willKeepValues_(source.willKeepValues_),
      isRowOriented_(source.isRowOriented_),
      haveStructure_(source.haveStructure_),
      haveValues_(source.haveValues_),
      hasUnitDiagonal_(source.hasUnitDiagonal_),
      hasClassicHbStructure_(source.ClassicHbStructure_),
      hasClassicHbValues_(sourceClassicHbValues._),
      numRows_(source.numRows_),
      numCols_(source.numCols_),
      numRC_(source.numRC_),
      numEntries_(source.numEntries_),
      values_(source.values_),
      indices_(source.indices_),
      profile_(source.profile_),
      allValues_(source.allValues_),
      allIndices_(source.allIndices_),
      costOfMatVec_(source.costOfMatVec_) {

    copyProfile();
    copyStructure();
    copyValues();
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyProfile() {

    if (profile_!=0) {
      OrdinalType * old_profiles = profiles_;
      profiles_ = new OrdinalType*[NumRC_];
      copyOrdinals(numRC_, old_profiles, profiles_);
    }
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyStructure() {

    if (indices_!=0) {
      OrdinalType ** tmp_indices =indices_;
      indices_ = new OrdinalType*[numRC_];
      if (willKeepStructure_) {
	for (OrdinalType i=0; i< numRC_; i++) indices_[i] = tmp_indices[i];
      }
      else {
	allIndices_ = new OrdinalType[numEntries_]; // Allocate storage for all entries at once
	OrdinalType offset = 0;
	for (OrdinalType i=0; i< numRC_; i++) {
	  indices_[i] = allIndices_+offset;
	  copyOrdinals(profiles_[i], tmp_indices_[i], indices_[i]);
	  offset += profiles_[i];
	}
	hasClassicHbStructure_ = true;
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyValues() {

    if (values_!=0) {
      ScalarType ** tmp_values =values_;
      values_ = new ScalarType*[numRC_];
      if (willKeepValues_) {
	for (OrdinalType i=0; i< numRC_; i++) values_[i] = tmp_values[i];
      }
      else {
	allValues_ = new ScalarType[numEntries_]; // Allocate storage for all entries at once
	OrdinalType offset = 0;
	for (OrdinalType i=0; i< numRC_; i++) {
	  values_[i] = allValues_+offset;
	  copyScalars(profiles_[i], tmp_values_[i], values_[i]);
	  offset += profiles_[i];
	}
	hasClassicHbValues_ = true;
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::deleteStructureAndProfile() {


    // If profile is present, then delete it
    // If indices are present and we allocated the storage, then delete it

    if (indices_!=0) {
      if (!willKeepStructure_ && allIndices_!=0) delete [] allIndices_;
      delete [] indices_;
      }
    if (profile_!=0) delete [] profile_;
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::deleteValues() {

    // If values are present and we allocated the storage, then delete it

    if (values_!=0) {
      if (!willKeepValues_ && allValues_!=0) delete [] allValues_;
      delete [] values_;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyOrdinals(OrdinalType len, 
								 OrdinalType * vecIn, 
								 OrdinalType * vecOut) {
    for (OrdinalType i=0; i<len; i++) vecOut[i] = vecIn[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseMultiply<OrdinalType, ScalarType>::copyScalars(OrdinalType len, 
								ScalarType * vecIn,
								ScalarType * vecOut) {
    for (OrdinalType i=0; i<len; i++) vecOut[i] = vecIn[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseMultiply<OrdinalType, ScalarType>::~BaseSparseMultiply(){

    deleteValues();
    deleteStructureAndProfile();

  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseMultiply<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
								       bool willKeepStructure) {


    if (haveStructure_) return(-1); // Can only call this one time!
    matrixForStructure_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i;
    willKeepStructure_ = willKeepStructure;
    isRowOriented_ = A.getIsRowOriented();
    numRows_ = A.getNumRows();
    numCols_ = A.getNumCols();
    numEntries_ = A.getNumEntries();
    hasUnitDiagonal_ = A.getHasImplicitUnitDiagonal();
    numRC_ = numCols_;
    if (isRowOriented_) numRC_ = numRows_;

    profile_ = new OrdinalType[numRC_];
    indices_ = new OrdinalType*[numRC_];
    OrdinalType numRCEntries;
    OrdinalType * indicesRC;

    if (willKeepStructure) {
      hasClassicHbStructure_ = true; // Assume packed storage, turn it off below if not
      for (i=0; i<numRC_; i++) {
	int ierr = A.getIndices(i, numRCEntries, indicesRC);
	if (ierr<0) return(ierr);
	profile_[i] = numRCEntries;
	indices_[i] = indicesRC;
	if (i>0) 
	  if (indicesRC - indices_[i-1] !=profile_[i-1]) hasClassicHbStructure_ = false;
      }
    }
    else { // If user will not keep structure, we must copy it
      
      allIndices_ = new OrdinalType[numEntries_]; // Allocate storage for all entries at once
      
      OrdinalType offset = 0;
      for (i=0; i< numRC_; i++) {
	int ierr = A.getIndices(i, numRCEntries, indicesRC);
	if (ierr<0) return(ierr);
	profile_[i] = numRCEntries;
	indices_[i] = allIndices_+offset;
	copyOrdinals(numRCEntries, indicesRC, indices_[i]);
	offset += numRCEntries;
      }
      hasClassicHbStructure_ = true; // True by construction
    }

    costOfMatVec_ = 2.0 * ((double) numEntries_);
    if (hasUnitDiagonal_) costOfMatVec_ += 2.0 * ((double) numRC_);
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

    if (willKeepValues_) {
      hasClassicHbValues_ = true; // Assume packed storage, turn it off below if not
      for (i=0; i<numRC_; i++) {
	int ierr = A.getValues(i, valuesRC);
	if (ierr<0) return(ierr);
	values_[i] = valuesRC;
	if (i>0) 
	  if (valuesRC - values_[i-1] !=profile_[i-1]) hasClassicHbValues_ = false;
      }
    }
    else { // If user will not keep values, we must copy it
      if (allValues_==0) allValues_ = new ScalarType[numEntries_]; // Allocate storage for all entries at once
      OrdinalType offset = 0;
      for (i=0; i< numRC_; i++) {
	int ierr = A.getValues(i, valuesRC);
	if (ierr<0) return(ierr);
	values_[i] = allValues_+offset;
	copyScalars(profile_[i], valuesRC, values_[i]);
	offset += profile_[i];
      }
      hasClassicHbValues_ = true; // True by construction
    }
    haveValues_ = true;

    if (checkStructure) { // Compare strucuture of current matrix to structure of matrix used in initializeStructure
      if (matrixForValues_==matrixForStructure_) return(0); // Pointing to the same matrix, assume structure is same
      OrdinalType numRCEntries;
      OrdinalType * indicesRC_ref;
      OrdinalType * indicesRC;

      for (i=0; i<numRC_; i++) {
	int ierr = matrixForValues_->getIndices(i, numRCEntries, indicesRC);
	if (ierr<0) return(-1);
	if (numRCEntries!=profile_[i]) return(-1);
	indicesRC_ref = indices_[i];
	for (j=0; j<numRCEntries; j++) if (indicesRC[j]!=indicesRC_ref[j]) return(-1);
	
      }
    }
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
	if (hasUnitDiagonal_) 
	  sum = xp[i];
	for(j = 0; j < curNumEntries; j++)
	  sum += curValues[j] * xp[curIndices[j]];
	yp[i] = sum;
      }
    }
    else {
      
      if (hasUnitDiagonal_) 
	for(i = 0; i < numRC_; i++)
	  yp[i] = xp[i]; // Initialize y for transpose multiply
      else
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
    updateFlops(costOfMatVec_);
    return(0);
  }


  //==============================================================================
  template<>
  int BaseSparseMultiply<int, double>::apply(const Vector<int, double>& x, 
						    Vector<int, double> & y,
						    bool transA, bool conjA) const {

    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA) return(-2); // Unsupported at this time
    if (x.getLength()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getLength()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x

    int i, j, curNumEntries;
    int * curIndices;
    double * curValues;

    int * profile = profile_;
    int ** indices = indices_;
    double ** values = values_;

    double * xp = x.getValues();
    double * yp = y.getValues();
    
    // Special case of classic HB matrix
    if (hasClassicHbStructure_ && hasClassicHbValues_) {
      int itrans = 1;
      if ((isRowOriented_ && !transA) ||
	  (!isRowOriented_ && transA)) itrans = 0;
      int udiag = 0;
      if (hasUnitDiagonal_) udiag = 1;
      KOKKOS_DCRSMV_F77(&itrans, &udiag, &numRows_, &numCols_, *values, *indices, profile, xp, yp);
    }
    else {// General case
      
      if ((isRowOriented_ && !transA) ||
	  (!isRowOriented_ && transA)) {
	
	for(i = 0; i < numRC_; i++) {
	  curNumEntries = *profile++;
	  curIndices = *indices++;
	  curValues  = *values++;
	  double sum = 0.0;
	  if (hasUnitDiagonal_) 
	    sum = xp[i];
	  for(j = 0; j < curNumEntries; j++)
	    sum += curValues[j] * xp[curIndices[j]];
	  yp[i] = sum;
	}
      }
      else {
	
	if (hasUnitDiagonal_) 
	  for(i = 0; i < numRC_; i++)
	    yp[i] = xp[i]; // Initialize y for transpose multiply
	else
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
    }
    updateFlops(costOfMatVec_);
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

    if ((isRowOriented_ && !transA) ||
	(!isRowOriented_ && transA)) {

      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	curIndices = *indices++;
	curValues  = *values++;
	for (k=0; k<numVectors; k++) {
	  ScalarType sum = 0.0;
	  if (hasUnitDiagonal_) 
	    sum = xp[k][i];
	  for(j = 0; j < curNumEntries; j++)
	    sum += curValues[j] * xp[k][curIndices[j]];
	  yp[k][i] = sum;
	}
      }
    }
    else {
      
      if (hasUnitDiagonal_) 
	for (k=0; k<numVectors; k++)
	  for(i = 0; i < numRC_; i++)
	    yp[k][i] = xp[k][i]; // Initialize y for transpose multiply
      else
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
    updateFlops(costOfMatVec_ * ((double) numVectors));
    return(0);
  }
  //==============================================================================
  template<>
  int BaseSparseMultiply<int, double>::apply(const MultiVector<int, double>& x, 
						    MultiVector<int, double> & y,
						    bool transA, bool conjA) const {
    if (!haveValues_) return(-1); // Can't compute without values!
    if (x.getNumRows()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getNumRows()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x
    int numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-5); // Not the same number of vectors in x and y

    int i, j, k, curNumEntries;
    int * curIndices;
    double * curValues;

    int * profile = profile_;
    int ** indices = indices_;
    double ** values = values_;

    double ** xp = x.getValues();
    double ** yp = y.getValues();

    // Special case of classic HB matrix with Fortran compatible multivectors
    if (hasClassicHbStructure_ && hasClassicHbValues_ && x.getIsStrided() && y.getIsStrided()
	&& x.getColInc()==1 && y.getColInc()==1) {
      int itrans = 1;
      if ((isRowOriented_ && !transA) ||
	  (!isRowOriented_ && transA)) itrans = 0;
      int udiag = 0;
      if (hasUnitDiagonal_) udiag = 1;
      int ldx = x.getRowInc();
      int ldy = y.getRowInc();
      KOKKOS_DCRSMM_F77( &itrans, &udiag, &numRows_, &numCols_, *values, *indices, profile, 
			 *xp, &ldx, *yp,&ldy , &numVectors );
    }
    else { // General case
      if ((isRowOriented_ && !transA) ||
	  (!isRowOriented_ && transA)) {
	
	for(i = 0; i < numRC_; i++) {
	  curNumEntries = *profile++;
	  curIndices = *indices++;
	  curValues  = *values++;
	  for (k=0; k<numVectors; k++) {
	    double sum = 0.0;
	    if (hasUnitDiagonal_) 
	      sum = xp[k][i];
	    for(j = 0; j < curNumEntries; j++)
	      sum += curValues[j] * xp[k][curIndices[j]];
	    yp[k][i] = sum;
	  }
	}
      }
      else {
	
	if (hasUnitDiagonal_) 
	  for (k=0; k<numVectors; k++)
	    for(i = 0; i < numRC_; i++)
	      yp[k][i] = xp[k][i]; // Initialize y for transpose multiply
	else
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
    }
    updateFlops(costOfMatVec_ * ((double) numVectors));
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_BASESPARSEMULTIPLY_H */
