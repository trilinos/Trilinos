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

#ifndef KOKKOS_BASESPARSESOLVE_H
#define KOKKOS_BASESPARSESOLVE_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_SparseOperation.hpp" 


namespace Kokkos {

//! Kokkos::BaseSparseSolve: A reference class for computing sparse matrix triangular solve operations.

/*! The Kokkos::BaseSparseSolve provide basic functionality for computing sparse triangular solves with one or more
    right-hand-side vectors.  This class is templated on the ordinal (integer) and 
    scalar (floating point) types, so it can compute using any reasonable data type.  It 
    implements the Kokkos::SparseOperation base class.

  <b>Constructing Kokkos::BaseSparseSolve objects</b>

  Constructing Kokkos::BaseSparseSolve objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::BaseSparseSolve instance:  The constructor takes no arguments.
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

  Each Kokkos::BaseSparseSolve object keeps track of the number
  of floating point operations performed using the specified object as the \e this argument
  to the function.  The getFlops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate  performance
  numbers.  The resetFlops() function resets the floating point counter.


*/    

  template<typename OrdinalType, typename ScalarType>
  class BaseSparseSolve: public virtual SparseOperation<OrdinalType, ScalarType> {
  public:

    //@{ \name Constructors/Destructor.
    //! BaseSparseSolve constuctor with variable number of indices per row.
    BaseSparseSolve();
  
    //! Copy constructor.
    BaseSparseSolve(const BaseSparseSolve& source);
	
    //! BaseSparseSolve Destructor
    virtual ~BaseSparseSolve();
    //@}
    //@{ \name Abstract Kokkos::CisMatrix Interface Initialization Methods
 
    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
             about the matrix can be obtained via this interface, including whether or not the matrix is upper or lower
	     triangular, and whether or not the diagonal is part of the structure, or should be implicitly assume to be
	     unit diagonal.
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
             about the matrix can be obtained via this interface, including whether or not the matrix is upper or lower
	     triangular, and whether or not the diagonal is part of the structure, or should be implicitly assume to be
	     unit diagonal.
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
	
    //! Returns the result of a Kokkos_BaseSparseSolve multiplied by a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to solve with.
      \param y (Out) A Kokkos::Vector containing results.  Note that any implementation must support x and y being 
             the same object.
      \param transA (In) If true, solve using the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, solve using the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;

    //! Returns the result of a Kokkos_BaseSparseSolve multiplied by multiple vectors in x, results in y.
    /*! 
      \param x (In) A Kokkos::MultiVector to solve with.
      \param y (Out) A Kokkos::MultiVector containing results.  Note that any implementation must support x and y being 
             the same object.
      \param transA (In) If true, solve using the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, solve using the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;
    //@}
	
    //@{ \name Operator attribute access methods.

    //! Returns true if this implementation of Kokkos::BaseSparseSolve can benefit from the user keeping the passed in structure.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        analysis of structure.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepStructure parameter 
	passed in to the initializeStructure() method.
    */
    virtual bool getCanUseStructure() const {return(true);};

    //! Returns true if this implementation of Kokkos::BaseSparseSolve can benefit from the user keeping the passed in values.
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
    bool isUpper_;
    bool hasUnitDiagonal_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numRC_;
    OrdinalType numEntries_;

    ScalarType ** values_;

    OrdinalType ** indices_;
    OrdinalType * profile_;
    ScalarType * allValues_;
    OrdinalType * allIndices_;
    double costOfSolve_;

  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseSolve<OrdinalType, ScalarType>::BaseSparseSolve() 
    : SparseOperation<OrdinalType, ScalarType>(),
      matrixForStructure_(0),
      matrixForValues_(0),
      willKeepStructure_(false),
      willKeepValues_(false),
      isRowOriented_(true),
      haveStructure_(false),
      haveValues_(false),
      isUpper_(false),
      hasUnitDiagonal_(false),
      numRows_(0),
      numCols_(0),
      numRC_(0),
      numEntries_(0),
      values_(0),
      indices_(0),
      profile_(0),
      allValues_(0),
      allIndices_(0),
      costOfSolve_(0.0) {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseSolve<OrdinalType, ScalarType>::BaseSparseSolve(const BaseSparseSolve<OrdinalType,
							    ScalarType> &source) 
    : SparseOperation<OrdinalType, ScalarType>(source),
      matrixForStructure_(source.matrixForStructure_),
      matrixForValues_(source.matrixForValues_),
      willKeepStructure_(source.willKeepStructure_),
      willKeepValues_(source.willKeepValues_),
      isRowOriented_(source.isRowOriented_),
      haveStructure_(source.haveStructure_),
      haveValues_(source.haveValues_),
      isUpper_(source.isUpper_),
      hasUnitDiagonal_(source.hasUnitDiagonal_),
      numRows_(source.numRows_),
      numCols_(source.numCols_),
      numRC_(source.numRC_),
      numEntries_(source.numEntries_),
      values_(source.values_),
      indices_(source.indices_),
      profile_(source.profile_),
      allValues_(source.allValues_),
      allIndices_(source.allIndices_),
      costOfSolve_(source.costOfSolve_) {

    copyProfile();
    copyStructure();
    copyValues();
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseSolve<OrdinalType, ScalarType>::copyProfile() {

    if (profile_!=0) {
      OrdinalType * old_profiles = profiles_;
      profiles_ = new OrdinalType*[NumRC_];
      copyOrdinals(numRC_, old_profiles, profiles_);
    }
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseSolve<OrdinalType, ScalarType>::copyStructure() {

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
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseSolve<OrdinalType, ScalarType>::copyValues() {

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
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseSolve<OrdinalType, ScalarType>::deleteStructureAndProfile() {


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
  void BaseSparseSolve<OrdinalType, ScalarType>::deleteValues() {

    // If values are present and we allocated the storage, then delete it

    if (values_!=0) {
      if (!willKeepValues_ && allValues_!=0) delete [] allValues_;
      delete [] values_;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseSolve<OrdinalType, ScalarType>::copyOrdinals(OrdinalType len, 
								 OrdinalType * vecIn, 
								 OrdinalType * vecOut) {
    for (OrdinalType i=0; i<len; i++) vecOut[i] = vecIn[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseSolve<OrdinalType, ScalarType>::copyScalars(OrdinalType len, 
								ScalarType * vecIn,
								ScalarType * vecOut) {
    for (OrdinalType i=0; i<len; i++) vecOut[i] = vecIn[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseSolve<OrdinalType, ScalarType>::~BaseSparseSolve(){

    deleteValues();
    deleteStructureAndProfile();

  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseSolve<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
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

    isUpper_ = A.getIsUpperTriangular();
    hasUnitDiagonal_ = A.getHasImplicitUnitDiagonal();
    if (isRowOriented_) numRC_ = numRows_;

    profile_ = new OrdinalType[numRC_];
    indices_ = new OrdinalType*[numRC_];
    OrdinalType numRCEntries;
    OrdinalType * indicesRC;

    if (willKeepStructure) {
      for (i=0; i<numRC_; i++) {
	int ierr = A.getIndices(i, numRCEntries, indicesRC);
	if (ierr<0) return(ierr);
	profile_[i] = numRCEntries;
	indices_[i] = indicesRC;
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
    }

    costOfSolve_ = 2.0 * ((double) numEntries_);
    if (hasUnitDiagonal_) costOfSolve_ += 2.0 * ((double) numRC_);
    haveStructure_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseSolve<OrdinalType, ScalarType>::initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, 
							       bool willKeepValues, bool checkStructure) {

    if (!haveStructure_) return(-1); // Must have structure first!

    matrixForValues_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    willKeepValues_ = willKeepValues;

    
    values_ = new ScalarType*[numRC_];

    ScalarType * valuesRC;

    if (willKeepValues_) {
      for (i=0; i<numRC_; i++) {
	int ierr = A.getValues(i, valuesRC);
	if (ierr<0) return(ierr);
	values_[i] = valuesRC;
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
  int BaseSparseSolve<OrdinalType, ScalarType>::apply(const Vector<OrdinalType, ScalarType>& x, 
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

      if ((!transA && isUpper_) || transA && !isUpper_) {
	profile += numRC_-1; // Point to end of structures
	indices += numRC_-1;
	values += numRC_-1;

	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--;
 	for(i = numRC_-1; i >=0; i--) {
	  curNumEntries = *profile--;
	  curIndices = *indices--;
	  curValues  = *values--;
	  ScalarType sum = 0.0;
	  for(j = j0; j < curNumEntries; j++)
	    sum += curValues[j] * yp[curIndices[j]];
	  if (hasUnitDiagonal_)
	    yp[i] = xp[i] - sum;
	  else
	    yp[i] = (xp[i] - sum)/curValues[0];
	}
      }
      else { // Lower triangular
	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--;
	for(i = 0; i < numRC_; i++) {
	  curNumEntries = *profile++ - j0;
	  curIndices = *indices++;
	  curValues  = *values++;
	  ScalarType sum = 0.0;
	  for(j = 0; j < curNumEntries; j++)
	    sum += curValues[j] * yp[curIndices[j]];
	  if (hasUnitDiagonal_)
	    yp[i] = xp[i] - sum;
	  else
	    yp[i] = (xp[i] - sum)/curValues[curNumEntries];
	}
      }
    }
    else { // ColOriented and no tranpose or RowOriented and transpose
      
      for(i = 0; i < numRC_; i++)
	yp[i] = xp[i]; // Initialize y for transpose multiply

      if ((!transA && !isUpper_) || transA && isUpper_) {
	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--;
	for(i = 0; i < numRC_; i++) {
	  curNumEntries = *profile++;
	  curIndices = *indices++;
	  curValues  = *values++;
	  if (!hasUnitDiagonal_) 
	    yp[i] = yp[i]/curValues[0];
	  for(j = j0; j < curNumEntries; j++)
	    yp[curIndices[j]] -= curValues[j] * yp[i];
	}
      }
      else { // Lower triangular

	profile += numRC_-1; // Point to end of structures
	indices += numRC_-1;
	values += numRC_-1;

	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--; // Include first term if no diagonal

	for(i = numRC_-1; i>=0; i--) {
	  curNumEntries = *profile-- - j0;
	  curIndices = *indices--;
	  curValues  = *values--;
	    if (!hasUnitDiagonal_) 
	      yp[i] = yp[i]/curValues[curNumEntries];
	    for(j = 0; j < curNumEntries; j++)
	      yp[curIndices[j]] -= curValues[j] * yp[i];
	}
      }
    }
    updateFlops(costOfSolve_);
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseSolve<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
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

      if ((!transA && isUpper_) || transA && !isUpper_) {
	profile += numRC_-1; // Point to end of structures
	indices += numRC_-1;
	values += numRC_-1;

	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--;
	for(i = numRC_-1; i >=0; i--) {
	  curNumEntries = *profile--;
	  curIndices = *indices--;
	  curValues  = *values--;
	  for (k=0; k<numVectors; k++) {
	    ScalarType sum = 0.0;
	    for(j = j0; j < curNumEntries; j++)
	      sum += curValues[j] * yp[k][curIndices[j]];
	    if (hasUnitDiagonal_)
	      yp[k][i] = xp[k][i] - sum;
	    else
	      yp[k][i] = (xp[k][i] - sum)/curValues[0];
	  }
	}
      }
      else { // Lower triangular
	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--;
	for(i = 0; i < numRC_; i++) {
	  curNumEntries = *profile++ - j0;
	  curIndices = *indices++;
	  curValues  = *values++;
	  for (k=0; k<numVectors; k++) {
	    ScalarType sum = 0.0;
	    for(j = 0; j < curNumEntries; j++)
	      sum += curValues[j] * yp[k][curIndices[j]];
	    if (hasUnitDiagonal_)
	      yp[k][i] = xp[k][i] - sum;
	    else
	      yp[k][i] = (xp[k][i] - sum)/curValues[curNumEntries];
	  }
	}
      }
    }
    else { // ColOriented and no tranpose or RowOriented and transpose
      
      for (k=0; k<numVectors; k++)
	for(i = 0; i < numRC_; i++)
	  yp[k][i] = xp[k][i]; // Initialize y
      
      if ((!transA && !isUpper_) || transA && isUpper_) {
	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--;
	for(i = 0; i < numRC_; i++) {
	  curNumEntries = *profile++;
	  curIndices = *indices++;
	  curValues  = *values++;
	  for (k=0; k<numVectors; k++) {
	    if (!hasUnitDiagonal_) 
	      yp[k][i] = yp[k][i]/curValues[0];
	    for(j = j0; j < curNumEntries; j++)
	      yp[k][curIndices[j]] -= curValues[j] * yp[k][i];
	  }
	}
      }
      else { // Lower triangular

	profile += numRC_-1; // Point to end of structures
	indices += numRC_-1;
	values += numRC_-1;

	OrdinalType j0 = 1;
	if (hasUnitDiagonal_) j0--; // Include first term if no diagonal

	for(i = numRC_-1; i>=0; i--) {
	  curNumEntries = *profile-- - j0;
	  curIndices = *indices--;
	  curValues  = *values--;
	  for (k=0; k<numVectors; k++) {
	    if (!hasUnitDiagonal_) 
	      yp[k][i] = yp[k][i]/curValues[curNumEntries];
	    for(j = 0; j < curNumEntries; j++)
	      yp[k][curIndices[j]] -= curValues[j] * yp[k][i];
	  }
	}
      }
    }
    updateFlops(costOfSolve_ * ((double) numVectors));
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_BASESPARSESOLVE_H */
