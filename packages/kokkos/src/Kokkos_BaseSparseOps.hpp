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

#ifndef KOKKOS_BASESPARSEOPS_H
#define KOKKOS_BASESPARSEOPS_H

#include "Kokkos_CompObject.hpp" 


namespace Kokkos {

//! Kokkos::BaseSparseOps: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Kokkos::BaseSparseOps enables the piecewise construction and use of real-valued double-precision sparse matrices
  where matrix entries are intended for row access.

  At this time, the primary function provided by Kokkos::BaseSparseOps is matrix times vector and matrix 
  times multi-vector multiplication.  It is also possible to extract matrix rows from a constructed matrix.

  <b>Constructing Kokkos::BaseSparseOps objects</b>

  Constructing Kokkos::BaseSparseOps objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::BaseSparseOps instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
  </ol>

  Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
  create new entries.

  <b> Counting Floating Point Operations </b>

  Each Kokkos::BaseSparseOps object keeps track of the number
  of \e serial floating point operations performed using the specified object as the \e this argument
  to the function.  The Flops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate parallel performance
  numbers.  The ResetFlops() function resets the floating point counter.


*/    

  template<typename OrdinalType, typename ScalarType>
  class BaseSparseOps: public CompObject {
  public:

    //@{ \name Constructors/Destructor.
    //! BaseSparseOps constuctor with variable number of indices per row.
    BaseSparseOps(void);
  
    //! Copy constructor.
    BaseSparseOps(const BaseSparseOps& operator);
	
    //! BaseSparseOps Destructor
    virtual ~BaseSparseOps();
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
	
    //! Returns the result of a Kokkos_BaseSparseOps multiplied by a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to multiply by.
      \param y (Out) A Kokkos::Vector containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;

    //! Returns the result of a Kokkos_BaseSparseOps multiplied by multiple vectors in x, results in y.
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

    //! Returns true if this implementation of Kokkos::BaseSparseOps can benefit from the user keeping the passed in structure.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        analysis of structure.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepStructure parameter 
	passed in to the initializeStructure() method.
    */
    virtual bool getCanUseStructure() const {return(true);} const;

    //! Returns true if this implementation of Kokkos::BaseSparseOps can benefit from the user keeping the passed in values.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        copying of values.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepValues parameter 
	passed in to the initializeValues() method.
    */
    virtual bool getCanUseValues() const {return(true);} const;

    virtual CisMatrix<Ordinaltype, ScalarType> * getMatrix() const;
	
    //@}
  
  protected:
    bool allocated() const {return(allocated_);};
    int setAllocated(bool flag) {allocated_ = flag; return(0);};
	
    void initializeDefaults();
    int allocate();

    bool allocated_;
    bool willKeepStructure_;
    bool willKeepValues_;
    bool isRowOriented_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numEntries_;

    ScalarType ** values_;
    ScalarType * allValues_;

    OrdinalType ** indices_;
    OrdinalType * allIndices_;

    OrdinalType * pntr_;
    OrdinalType * profile_;

  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseOps<OrdinalType, ScalarType>::BaseSparseOps() 
    : CompObject(),
      allocated_(false),
      willKeepStructure_(false),
      willKeepValues_(false),
      isRowOriented_(true),
      numRows_(0),
      numCols_(0),
      numEntries_(0),
      values_(0),
      allValues_(0),
      indices_(0),
      allIndices_(0),
      pntr_(0),
      profile_(0) {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseOps<OrdinalType, ScalarType>::BaseSparseOps(const BaseSparseOps<OrdinalType, ScalarType> &operator) 
    : CompObject(operator),
      allocated_(operator.allocated_),
      willKeepStructure_(operator.willKeepStructure_),
      willKeepValues_(operator.willKeepValues_),
      isRowOriented_(operator.isRowOriented_),
      numRows_(operator.numRows_),
      numCols_(operator.numCols_),
      numEntries_(operator.numEntries_),
      values_(operator.values_),
      allValues_(operator.allValues_),
      indices_(operator.indices_),
      allIndices_(operator.allIndices_),
      pntr_(operator.pntr_),
      profile_(operator.profile_) {

    if (!willKeepStructure_) copyStructure();
    if (!willKeepValues_) copyValues()
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseOps<OrdinalType, ScalarType>::copyStructure() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // pntr_, profile_, allIndices_, indices_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;
    OrdinalType numRC = numCols_;
    if (isRowOriented_) numRC = numRows_;

    if (pntr_!=0) copyOrdinals(pntr_, numRC+1);
    if (profile_!=0) copyOrdinals(profiles_, numRC);
    if (allIndices_!=0) copyOrdinals(allIndices_, numEntries_);
    if (indices_!=0) {
      OrdinalType ** tmp_indices =indices_;
      indices_ = new OrdinalType*[numRC];
      for (i=0; i< numRC; i++) {
	indices_[i] = tmp_indices[i];
	copyOrdinals(indices_[i], profiles_[i]);
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseOps<OrdinalType, ScalarType>::copyValues() {

    // Note:  In order to call this function, the following attributes must be set, and pointers
    //        must be set to the data that is being copied:
    // numCols_, numRows_, isRowOriented_ must be set
    // allValues_, values_ must be set to data that will be copied.
    // These pointers will be set to new memory that is a copy of these arrays (except 
    // for pointers that are null.

    OrdinalType i;
    OrdinalType numRC = numCols_;
    if (isRowOriented_) numRC = numRows_;

    if (allValues_!=0) copyScalars(allValues_, numEntries_);
    if (values_!=0) {
      ScalarType ** tmp_values =values_;
      values_ = new ScalarType*[numRC];
      for (i=0; i< numRC; i++) {
	values_[i] = tmp_values[i];
	copyScalars(values_[i], profiles_[i]);
      }
    }
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseOps<OrdinalType, ScalarType>::copyOrdinals(OrdinalType * vec, OrdinalType len) {

    OrdinalType i;
    OrdinalType * tmp_pntr = vec;
    vec = new OrdinalType[len];
    for (i=0; i<len; i++) vec[i] = tmp_pntr[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseOps<OrdinalType, ScalarType>::copyScalars(ScalarType * vec, OrdinalType len) {

    OrdinalType i;
    ScalarType * tmp_vals = vec;
    vec = new OrdinalType[len];
    for (i=0; i<len; i++) vec[i] = tmp_vals[i];
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  BaseSparseOps<OrdinalType, ScalarType>::~BaseSparseOps(){}

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void BaseSparseOps<OrdinalType, ScalarType>::initializeDefaults() {
    return;
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseOps<OrdinalType, ScalarType>::allocate() {
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseOps<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
								  bool willKeepStructure) {

    OrdinalType i;
    willKeepStructure_ = willKeepStructure;
    isRowOriented_ = A.isRowOriented();
    numRows_ = A.numRows();
    numCols_ = A.numCols();
    numEntries_ = A.numEntries();
    OrdinalType numRC = numCols_;
    if (isRowOriented_) numRC = numRows_;

    OrdinalType * profiles = new[numRC+1];  // Make this longer in case we can use it for pntr data
    indices_ = new[numRC];

    OrdinalType numRCEntries;
    OrdinalType * indicesRC;

    bool contiguousStructure = true; // Assume structure is contiguously stored.  Check in the loop below.
    A.getIndices(0, numRCEntries, indicesRC);
    profiles[0] = numRCEntries;
    indices_[0] = indicesRC;
    i = 1;
    while (contiguousStructure && i<numRC) {
      A.getIndices(i, numRCEntries, indicesRC);
      profiles[i] = numRCEntries;
      indices_[i] = indicesRC;
      if (indices_[i-1]+profiles_[i]!=indices_[i]) contiguousStructure = false; // Not contiguous
      i++;
    }

    if (contiguousStructure) {
      allIndices_ = indices_[0];
      delete [] indices_;
      indices_ = 0;
    }

    copyStructure();

    if (contiguousStructure) {
      pntr_ = profile;
      OrdinalType numEntries = profile[0];
      pntr_[0] = 0;
      for (i=1; i< numRC+1; i++) {
	numEntries = profile[i];
	pntr_[i] = pntr_[i-1] + numEntries;
      }
    }


    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseOps<OrdinalType, ScalarType>::initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, 
							       bool willKeepValues, bool checkStructure) {
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseOps<OrdinalType, ScalarType>::apply(const Vector<OrdinalType, ScalarType>& x, 
						    Vector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int BaseSparseOps<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
						    MultiVector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_BASESPARSEOPS_H */
