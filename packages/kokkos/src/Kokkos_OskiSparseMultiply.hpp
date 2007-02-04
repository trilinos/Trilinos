//@HEADER
// ************************************************************************
// 
//                 Kokkos: A Fast Kernel Package
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

#ifndef KOKKOS_OSKISPARSEMULTIPLY_H
#define KOKKOS_OSKISPARSEMULTIPLY_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_OskiMatrix.hpp" 
#include "Kokkos_Vector.hpp"
#include "Kokkos_Permutation.hpp"
#include "Kokkos_SparseOperation.hpp" 
extern "C" {
  #include <oski/oski.h>
}

namespace Kokkos {

// ** Currently we are using CisMatrix and MultiVector in places where we 
// actually require OSKI objects.  We need to enforce that the types are right!

//! Kokkos::OskiSparseMultiply: A class for computing sparse matrix multiplication operations using functions provided by OSKI.

/*! The Kokkos::OskiOskiSparseMultiply provide basic functionality for computing sparse matrix times vector, or
    sparse matrix times multivector operations.  This class is templated on the ordinal (integer) and scalar (floating
    point) types, so it can compute using any reasonable data type.  Support
    for those data types must have been compiled into OSKI.

  <b>Constructing Kokkos::OskiSparseMultiply objects</b>

  Constructing Kokkos::OskiSparseMultiply objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::OskiSparseMultiply instance:  The constructor takes no arguments.
  <li> Register the structure of a Kokkos::OskiMatrix object using initializeStructure(): 
       We provide this method so that derived implementations can
       take advantage of multiple problems that have the same structure.  In this situation, initializeStructure() would
       be called once and then initializeValues() would be called repeatedly, amortizing the cost of setting up the structure.
       This method may be called only once.
  <li> Register the values of a Kokkos::OskiMatrix object using initializeValues(): This method is used to pass values to the
       multiply class.  It can be called repeatedly if multiple matrices have the same structure.
  </ol>

  <b> Counting Floating Point Operations </b>

  Each Kokkos::OskiSparseMultiply object keeps track of the number
  of floating point operations performed using the specified object as the \e this argument
  to the function.  The getFlops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate  performance
  numbers.  The resetFlops() function resets the floating point counter.

*/    

  template<typename OrdinalType, typename ScalarType>
  class OskiSparseMultiply: public virtual SparseOperation<OrdinalType, ScalarType> {
  public:

    //! @name Constructors/Destructor

    //@{

    //! OskiSparseMultiply constuctor with variable number of indices per row.
    OskiSparseMultiply();
  
    //! Copy constructor.
    OskiSparseMultiply(const OskiSparseMultiply& source);
	
    //! OskiSparseMultiply Destructor
    virtual ~OskiSparseMultiply();
    //@}
    //! @name Abstract Kokkos::OskiMatrix Interface Initialization Methods

    //@{ 

    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the Kokkos::OskiMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::OskiMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepStructure (In) This argument is unused by this implementation of the SparseOperation
             class since structure and values will be shared.
      \return Integer error code, set to 0 if successful.
    */
    virtual int initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepStructure = true);
 
    //! Initialize values of matrix
    /*!
      This interface supports matrices that implement the Kokkos::OskiMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::OskiMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepValues (In) This argument is unused by this implementation of the SparseOperation 
             class since structure and values will be shared.
      \param checkStructure (In) If set to true, the structure of A will be checked against the structure of
      the matrix passed in to the initializeStructure() methods.  This parameter is false by default.

      \return Integer error code, set to 0 if successful, returns - 1 if checkStructure is true and structure is changed.
    */
    virtual int initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepValues = true,
				 bool checkStructure = false);
 
    //@}

    //! @name Computational methods

    //@{

// We won't use a Vector version of this because OskiVector inherits from 
// OskiMultiVector - we need this to satisfy the requirement that certain
// methods be implemented.

    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y,
                      bool transA = false, bool conjA = false) const;
	
    //! Returns the result of a Kokkos_OskiSparseMultiply multiplied by multiple vectors in x, results in y.
    /*! 
      \param x (In) A Kokkos::OskiMultiVector to multiply by.
      \param y (Out) A Kokkos::OskiMultiVector containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;
    //@}
	
    //! @name Operator attribute access methods

    //@{

    //! Returns true for this implementation.
    /*! This implementation will use the user's copy of the matrix structure.
    */
    virtual bool getCanUseStructure() const {return(true);};

    //! Returns true for this implementation.
    /*! This implementation will use the user's copy of the matrix values.
    */
    virtual bool getCanUseValues() const {return(true);};

    //! Returns a reference to the most recent Kokkos::OskiMatrix that was passed into the \e this object.
    virtual const CisMatrix<OrdinalType, ScalarType> & getMatrix() const {
      if (matrixForValues_==0) return(*matrixForStructure_);
      else return(*matrixForValues_);
    };
		
    //@}

    //! Returns a reference to the left Kokkos::Permutation object, which is the identity for this implementation.
    virtual const Permutation<OrdinalType, ScalarType> & getLeftPermutation() const {
      return(leftPermutation_);};

    //! Returns a reference to the right Kokkos::Permutation object, which is the identity for this implementation.
    virtual const Permutation<OrdinalType, ScalarType> & getRightPermutation() const {
      return(rightPermutation_);};

  protected:

    void copyEntries();
    void deleteStructureAndValues();

    struct EntryStruct {
      OrdinalType index;
      ScalarType value;
    }; 
    typedef struct EntryStruct Entry;
    
    CisMatrix<OrdinalType, ScalarType> * matrixForStructure_;
    CisMatrix<OrdinalType, ScalarType> * matrixForValues_;
    Permutation<OrdinalType, ScalarType> leftPermutation_;
    Permutation<OrdinalType, ScalarType> rightPermutation_;

    bool isRowOriented_; //Some of these can be removed!
    bool haveStructure_;
    bool haveValues_;
    bool hasUnitDiagonal_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numRC_;
    OrdinalType numEntries_;

    OrdinalType * profile_;
    double costOfMatVec_;
    Entry * allEntries_;
  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  OskiSparseMultiply<OrdinalType, ScalarType>::OskiSparseMultiply() 
    : matrixForStructure_(0),
      matrixForValues_(0),
      isRowOriented_(true),
      haveStructure_(false),
      haveValues_(false),
      hasUnitDiagonal_(false),
      numRows_(0),
      numCols_(0),
      numRC_(0),
      numEntries_(0),
      profile_(0),
      costOfMatVec_(0.0),
      allEntries_(0) {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  OskiSparseMultiply<OrdinalType, ScalarType>::OskiSparseMultiply(const OskiSparseMultiply<OrdinalType, ScalarType> &source) 
    : matrixForStructure_(source.matrixForStructure_),
      matrixForValues_(source.matrixForValues_),
      leftPermutation_(source.leftPermutation_),
      rightPermutation_(source.rightPermutation_),
      isRowOriented_(source.isRowOriented_),
      haveStructure_(source.haveStructure_),
      haveValues_(source.haveValues_),
      hasUnitDiagonal_(source.hasUnitDiagonal_),
      numRows_(source.numRows_),
      numCols_(source.numCols_),
      numRC_(source.numRC_),
      numEntries_(source.numEntries_),
      profile_(source.profile_),
      costOfMatVec_(source.costOfMatVec_),
      allEntries_(source.allEntries_) {

    copyEntries();
  }

  //==============================================================================
// **Not sure if I will need this method or not
  template<typename OrdinalType, typename ScalarType>
  void OskiSparseMultiply<OrdinalType, ScalarType>::copyEntries() {

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
  void OskiSparseMultiply<OrdinalType, ScalarType>::deleteStructureAndValues() {


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
  OskiSparseMultiply<OrdinalType, ScalarType>::~OskiSparseMultiply(){     

    deleteStructureAndValues();

  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int OskiSparseMultiply<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
								       bool willKeepStructure) {

// It we just have to pass the Oski Objects, we may not need this.
// It could be as simple as storing a pointer to the appropriate OskiMatrix
    if (haveStructure_) return(-1); // Can only call this one time!

    matrixForStructure_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    isRowOriented_ = A.getIsRowOriented();
    hasUnitDiagonal_ = A.getHasImplicitUnitDiagonal();
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
      if (ierr<0) return(ierr);
      profile_[i] = numRCEntries;
      Entry * curRC = allEntries_+offset;
      for (j=0; j<numRCEntries; j++) curRC[j].index = indicesRC[j];
      offset += numRCEntries;
    }

    costOfMatVec_ = 2.0 * ((double) numEntries_);
    if (hasUnitDiagonal_) costOfMatVec_ += 2.0 * ((double) numRC_);
    haveStructure_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int OskiSparseMultiply<OrdinalType, ScalarType>::initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, 
							       bool willKeepValues, bool checkStructure) {
// **Again, we may just have to store the OskiMatrix.  If this is true, we may
// have to extract some info from the mat vec objects anyway

    if (!haveStructure_) return(-1); // Must have structure first!

    matrixForValues_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;

    ScalarType * valuesRC;

    OrdinalType offset = 0;
    for (i=0; i<numRC_; i++) {
      int ierr = A.getValues(i, valuesRC);
      if (ierr<0) return(ierr);
      Entry * curRC = allEntries_+offset;
      OrdinalType numRCEntries = profile_[i];
      for (j=0; j<numRCEntries; j++) curRC[j].value = valuesRC[j];
      offset += numRCEntries;
    }
    haveValues_ = true;
    return(0);
  }

  //==============================================================================

// Dummy method

  template<typename OrdinalType, typename ScalarType>
  int OskiSparseMultiply<OrdinalType, ScalarType>::apply(const Vector<OrdinalType, ScalarType>& x,
                                                    Vector<OrdinalType, ScalarType> & y,
                                                    bool transA, bool conjA) const {return(-1);}

//=============================================================================
  template<typename OrdinalType, typename ScalarType>
  int OskiSparseMultiply<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
						    MultiVector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {
// ** Just call the wrapped Oski function here

    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA && !transA) return(-6); // Unsupported at this time
    //OSKI verifies dimensions
    //if (x.getNumRows()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    //if (y.getNumRows()!=numRows_) return(-4); // Number of rows in A not same as number of rows in y
    OrdinalType numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-5); // Not the same number of vectors in x and y

    OrdinalType * profile = profile_;

    int oskiReturn = 0;

    OskiMatrix<OrdinalType,ScalarType> *A_tunable = dynamic_cast<Kokkos::OskiMatrix<OrdinalType,ScalarType> *>(matrixForValues_);

    const OskiMultiVector<OrdinalType,ScalarType> * x_view = dynamic_cast<Kokkos::OskiMultiVector<OrdinalType,ScalarType> const *>(&x);

    OskiMultiVector <OrdinalType,ScalarType> *y_view = dynamic_cast<Kokkos::OskiMultiVector<OrdinalType,ScalarType> *>(&y);


    if (!conjA && !transA) {

	oskiReturn = oski_MatMult(A_tunable->getA_tunable(), OP_NORMAL, 1.0, x_view->getX_view(), 0.0, y_view->getX_view());
    }
    else if (!conjA && transA) {

	oskiReturn = oski_MatMult(A_tunable->getA_tunable(), OP_TRANS, 1.0, x_view->getX_view(), 0.0, y_view->getX_view());
    }
    else {//conjA && transA

	oskiReturn = oski_MatMult(A_tunable->getA_tunable(), OP_CONJ_TRANS, 1.0, x_view->getX_view(), 0.0, y_view->getX_view());
    }
    updateFlops(this->costOfMatVec_ * ((double) numVectors));
    return(oskiReturn);
  }

} // namespace Kokkos
#endif /* KOKKOS_OSKISPARSEMULTIPLY_H */
