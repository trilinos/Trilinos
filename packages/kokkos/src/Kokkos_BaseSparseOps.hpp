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

namespace Kokkos {

  template<typename OrdinalType, typename ScalarType>
  class BaseSparseOps: public CompObject {
  public:

    //@{ \name Constructors/Destructor.
    //! BaseSparseOps constuctor with variable number of indices per row.
    BaseSparseOps(void);
  
    //! Copy constructor.
    BaseSparseOps(const BaseSparseOps& Matrix);
	
    //! BaseSparseOps Destructor
    virtual ~BaseSparseOps();
    //@}
    //@{ \name Abstract Kokkos::CisMatrix Interface Initialization Methods
 
    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
                     about the matrix can be obtained via this interface.
      \return Integer error code, set to 0 if successful.
    */
    int initializeStructure(CisMatrix<OrdinalType, ScalarType> const& A);
 
    //! Initialize values of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
                     about the matrix can be obtained via this interface.
      \param checkStructure (In) If set to true, the structure of A will be checked against the structure of
                                 the matrix passed in to the initializeStructure() methods.  This parameter is false
				 by default.
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(CisMatrix<OrdinalType, ScalarType> const& A, bool checkStructure = false);
 
    //@}

    //@{ \name Computational methods.
	
    //! Returns the result of a Kokkos_BaseSparseOps multiplied by a vector x in y.
    /*! 
      \param xLength (In) Length of vector x.
      \param x (In) A Vector to multiply by.
      \param yLength (In) Length of vector y.
      \param y (Out) A Vector containing result.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    int apply(OrdinalType xLength, ScalarType * x, OrdinalType yLength, ScalarType * y, 
	      bool transA = false, bool conjA = false) const;

    //! Returns the result of a Kokkos_BaseSparseOps multiplied by multiple vectors in x, results in y.
    /*! 
      \param numVectors (In) Number of vectors.
      \param xLength (In) Length of vectors in x.
      \param x (In) An array of pointers to vectors to multiply by.
      \param yLength (In) Length of vector y.
      \param y (Out) A array of pointers to vectors containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    int apply(OrdinalType numVectors, OrdinalType xLength, ScalarType ** x, OrdinalType yLength, ScalarType ** y, 
	      bool transA = false, bool conjA = false) const;


 
    //! Returns the result of a triangular solve of vector x in y.
    /*!
      \param xLength (In) Length of vector x.
      \param x (In) A Vector to multiply by.
      \param yLength (In) Length of vector y.
      \param y (Out) A Vector containing result.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
      \param upper (In) If true, solve Uy = x, otherwise solve Ly = x.
      \param unitDiagonal (In) If true, assume diagonal is unit (whether it's stored or not).
 
      \return Integer error code, set to 0 if successful.

      \warning Any implementation of this method must work for the case when x and y are the same vector!
    */
    int applyInverse(OrdinalType xLength, ScalarType * x, OrdinalType yLength, ScalarType * y,
		     bool transA = false, bool conjA = false, bool upper = true, bool unitDiagonal = true) const;
 
    //! Returns the result of a triangular solve of multiple vectors in x, results in y.
    /*!
      \param numVectors (In) Number of vectors.
      \param xLength (In) Length of vectors in x.
      \param x (In) An array of pointers to vectors to multiply by.
      \param yLength (In) Length of vector y.
      \param y (Out) A array of pointers to vectors containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
      \param upper (In) If true, solve Uy = x, otherwise solve Ly = x.
      \param unitDiagonal (In) If true, assume diagonal is unit (whether it's stored or not).
 
      \return Integer error code, set to 0 if successful.

      \warning Any implementation of this method must work for the case when x and y are the same vector!
    */
    int applyInverse(OrdinalType numVectors, OrdinalType xLength, ScalarType ** x, OrdinalType yLength, ScalarType ** y,
		     bool transA = false, bool conjA = false, bool upper = true, bool unitDiagonal = true) const;
 
    //@}
	
    //@{ \name Matrix Attribute access methods.
	
    //! Number of rows
    OrdinalType numRows() const {return(numRows_);};
	
    //! Number of columns
    OrdinalType numCols() const {return(numCols_);};
	
    //! Number of matrix entries
    OrdinalType numEntries() const {return(numEntries_);};
	
    //@}
  
  protected:
    bool allocated() const {return(allocated_);};
    int setAllocated(bool flag) {allocated_ = flag; return(0);};
	
    void initializeDefaults();
    int allocate();

    bool allocated_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numEntries_;

    ScalarType ** values_;
    ScalarType * allValues_;

    OrdinalType ** indices_;
    OrdinalType * allIndices_;

    OrdinalType * pntr_;
    OrdinalType ** profile_;

  };

//==============================================================================
template<typename OrdinalType, typename ScalarType>
 BaseSparseOps<OrdinalType, ScalarType>::BaseSparseOps() 
  : CompObject(),
    allocated_(false),
    numRows_(0),
    numCols_(0),
    numEntries_(0),
    values_(0),
    allValues_(0),
    indices_(0),
    allIndices_(0),
    pntr_(0),
    profile_(0)
{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
BaseSparseOps<OrdinalType, ScalarType>::BaseSparseOps(const BaseSparseOps<OrdinalType, ScalarType> &matrix) 
  : CompObject(matrix),
    allocated_(matrix.allocated_),
    numRows_(matrix.numRows_),
    numCols_(matrix.numCols_),
    numEntries_(matrix.numEntries_),
    values_(matrix.values_),
    allValues_(matrix.allValues_),
    indices_(matrix.indices_),
    allIndices_(matrix.allIndices_),
    pntr_(matrix.pntr_),
    profile_(matrix.profile_)

{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
BaseSparseOps<OrdinalType, ScalarType>::~BaseSparseOps(){}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
void BaseSparseOps<OrdinalType, ScalarType>::initializeDefaults() { // Initialize all attributes that have trivial default values
  return;
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int BaseSparseOps<OrdinalType, ScalarType>::allocate() { // Initialize all attributes that have trivial default values
  return;
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int BaseSparseOps<OrdinalType, ScalarType>::apply(OrdinalType xLength, ScalarType * x, 
						    OrdinalType yLength, ScalarType * y,
						    bool transA, bool conjA) const {
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int BaseSparseOps<OrdinalType, ScalarType>::apply(OrdinalType numVectors, 
						    OrdinalType xLength, ScalarType ** x, 
						    OrdinalType yLength, ScalarType ** y,
						    bool transA, bool conjA) const {
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int BaseSparseOps<OrdinalType, ScalarType>::applyInverse(OrdinalType xLength, ScalarType * x, 
							   OrdinalType yLength, ScalarType * y,
							   bool transA, bool conjA, bool upper, 
							   bool unitDiagonal) const {
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int BaseSparseOps<OrdinalType, ScalarType>::applyInverse(OrdinalType numVectors, 
							   OrdinalType xLength, ScalarType ** x, 
							   OrdinalType yLength, ScalarType ** y,
							   bool transA, bool conjA, bool upper, 
							   bool unitDiagonal) const {
}

} // namespace Kokkos
#endif /* KOKKOS_BASESPARSEOPS_H */
