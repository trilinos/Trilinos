
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

#ifndef KOKKOS_OPTSPARSEMATRIX_H
#define KOKKOS_OPTSPARSEMATRIX_H

#include "Kokkos_CompObject.hpp" 

//! Kokkos::OptSparseMatrix: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Kokkos::OptSparseMatrix enables the piecewise construction and use of real-valued double-precision sparse matrices
  where matrix entries are intended for row access.

  At this time, the primary function provided by Kokkos::OptSparseMatrix is matrix times vector and matrix 
  times multi-vector multiplication.  It is also possible to extract matrix rows from a constructed matrix.

  <b>Constructing Kokkos::OptSparseMatrix objects</b>

  Constructing Kokkos::OptSparseMatrix objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::OptSparseMatrix instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
  </ol>

  Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
  create new entries.

  <b> Counting Floating Point Operations </b>

  Each Kokkos::OptSparseMatrix object keeps track of the number
  of \e serial floating point operations performed using the specified object as the \e this argument
  to the function.  The Flops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate parallel performance
  numbers.  The ResetFlops() function resets the floating point counter.


*/    

namespace Kokkos {

  template<typename OrdinalType, typename ScalarType>
  class OptSparseMatrix: public CompObject {
  public:

    //@{ \name Constructors/Destructor.
    //! OptSparseMatrix constuctor with variable number of indices per row.
    OptSparseMatrix(void);
  
    //! Copy constructor.
    OptSparseMatrix(const OptSparseMatrix& Matrix);
	
    //! OptSparseMatrix Destructor
    virtual ~OptSparseMatrix();
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

    //@{ \name Computational methods.
	
    //! Returns the result of a Kokkos_OptSparseMatrix multiplied by a vector x in y.
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

    //! Returns the result of a Kokkos_OptSparseMatrix multiplied by multiple vectors in x, results in y.
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

} // namespace Kokkos
#endif /* KOKKOS_OPTSPARSEMATRIX_H */
