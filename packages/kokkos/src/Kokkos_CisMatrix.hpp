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

#ifndef KOKKOS_CISMATRIX_H
#define KOKKOS_CISMATRIX_H


namespace Kokkos {

//! Kokkos::CisMatrix: Kokkos compressed index sparse matrix base class.

/*! The Kokkos::CisMatrix specifies the interface that any sparse matrix interfacing to the Kokkos 
  Operators classes must implement.

  At this time, the primary function provided by Kokkos::CisMatrix is matrix times vector and matrix 
  times multi-vector multiplication and solves with single or multiple right-hand-sides.

*/    

  template<typename OrdinalType, typename ScalarType>
  class CisMatrix {
  public:

    //@{ \name Constructors/Destructor.

    //! CisMatrix Destructor
    virtual ~CisMatrix();
    //@}

    //@{ \name Matrix entry access methods.

    //! Returns number of entries in ith row/column, and pointer to an array of these indices.
    /*! Extract the number of entries and a pointer to the indices in the ith row/column of the matrix.  Note that
        the indices are not copied by this method.  Memory allocation is handled by the matrix object itself.

	\param i (In) The row (if isRowOriented() is true) or column that should be returned.
	\param numEntries (Out) The number of entries in the ith row/column.
	\param indices (Out) A pointer to the list of indices in the ith row/column.

	\return Integer error code, set to 0 if successful.
    */
    virtual int getIndices(OrdinalType i, OrdinalType & numEntries, OrdinalType *& indices) const = 0;

    //! Returns a pointer to an array of values for the ith row/column.
    /*! Extract the values in the ith row/column of the matrix.  Note that
        the values are not copied by this method.  Memory allocation is handled by the matrix object itself.

	\param i (In) The row (if isRowOriented() is true) or column that should be returned.
	\param numEntries (Out) The number of entries in the ith row/column.
	\param indices (Out) A pointer to the list of indices in the ith row/column.

	\return Integer error code, set to 0 if successful.
    */
    virtual int getValues(OrdinalType i, ScalarType *& values) const = 0;
	
	
    //@}

    //@{ \name Matrix Attribute access methods.
	
    //! Returns true if the compressed index matrix should be interpreted as a row matrix.
    virtual bool isRowOriented() const = 0;
	
    //! Number of rows
    virtual OrdinalType numRows() const = 0;
	
    //! Number of columns
    virtual OrdinalType numCols() const = 0;
	
    //! Number of matrix entries
    virtual OrdinalType numEntries() const = 0;
	
    //@}
  };

} // namespace Kokkos
#endif /* KOKKOS_CISMATRIX_H */
