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

#ifndef KOKKOS_VECTOR_H
#define KOKKOS_VECTOR_H


namespace Kokkos {

//! Kokkos::Vector: Kokkos compressed index sparse matrix base class.

/*! The Kokkos::Vector specifies the interface that any vector class interfacing to the Kokkos 
  Operators classes must implement.

  At this time, the primary function provided by Kokkos::Vector is matrix times vector and matrix 
  times multi-vector multiplication and solves with single or multiple right-hand-sides.

*/    

  template<typename OrdinalType, typename ScalarType>
  class Vector {
  public:

    //@{ \name Constructors/Destructor.

    //! Vector Destructor
    virtual ~Vector();
    //@}

    //@{ \name Vector access methods.

    //! Returns a pointer to an array of values in the vector.
    /*! Extract the values in the ith row/column of the matrix.  Note that
        the values are not copied by this method.  Memory allocation is handled by the matrix object itself.

	\param i (In) The row (if isRowOriented() is true) or column that should be returned.
	\param numEntries (Out) The number of entries in the ith row/column.
	\param indices (Out) A pointer to the list of indices in the ith row/column.

	\return Integer error code, set to 0 if successful, set to -1 if not implemented.
    */
    virtual int getValues(ScalarType *& values) const = 0;
	
	
    //@}

    //@{ \name Vector Attribute access methods.
	
    //! Length of vector
    virtual OrdinalType length() const = 0;
	
    //@}
  };

} // namespace Kokkos
#endif /* KOKKOS_VECTOR_H */
