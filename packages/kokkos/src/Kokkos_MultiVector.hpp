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

#ifndef KOKKOS_MULTIVECTOR_H
#define KOKKOS_MULTIVECTOR_H

#include "Kokkos_ConfigDefs.hpp"

namespace Kokkos {

//! Kokkos::MultiVector: Kokkos multivector base class, used by Kokkos operator classes.

/*! The Kokkos::MultiVector specifies the interface that any multivector interfacing to the Kokkos 
  Operators classes must implement.

  At this time, the primary function provided by Kokkos::MultiVector is access to 
  multivector entries.  Two basic categories of data structures are supported:
  <ol>
  <li> MultiVector is described by an array of pointers:  
       In this situation, the ith entry of this array of pointers
       is the starting address of a contiguous array of values for 
       the ith vector in the multivector.  The storage mode
       will be assumed if the getIsStrided() method returns false.
  <li> MultiVector is a regular strided two-dimensional array of values:  
       In this situation, the increment between 
       elements in the row and column dimensions is specified by the 
       getRowInc() and getColInc() methods. This is
       a very general mechanism for describing strided access.  Typical situations include:
       <ul>
       <li> getRowInc() = getNumCols(), getColInc() = 1 - column entries are contiguous.  
       <li> getRowInc() = 1, getColInc() = getNumRows() - row entries are contiguous.  
       </ul>
       However, this mechanism also allows extraction of array subsections, 
       real or imaginary parts from complex-valued arrays.

       This storage mode will be assumed if getIsStrided() returns true.  
       The base address for the 2D array
       will be obtain by call getValues() with the argument equal to 0.

  </ol>

*/    

  template<typename OrdinalType, typename ScalarType>
  class MultiVector {
  public:

    //@{ \name Constructors/Destructor.

    //! MultiVector Destructor
    virtual ~MultiVector(){};
    //@}

    //@{ \name Multivector entry access methods.

    //! Returns an array of pointers such that the ith pointer points to an array of values in the ith column of the multivector.
    /*! Extract an array of pointers such that the ith pointer points to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the multivector object itself.  
    */
    virtual ScalarType ** getValues() const = 0;

    //! Returns a pointer to an array of values in the ith column of the multivector.
    /*! Extract a pointer to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the multivector object itself.  Also, if the getIsStrided() method returns
	true, then the getColInc() should be used to access values in the ith column
	of the multivector, especially if getColInc() != 1.

	\param i (In) The column that should be returned.

    */
    virtual ScalarType * getValues(OrdinalType i) const = 0;
	
	
    //@}

    //@{ \name MultiVector Attribute access methods.
	
    //! Number of rows
    virtual OrdinalType getNumRows() const = 0;
	
    //! Number of columns
    virtual OrdinalType getNumCols() const = 0;
	
    //! Indicates whether or not array is strided
    virtual bool getIsStrided() const = 0;
	
    //! Increment between entries in a row of the multivector, normally = numRows().
    virtual OrdinalType getRowInc() const = 0;
	
    //! Increment between entries in a column of the multivector, normally = 1.
    virtual OrdinalType getColInc() const = 0;
	
    //@}
  };

} // namespace Kokkos
#endif /* KOKKOS_MULTIVECTOR_H */
