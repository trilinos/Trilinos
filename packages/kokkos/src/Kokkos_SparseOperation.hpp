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

#ifndef KOKKOS_SPARSEOPERATION_H
#define KOKKOS_SPARSEOPERATION_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CompObject.hpp"

namespace Kokkos {

//! Kokkos::SparseOperation: A virtual class that defines the interface for sparse matrix operations.

/*! The Kokkos::SparseOperation specifies the basic functionality that any Kokkos sparse operation must satisfy. 
    This class is templated on the ordinal (integer) and 
    scalar (floating point) types, so it can interface with any reasonable data type.

  <b>Constructing Kokkos::SparseOperation objects</b>

  Constructing Kokkos::SparseOperation objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create Kokkos::SparseOperation instance:  There are several classes that implement this interface.
       The constructor from one of these classes must be called.  However, once constructed, any of these
       objected can be treated from thence on as a member of the Kokkos::SparseOperation base class.
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

  Each Kokkos::SparseOperation object keeps track of the number
  of floating point operations performed using the specified object as the \e this argument
  to the function.  The getFlops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Kokkos::Time class, one can get accurate  performance
  numbers.  The resetFlops() function resets the floating point counter.


*/    

  template<typename OrdinalType, typename ScalarType>
  class SparseOperation: public CompObject {
  public:

    //@{ \name Constructors/Destructor.

    //! SparseOperation Constructor
    SparseOperation(): CompObject(){};

    //! SparseOperation Destructor
    virtual ~SparseOperation(){};
    //@}
    //@{ \name Abstract Kokkos::CisMatrix Interface Initialization Methods
 
    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
             about the matrix can be obtained via this interface.
      \param willKeepStructure (In) If set to true, the user is asserting that the strucuture of the matrix, as
             defined in the getIndices() method of the CisMatrix object A will be kept.  Specifically, the pointer to an 
	     array of indices returned for each i in that method will continue to point to valid index data.

      \return Integer error code, set to 0 if successful.
    */
    virtual int initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepStructure) = 0;
 
    //! Initialize values of matrix
    /*!
      This interface supports matrices that implement the Kokkos::CisMatrix matrix interface
.
      \param A (In)  An instance of a class that implements the Kokkos::CisMatrix.  All necessary information
             about the matrix can be obtained via this interface.
      \param willKeepValues (In) If set to true, the user is asserting that the strucuture of the matrix, as
             defined in the getIndices() method of the CisMatrix object A will be kept.  Specifically, the pointer to an 
	     array of indices returned for each i in that method will continue to point to valid index data
.
      \param checkStructure (In) If set to true, the structure of A will be checked against the structure of
             the matrix passed in to the initializeStructure() methods.

      \return Integer error code, set to 0 if successful, returns - 1 if checkStructure is true and structure is changed.
    */
    virtual int initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepValues,
				 bool checkStructure) = 0;
 
    //@}

    //@{ \name Computational methods.
	
    //! Returns the result of a Kokkos_SparseOperation multiplied by a vector x in y.
    /*! 
      \param x (In) A Kokkos::Vector to operate on.
      \param y (Out) A Kokkos::Vector containing results.
      \param transA (In) If true, apply the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, applythe conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType>& y, 
		      bool transA, bool conjA) const = 0;

    //! Returns the result of a Kokkos_SparseOperation multiplied by multiple vectors in x, results in y.
    /*! 
      \param x (In) A Kokkos::MultiVector to operate on.
      \param y (Out) A Kokkos::MultiVector containing results.
      \param transA (In) If true, apply the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, apply the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y, 
		      bool transA, bool conjA) const = 0;
    //@}
	
    //@{ \name Operator attribute access methods.

    //! Returns true if this implementation of Kokkos::SparseOperation can benefit from the user keeping the passed in structure.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        analysis of structure.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepStructure parameter 
	passed in to the initializeStructure() method.
    */
    virtual bool getCanUseStructure() const = 0;

    //! Returns true if this implementation of Kokkos::SparseOperation can benefit from the user keeping the passed in values.
    /*! Some implementations of optimized kernels do not rely on the user's data except for the initial 
        copying of values.  Other implementations, in order to reduce memory requirements, may find it
	beneficial to rely on the user's data.  Since it is very possible that the user would retain this
	data anyway, we want to allow for this possibility.  This method is related to the willKeepValues parameter 
	passed in to the initializeValues() method.
    */
    virtual bool getCanUseValues() const = 0;

    //! Returns a reference to the most recent Kokkos::CisMatrix that was passed into the \e this object.
    virtual const CisMatrix<OrdinalType, ScalarType> & getMatrix() const = 0;
	

    //! Returns a reference to the left Kokkos::Permutation object.
    virtual const Permutation<OrdinalType, ScalarType> & getLeftPermutation() const = 0;

    //! Returns a reference to the right Kokkos::Permutation object, which is the identity for this implementation.
    virtual const Permutation<OrdinalType, ScalarType> & getRightPermutation() const = 0;

    //@}
  
  };


} // namespace Kokkos
#endif /* KOKKOS_SPARSEOPERATION_H */
