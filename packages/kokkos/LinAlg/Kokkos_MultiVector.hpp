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

#ifndef KOKKOS_MULTIVECTOR_H
#define KOKKOS_MULTIVECTOR_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>

namespace Kokkos {

//! Kokkos::MultiVector: Kokkos default implementation of the abstract Kokkos::MultiVector.

/*! The Kokkos::MultiVector provides a default implementation of Kokkos::MultiVector interface.

  At this time, the primary functions provided by Kokkos::MultiVector is wrapping of a
  set of dense vector and providing access to multivector entries.  Two basic
  categories of data structures are supported:
  <ol>
  <li> MultiVector is described by an array of pointers:  In this situation, the 
       ith entry of this array of pointers
       is the starting address of a contiguous array of values for the ith vector 
       in the multivector.  The storage mode
       will be assumed if the getIsStrided() method returns false.
  <li> MultiVector is a regular strided two-dimensional array of values:  In this situation, 
       the increment between 
       elements in the row and column dimensions is specified by the getRowInc() 
       and getColInc() methods. This is
       a very general mechanism for describing strided access.  Typical situations include:
       <ul>
       <li> getRowInc() = getNumCols(), getColInc() = 1 - column entries are contiguous.  
       <li> getRowInc() = 1, getColInc() = getNumRows() - row entries are contiguous.  
       </ul>
       However, this mechanism also allows extraction of array subsections, real or
       imaginary parts from 
       complex-valued arrays.

       This storage mode will be assumed if getIsStrided() returns true.  The base 
       address for the 2D array
       will be obtain by call getValues() with the argument equal to 0.

  </ol>

*/    

  template<class Scalar, class Node = DefaultNode::DefaultNodeType>
  class MultiVector {
    public:
      typedef Scalar  ScalarType;
      typedef Node    NodeType;

      //! @name Constructors/Destructor

      //@{

      //! Default constructor
      MultiVector(Node &node = DefaultNode::getDefaultNode())
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0) {}

      //! Copy constructor.
      MultiVector(const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_) {}

      //! MultiVector Destructor
      ~MultiVector() {}

      //@}

      //! @name Initialization methods

      //@{

      //! Initialize using a two-dimensional array
      /*!
        This interface supports multivectors that are stored as 2D arrays, or subsections of one.
        \param numRows (In)  Number of rows in multivector (length of each vector).
        \param numCols (In)  Number of columns in multivector (number of vectors).
        \param values (In)  Pointer to the first entry in the multivector.  Subsequent column 
        entries are spaced a distance of getColInc().  Subsequent row entries
        are spaced by getRowInc() increments.
        \param rowInc (In) The increment between two elements in a row of the multivector.  
        Typically this value should be set to numRows.
        \param colInc (In) The increment between two elements in a column of the multivector.  
        Typically this value should be set to 1, which is the default value.

        \return Integer error code, set to 0 if successful.
        */
      int initializeValues(size_type numRows, size_type numCols, 
                           Teuchos::ArrayRCP<Scalar> values,
                           size_type stride)
      {
        numRows_ = numRows;
        numCols_ = numCols;
        stride_ = stride;
        contigValues_ = values;
        return(0);
      };

      //@}

      //! @name Multivector entry access methods

      //@{

      //! Returns a copy of the ArrayRCP passed to initializeValues().
      Teuchos::ArrayRCP<Scalar>
      getValuesNonConst() {
        return contigValues_;
      }

      //! Returns a copy of the ArrayRCP passed to initializeValues().
      Teuchos::ArrayRCP<const Scalar>
      getValues() const {
        return contigValues_;
      }

      //! Returns a pointer to an array of values for the ith column of the multivector.
      /*! Extract a pointer to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
        handled by the multivector object itself.  Also, if the getIsStrided() method returns
        true, then the getColInc() should be used to access values in the ith column
        of the multivector, especially if getColInc() != 1.

        \param i (In) The column that should be returned.
        */
      Teuchos::ArrayRCP<Scalar>
      getValuesNonConst(size_type i) {
        TEST_FOR_EXCEPTION((contigValues_ == Teuchos::null) || // No data to return
                           i < 0 || i >= numRows_, // Out of range
                           std::runtime_error, 
                           Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
        return contigValues_.persistingView(stride_*i,numRows_);
      };

      //! Returns a pointer to an array of values for the ith column of the multivector.
      /*! Extract a pointer to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
        handled by the multivector object itself.  Also, if the getIsStrided() method returns
        true, then the getColInc() should be used to access values in the ith column
        of the multivector, especially if getColInc() != 1.

        \param i (In) The column that should be returned.
        */
      Teuchos::ArrayRCP<const Scalar>
      getValues(size_type i) const {
        TEST_FOR_EXCEPTION((contigValues_ == Teuchos::null) || // No data to return
                           i < 0 || i >= numRows_, // Out of range
                           std::runtime_error, 
                           Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
        return contigValues_.persistingView(stride_*i,numRows_);
      };

      //@}

      //! @name MultiVector Attribute access methods

      //@{

      Node & getNode() const {return node_;}

      //! Number of rows
      size_type getNumRows() const {return(numRows_);};

      //! Number of columns
      size_type getNumCols() const{return(numCols_);};

      //! Increment between entries in a row of the multivector, normally = numRows().
      size_type getStride() const {return(stride_);};

      //@}

    protected:
      Node &node_;

      Teuchos::ArrayRCP<Scalar> contigValues_;

      bool dataInitialized_;
      size_type numRows_, numCols_;
      size_type stride_;
  };

} // namespace Kokkos

#endif /* KOKKOS_MULTIVECTOR_H */
