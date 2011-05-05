//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>

namespace Kokkos {

  /** \brief Data structure for vector and multivector data.
  
    The primary functions provided by MultiVector is wrapping of a
    set of dense vectors and providing access to multivector entries.
  
    MultiVector is a regular strided two-dimensional array of values.  The
    increment between elements in the row and column dimensions is specified by
    getStride().
  */    
  template<class Scalar, class Node = DefaultNode::DefaultNodeType>
  class MultiVector {
    public:
      typedef Scalar  ScalarType;
      typedef Node    NodeType;

      //! @name Constructors/Destructor

      //@{

      //! Default constructor
      MultiVector(RCP<Node> node)
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0) {
      }

      //! Copy constructor.
      MultiVector(const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_) {
      }

      //! MultiVector Destructor
      ~MultiVector() {
      }

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
      void initializeValues(size_t numRows, size_t numCols, 
                            const ArrayRCP<Scalar> &values,
                            size_t stride) {
        numRows_ = numRows;
        numCols_ = numCols;
        stride_ = stride;
        contigValues_ = values;
      };

      //@}

      //! @name Multivector entry access methods

      //@{

      //! Returns a copy of the ArrayRCP passed to initializeValues().
      ArrayRCP<Scalar>
      getValuesNonConst() {
        return contigValues_;
      }

      //! Returns a copy of the ArrayRCP passed to initializeValues().
      ArrayRCP<const Scalar>
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
      ArrayRCP<Scalar>
      getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOS_DEBUG
        TEST_FOR_EXCEPTION( !( (contigValues_ != null) &&  // Data to return
                               ( (i > 0 || i == 0) && i < numCols_)    // In range
                             ), std::runtime_error, 
                             Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
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
      ArrayRCP<const Scalar>
      getValues(size_t i) const {
#ifdef HAVE_KOKKOS_DEBUG
        TEST_FOR_EXCEPTION( !( (contigValues_ != null) &&  // Data to return
                               ( (i > 0 || i == 0) && i < numCols_)    // In range
                             ), std::runtime_error, 
                             Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
        return contigValues_.persistingView(stride_*i,numRows_);
      };

      //@}

      //! @name MultiVector Attribute access methods

      //@{

      RCP<Node> getNode() const {return node_;}

      //! Number of rows
      size_t getNumRows() const {return(numRows_);};

      //! Number of columns
      size_t getNumCols() const{return(numCols_);};

      //! Increment between entries in a row of the multivector, normally = numRows().
      size_t getStride() const {return(stride_);};

      //@}

    protected:
      RCP<Node> node_;

      ArrayRCP<Scalar> contigValues_;

      bool dataInitialized_;
      size_t numRows_, numCols_, stride_;
  };

} // namespace Kokkos

#endif /* KOKKOS_MULTIVECTOR_H */
