//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
#include <Teuchos_Assert.hpp>
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
      , stride_(0)
      , origNumRows_(0)
      , origNumCols_(0) {
      }

      //! Copy constructor.
      MultiVector(const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_)
      , origNumRows_(source.origNumRows_)
      , origNumCols_(source.origNumCols_) {
      }

      //! MultiVector Destructor
      ~MultiVector() {
      }

      //@}

      //! @name Initialization methods

      //@{

      /// \brief Initialize a multivector that is not a view of another multivector.
      ///
      /// \param numRows [in] Number of rows in the multivector.
      /// \param numCols [in] Number of columns in the multivector.
      /// \param values [in] Array of the multivector's entries,
      ///   stored in column-major order with stride \c stride between
      ///   columns.  If you are familiar with the BLAS or LAPACK,
      ///   <tt>stride</tt> here corresponds to "LDA" (<i>l</i>eading
      ///   <i>d</i>imension of the matrix A).
      /// \param stride [in] The stride (number of entries between)
      ///   adjacent columns of the multivector.
      void initializeValues(size_t numRows, size_t numCols, 
                            const ArrayRCP<Scalar> &values,
                            size_t stride) {
        numRows_ = numRows;
        numCols_ = numCols;
        stride_ = stride;
        contigValues_ = values;
	origNumRows_ = numRows;
	origNumCols_ = numCols;
      }

      /// \brief Initialize a multivector that <i>is</i> a view of another multivector.
      ///
      /// \param numRows [in] Number of rows in the multivector.
      /// \param numCols [in] Number of columns in the multivector.
      /// \param values [in] Array of the multivector's entries,
      ///   stored in column-major order with stride \c stride between
      ///   columns.  If you are familiar with the BLAS or LAPACK,
      ///   <tt>stride</tt> here corresponds to "LDA" (<i>l</i>eading
      ///   <i>d</i>imension of the matrix A).
      /// \param stride [in] The stride (number of entries between)
      ///   adjacent columns of the multivector.
      /// \param origNumRows [in] Number of rows in the "original"
      ///   multivector (of which this multivector will be a view).
      /// \param origNumCols [in] Number of columns in the "original"
      ///   multivector (of which this multivector will be a view).
      void 
      initializeValues (size_t numRows, 
			size_t numCols,
			const ArrayRCP<Scalar> &values,
			size_t stride,
			size_t origNumRows,
			size_t origNumCols)
      {
        numRows_ = numRows;
        numCols_ = numCols;
        stride_ = stride;
        contigValues_ = values;
	origNumRows_ = origNumRows;
	origNumCols_ = origNumCols;
      }
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

      //! Return a nonconst view of the data in the i-th column of the multivector.
      ArrayRCP<Scalar>
      getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION( !( (contigValues_ != null) &&  // Data to return
                               ( (i > 0 || i == 0) && i < numCols_)    // In range
                             ), std::runtime_error, 
                             Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
        return contigValues_.persistingView(stride_*i,numRows_);
      };

      //! Return a const view of the data in the i-th column of the multivector.
      ArrayRCP<const Scalar>
      getValues(size_t i) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION( !( (contigValues_ != null) &&  // Data to return
                               ( (i > 0 || i == 0) && i < numCols_)    // In range
                             ), std::runtime_error, 
                             Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
        return contigValues_.persistingView(stride_*i,numRows_);
      };

      //@}
      //! @name View "constructors"
      //@{

      /// \brief A const offset view of the multivector.
      ///
      /// \param newNumRows [in] Number of rows in the view.
      /// \param newNumCols [in] Number of columns in the view.
      /// \param offsetRow [in] Zero-based index of the starting row of the view.
      /// \param offsetCol [in] Zero-based index of the starting column of the view.
      const MultiVector<Scalar,Node>
      offsetView (size_t newNumRows,
		  size_t newNumCols,
		  size_t offsetRow,
		  size_t offsetCol) const
      {
	MultiVector<Scalar,Node> B (this->getNode ());

	TEUCHOS_TEST_FOR_EXCEPTION(
          offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
	  std::invalid_argument,
	  Teuchos::typeName (*this) << "::offsetView: offset row or column are out of bounds.  "
	  "The original multivector has dimensions " << this->getOrigNumRows () << " x " << this->getOrigNumCols () 
	  << ", but your requested offset row and column are " << offsetRow << ", " << offsetCol << ".");
	TEUCHOS_TEST_FOR_EXCEPTION(
          newNumRows >= this->getOrigNumRows () || newNumCols >= this->getOrigNumCols (),
	  std::invalid_argument,
	  Teuchos::typeName (*this) << "::offsetView: new dimensions are out of bounds.  "
	  "The original multivector has dimensions " << this->getOrigNumRows () << " x " << this->getOrigNumCols () 
	  << ", but your requested new dimensions are " << newNumRows << " x " << newNumCols << ".");

	// Starting position of the view of the data.
	const size_t startPos = offsetRow + this->getStride () * offsetCol;
	// Length of the view of the data.
	const size_t len = (newNumCols > 0) ? (this->getStride () * newNumCols - offsetRow) : 0;

	B.initializeValues (newNumRows,
			    newNumCols,
			    contigValues_.persistingView (startPos, len),
			    this->getStride (),
			    this->getOrigNumRows (),
			    this->getOrigNumCols ());
	return B;
      }

      /// \brief A nonconst offset view of the multivector.
      ///
      /// \param newNumRows [in] Number of rows in the view.
      /// \param newNumCols [in] Number of columns in the view.
      /// \param offsetRow [in] Zero-based index of the starting row of the view.
      /// \param offsetCol [in] Zero-based index of the starting column of the view.
      MultiVector<Scalar,Node>
      offsetViewNonConst (size_t newNumRows,
			  size_t newNumCols,
			  size_t offsetRow,
			  size_t offsetCol)
      {
	MultiVector<Scalar,Node> B (this->getNode ());

	TEUCHOS_TEST_FOR_EXCEPTION(
          offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
	  std::invalid_argument,
	  Teuchos::typeName (*this) << "::offsetViewNonConst: offset row or column are out of bounds.  "
	  "The original multivector has dimensions " << this->getOrigNumRows () << " x " << this->getOrigNumCols () 
	  << ", but your requested offset row and column are " << offsetRow << ", " << offsetCol << ".");
	TEUCHOS_TEST_FOR_EXCEPTION(
          newNumRows >= this->getOrigNumRows () || newNumCols >= this->getOrigNumCols (),
	  std::invalid_argument,
	  Teuchos::typeName (*this) << "::offsetViewNonConst: new dimensions are out of bounds.  "
	  "The original multivector has dimensions " << this->getOrigNumRows () << " x " << this->getOrigNumCols () 
	  << ", but your requested new dimensions are " << newNumRows << " x " << newNumCols << ".");

	// Starting position of the view of the data.
	const size_t startPos = offsetRow + this->getStride () * offsetCol;
	// Length of the view of the data.
	const size_t len = (newNumCols > 0) ? (this->getStride () * newNumCols - offsetRow) : 0;

	B.initializeValues (newNumRows,
			    newNumCols,
			    contigValues_.persistingView (startPos, len),
			    this->getStride (),
			    this->getOrigNumRows (),
			    this->getOrigNumCols ());
	return B;
      }

      //@}
      //! @name MultiVector Attribute access methods
      //@{

      //! Node accessor
      RCP<Node> getNode() const {return node_;}

      //! Number of rows
      size_t getNumRows() const {return(numRows_);};

      //! Number of columns
      size_t getNumCols() const{return(numCols_);};

      //! Increment between entries in a row of the multivector, normally = numRows().
      size_t getStride() const {return(stride_);};

      /// \brief "Original" number of rows (of the multivector of
      ///   which <tt>*this</tt> is a view).
      ///
      /// If this multivector is <i>not</i> a view of another
      /// multivector, then this method just returns the number of
      /// rows.
      size_t getOrigNumRows() const {return(origNumRows_);};

      /// \brief "Original" number of columns (of the multivector of
      ///   which <tt>*this</tt> is a view).
      ///
      /// If this multivector is <i>not</i> a view of another
      /// multivector, then this method just returns the number of
      /// columns.
      size_t getOrigNumCols() const{return(origNumCols_);};

      //@}

    protected:
      RCP<Node> node_;

      ArrayRCP<Scalar> contigValues_;

      bool dataInitialized_;
      size_t numRows_, numCols_, stride_, origNumRows_, origNumCols_;
  };

} // namespace Kokkos

#endif /* KOKKOS_MULTIVECTOR_H */
