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

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TypeNameTraits.hpp"

namespace KokkosClassic {

  /** \brief Data structure for local vector and multivector data.

    \tparam Scalar The type of entries in the (multi)vector.
    \tparam Node The Kokkos Node type.

    MultiVector contains the local data for a distributed dense vector
    (Tpetra::Vector) or multivector (Tpetra::MultiVector).  (A
    multivector is a set of one or more vectors, all of which have the
    same number of entries and parallel distribution.)

    This class stores its data in column-major order, just like the
    BLAS or LAPACK.  getStride() specifies the stride (number of
    entries) between columns.  Consecutive elements in the same column
    are stored contiguously.

    Different MultiVector instances may view the same data.  ArrayRCP
    arbitrates ownership of that shared data, so that the array is not
    deallocated until all views have gone out of scope.
  */
  template<class Scalar, class Node = DefaultNode::DefaultNodeType>
  class MultiVector {
  public:
    //! @name Typedefs
    //@{

    //! The type of entries in the (multi)vector.
    typedef Scalar  ScalarType;
    //! The Kokkos Node type.
    typedef Node    NodeType;

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Default constructor
    MultiVector (const Teuchos::RCP<Node>& node)
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0)
      , origNumRows_(0)
      , origNumCols_(0) {
    }

    //! Copy constructor (shallow copy only).
    MultiVector (const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_)
      , origNumRows_(source.origNumRows_)
      , origNumCols_(source.origNumCols_) {
    }

    //! Destructor
    ~MultiVector() {
    }

    //@}
    //! @name Initialization methods
    //@{

    /// \brief Initialize a multivector that is not a view of another multivector.
    ///
    /// \param numRows [in] Number of rows in the multivector.
    /// \param numCols [in] Number of columns in the multivector.
    /// \param values [in] Array of the multivector's entries, stored
    ///   in column-major order with stride \c stride between columns.
    ///   The array is not copied; this method just keeps a reference.
    /// \param stride [in] The column stride, that is, the number of
    ///   entries between adjacent columns of the multivector.  If you
    ///   are familiar with the BLAS or LAPACK, <tt>stride</tt> here
    ///   corresponds to "LDA" (<i>l</i>eading <i>d</i>imension of the
    ///   matrix A).
    void
    initializeValues (size_t numRows,
                      size_t numCols,
                      const Teuchos::ArrayRCP<Scalar> &values,
                      size_t stride)
    {
      numRows_ = numRows;
      numCols_ = numCols;
      stride_ = stride;
      contigValues_ = values;
      origNumRows_ = numRows;
      origNumCols_ = numCols;
    }

    /// \brief Initialize a multivector that <i>is</i> a view of another multivector.
    ///
    /// \param numRows [in] Number of rows in the view.
    /// \param numCols [in] Number of columns in the view.
    /// \param values [in] Array of the multivector's entries, stored
    ///   in column-major order with stride \c stride between columns.
    ///   The array is not copied; this method just keeps a reference.
    /// \param stride [in] The column stride, that is, the number of
    ///   entries between adjacent columns of the multivector.  If you
    ///   are familiar with the BLAS or LAPACK, <tt>stride</tt> here
    ///   corresponds to "LDA" (<i>l</i>eading <i>d</i>imension of the
    ///   matrix A).
    /// \param origNumRows [in] Number of rows in the "original"
    ///   multivector (of which this multivector will be a view).
    /// \param origNumCols [in] Number of columns in the "original"
    ///   multivector (of which this multivector will be a view).
    void
    initializeValues (size_t numRows,
                      size_t numCols,
                      const Teuchos::ArrayRCP<Scalar> &values,
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

    //! The multivector data, as a nonconst array.
    Teuchos::ArrayRCP<Scalar> getValuesNonConst () const {
      return contigValues_;
    }

    //! The multivector data, as a const array.
    Teuchos::ArrayRCP<const Scalar> getValues () const {
      return contigValues_;
    }

    //! Return a nonconst view of the data in the i-th column of the multivector.
    Teuchos::ArrayRCP<Scalar>
    getValuesNonConst(size_t i) {
#ifdef HAVE_TPETRACLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    //! Return a const view of the data in the i-th column of the multivector.
    Teuchos::ArrayRCP<const Scalar>
    getValues (size_t i) const {
#ifdef HAVE_TPETRACLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    //@}
    //! @name Attribute access methods
    //@{

    //! Node accessor
    Teuchos::RCP<Node> getNode() const { return node_; }

    //! Number of rows in the multivector (view).
    size_t getNumRows() const { return numRows_; }

    //! Number of columns
    size_t getNumCols() const { return numCols_; }

    //! Stride between adjacent columns of the multivector.
    size_t getStride() const { return stride_; }

    /// \brief "Original" number of rows (of the multivector of
    ///   which <tt>*this</tt> is a view).
    ///
    /// If this multivector is <i>not</i> a view of another
    /// multivector, then this method just returns the number of
    /// rows.
    size_t getOrigNumRows() const { return origNumRows_; }

    /// \brief "Original" number of columns (of the multivector of
    ///   which <tt>*this</tt> is a view).
    ///
    /// If this multivector is <i>not</i> a view of another
    /// multivector, then this method just returns the number of
    /// columns.
    size_t getOrigNumCols() const { return origNumCols_; }

    //@}

  private:
    Teuchos::RCP<Node> node_;
    Teuchos::ArrayRCP<Scalar> contigValues_;
    size_t numRows_;
    size_t numCols_;
    size_t stride_;
    size_t origNumRows_;
    size_t origNumCols_;
  };

} // namespace KokkosClassic

#endif /* KOKKOS_MULTIVECTOR_H */
