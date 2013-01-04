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
    MultiVector (RCP<Node> node)
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
                      const ArrayRCP<Scalar> &values,
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

    //! The multivector data, as a nonconst array.
    ArrayRCP<Scalar> getValuesNonConst () {
      return contigValues_;
    }

    //! The multivector data, as a const array.
    ArrayRCP<const Scalar> getValues () const {
      return contigValues_;
    }

    //! Return a nonconst view of the data in the i-th column of the multivector.
    ArrayRCP<Scalar>
    getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    //! Return a const view of the data in the i-th column of the multivector.
    ArrayRCP<const Scalar>
    getValues (size_t i) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    //@}
    //! @name Element update / replace methods (VERY SLOW on GPU Nodes)
    //@{

    /// \brief <tt>X(i,j) = X(i,j) + newVal</tt>, where X is this multivector.
    ///
    /// \warning This method will be VERY SLOW if Node is a GPU node.
    ///   That's because it has to copy the one entry from the GPU to
    ///   the CPU, update it on the CPU, and then copy it back to the
    ///   GPU.
    ///
    /// \param i [in] (Local) row index of the entry to update.
    /// \param j [in] Column index of the entry to update.
    /// \param newVal [in] The value to use when updating the (i,j)
    ///   entry of this multivector.
    void
    sumIntoLocalValue (size_t i, size_t j, const Scalar& newVal) {
      using Teuchos::ArrayRCP;

      const size_t offset = i + j * getStride ();
      const size_t numEntries = 1; // We're only getting one entry.

      // Device view of the one entry.
      Teuchos::ArrayRCP<Scalar> deviceView =
        contigValues_.persistingView (offset, numEntries);
      // Get a read-write host view of the one entry.
      Teuchos::ArrayRCP<Scalar> hostView =
        this->getNode ()->template viewBufferNonConst<Scalar> (Kokkos::ReadWrite,
                                                               numEntries, deviceView);
      // Update the entry.  When hostView falls out of scope at the
      // end of this method, it will copy the updated entry back to
      // the device.
      hostView[0] += newVal;
    }

    /// \brief <tt>X(i,j) = newVal</tt>, where X is this multivector.
    ///
    /// \warning This method will be VERY SLOW if Node is a GPU node.
    ///   That's because it has to copy the new value of the entry
    ///   from the CPU to the GPU.
    ///
    /// \param i [in] (Local) row index of the entry to replace.
    /// \param j [in] Column index of the entry to replace.
    /// \param newVal [in] The value to use when replacing the (i,j)
    ///   entry of this multivector.
    void
    replaceLocalValue (size_t i, size_t j, const Scalar& newVal) {
      const size_t offset = i + j * getStride ();
      const size_t numEntries = 1; // We're only getting one entry.

      // Device view of the one entry.
      Teuchos::ArrayRCP<Scalar> deviceView =
        contigValues_.persistingView (offset, numEntries);
      // Get a host view of the one entry.  A write-only view
      //  suffices, since we're just replacing the value.
      Teuchos::ArrayRCP<Scalar> hostView =
        this->getNode ()->template viewBufferNonConst<Scalar> (Kokkos::WriteOnly,
                                                               numEntries, deviceView);
      // Update the entry.  When hostView falls out of scope at the
      // end of this method, it will copy the value back to the
      // device.
      hostView[0] = newVal;
    }

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

      const size_t origNumRows = this->getOrigNumRows ();
      const size_t origNumCols = this->getOrigNumCols ();
      // mfh 04 Jan 2013: Thanks to Jonathan Hu for pointing out the
      // necessary test to allow a view of a 0 x C or R x 0
      // multivector when the offset(s) corresponding to the zero
      // dimension(s) is/are also zero.
      TEUCHOS_TEST_FOR_EXCEPTION(
        (origNumRows == 0 && offsetRow > 0) ||
        (origNumRows > 0 && offsetRow >= origNumRows) ||
        (origNumCols == 0 && offsetCol > 0) ||
        (origNumCols > 0 && offsetCol >= origNumCols),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: offset row or column are "
        "out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: new dimensions are out "
        "of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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
        Teuchos::typeName (*this) << "::offsetViewNonConst: offset row or "
        "column are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: new dimensions "
        "are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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
    //! @name Attribute access methods
    //@{

    //! Node accessor
    RCP<Node> getNode() const { return node_; }

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
    RCP<Node> node_;
    ArrayRCP<Scalar> contigValues_;
    size_t numRows_;
    size_t numCols_;
    size_t stride_;
    size_t origNumRows_;
    size_t origNumCols_;
  };

  //
  // Partial specializations for various CPU Node types.  The only
  // methods that are different thus far are sumIntoLocalValue and
  // replaceLocalValue.  The specializations will probably make the
  // methods faster on CPU Nodes.
  //

  // Partial specialization for SerialNode.
  template<class Scalar>
  class MultiVector<Scalar, SerialNode> {
  public:
    typedef Scalar ScalarType;
    typedef SerialNode NodeType;

    MultiVector (RCP<SerialNode> node)
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0)
      , origNumRows_(0)
      , origNumCols_(0) {
    }

    MultiVector (const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_)
      , origNumRows_(source.origNumRows_)
      , origNumCols_(source.origNumCols_) {
    }

    ~MultiVector() {
    }

    void
    initializeValues (size_t numRows,
                      size_t numCols,
                      const ArrayRCP<Scalar> &values,
                      size_t stride)
    {
      numRows_ = numRows;
      numCols_ = numCols;
      stride_ = stride;
      contigValues_ = values;
      origNumRows_ = numRows;
      origNumCols_ = numCols;
    }

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

    ArrayRCP<Scalar> getValuesNonConst () {
      return contigValues_;
    }

    ArrayRCP<const Scalar> getValues () const {
      return contigValues_;
    }

    ArrayRCP<Scalar>
    getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    ArrayRCP<const Scalar>
    getValues (size_t i) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    void
    sumIntoLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] += newVal;
    }

    void
    replaceLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] = newVal;
    }

    const MultiVector<Scalar,SerialNode>
    offsetView (size_t newNumRows,
                size_t newNumCols,
                size_t offsetRow,
                size_t offsetCol) const
    {
      MultiVector<Scalar,SerialNode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: offset row or column are "
        "out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: new dimensions are out "
        "of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    MultiVector<Scalar,SerialNode>
    offsetViewNonConst (size_t newNumRows,
                        size_t newNumCols,
                        size_t offsetRow,
                        size_t offsetCol)
    {
      MultiVector<Scalar,SerialNode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: offset row or "
        "column are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: new dimensions "
        "are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    RCP<SerialNode> getNode() const { return node_; }
    size_t getNumRows() const { return numRows_; }
    size_t getNumCols() const { return numCols_; }
    size_t getStride() const { return stride_; }
    size_t getOrigNumRows() const { return origNumRows_; }
    size_t getOrigNumCols() const { return origNumCols_; }

  private:
    RCP<SerialNode> node_;
    ArrayRCP<Scalar> contigValues_;
    size_t numRows_;
    size_t numCols_;
    size_t stride_;
    size_t origNumRows_;
    size_t origNumCols_;
  };

#if defined (HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)
  // Partial specialization for TPINode.
  template<class Scalar>
  class MultiVector<Scalar, TPINode> {
  public:
    typedef Scalar ScalarType;
    typedef TPINode NodeType;

    MultiVector (RCP<TPINode> node)
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0)
      , origNumRows_(0)
      , origNumCols_(0) {
    }

    MultiVector (const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_)
      , origNumRows_(source.origNumRows_)
      , origNumCols_(source.origNumCols_) {
    }

    ~MultiVector() {
    }

    void
    initializeValues (size_t numRows,
                      size_t numCols,
                      const ArrayRCP<Scalar> &values,
                      size_t stride)
    {
      numRows_ = numRows;
      numCols_ = numCols;
      stride_ = stride;
      contigValues_ = values;
      origNumRows_ = numRows;
      origNumCols_ = numCols;
    }

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

    ArrayRCP<Scalar> getValuesNonConst () {
      return contigValues_;
    }

    ArrayRCP<const Scalar> getValues () const {
      return contigValues_;
    }

    ArrayRCP<Scalar>
    getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    ArrayRCP<const Scalar>
    getValues (size_t i) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    void
    sumIntoLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] += newVal;
    }

    void
    replaceLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] = newVal;
    }

    const MultiVector<Scalar,TPINode>
    offsetView (size_t newNumRows,
                size_t newNumCols,
                size_t offsetRow,
                size_t offsetCol) const
    {
      MultiVector<Scalar,TPINode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: offset row or column are "
        "out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: new dimensions are out "
        "of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    MultiVector<Scalar,TPINode>
    offsetViewNonConst (size_t newNumRows,
                        size_t newNumCols,
                        size_t offsetRow,
                        size_t offsetCol)
    {
      MultiVector<Scalar,TPINode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: offset row or "
        "column are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: new dimensions "
        "are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    RCP<TPINode> getNode() const { return node_; }
    size_t getNumRows() const { return numRows_; }
    size_t getNumCols() const { return numCols_; }
    size_t getStride() const { return stride_; }
    size_t getOrigNumRows() const { return origNumRows_; }
    size_t getOrigNumCols() const { return origNumCols_; }

  private:
    RCP<TPINode> node_;
    ArrayRCP<Scalar> contigValues_;
    size_t numRows_;
    size_t numCols_;
    size_t stride_;
    size_t origNumRows_;
    size_t origNumCols_;
  };
#endif // defined (HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)


#if defined (HAVE_KOKKOSCLASSIC_DEFAULTNODE_TBBNODE)
  // Partial specialization for TBBNode.
  template<class Scalar>
  class MultiVector<Scalar, TBBNode> {
  public:
    typedef Scalar ScalarType;
    typedef TBBNode NodeType;

    MultiVector (RCP<TBBNode> node)
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0)
      , origNumRows_(0)
      , origNumCols_(0) {
    }

    MultiVector (const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_)
      , origNumRows_(source.origNumRows_)
      , origNumCols_(source.origNumCols_) {
    }

    ~MultiVector() {
    }

    void
    initializeValues (size_t numRows,
                      size_t numCols,
                      const ArrayRCP<Scalar> &values,
                      size_t stride)
    {
      numRows_ = numRows;
      numCols_ = numCols;
      stride_ = stride;
      contigValues_ = values;
      origNumRows_ = numRows;
      origNumCols_ = numCols;
    }

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

    ArrayRCP<Scalar> getValuesNonConst () {
      return contigValues_;
    }

    ArrayRCP<const Scalar> getValues () const {
      return contigValues_;
    }

    ArrayRCP<Scalar>
    getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    ArrayRCP<const Scalar>
    getValues (size_t i) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    void
    sumIntoLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] += newVal;
    }

    void
    replaceLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] = newVal;
    }

    const MultiVector<Scalar,TBBNode>
    offsetView (size_t newNumRows,
                size_t newNumCols,
                size_t offsetRow,
                size_t offsetCol) const
    {
      MultiVector<Scalar,TBBNode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: offset row or column are "
        "out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: new dimensions are out "
        "of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    MultiVector<Scalar,TBBNode>
    offsetViewNonConst (size_t newNumRows,
                        size_t newNumCols,
                        size_t offsetRow,
                        size_t offsetCol)
    {
      MultiVector<Scalar,TBBNode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: offset row or "
        "column are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: new dimensions "
        "are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    RCP<TBBNode> getNode() const { return node_; }
    size_t getNumRows() const { return numRows_; }
    size_t getNumCols() const { return numCols_; }
    size_t getStride() const { return stride_; }
    size_t getOrigNumRows() const { return origNumRows_; }
    size_t getOrigNumCols() const { return origNumCols_; }

  private:
    RCP<TBBNode> node_;
    ArrayRCP<Scalar> contigValues_;
    size_t numRows_;
    size_t numCols_;
    size_t stride_;
    size_t origNumRows_;
    size_t origNumCols_;
  };
#endif // defined (HAVE_KOKKOSCLASSIC_DEFAULTNODE_TBBNODE)


#if defined (HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPNODE)
  // Partial specialization for OpenMPNode.
  template<class Scalar>
  class MultiVector<Scalar, OpenMPNode> {
  public:
    typedef Scalar ScalarType;
    typedef OpenMPNode NodeType;

    MultiVector (RCP<OpenMPNode> node)
      : node_(node)
      , numRows_(0)
      , numCols_(0)
      , stride_(0)
      , origNumRows_(0)
      , origNumCols_(0) {
    }

    MultiVector (const MultiVector& source)
      : node_(source.node_)
      , contigValues_(source.contigValues_)
      , numRows_(source.numRows_)
      , numCols_(source.numCols_)
      , stride_(source.stride_)
      , origNumRows_(source.origNumRows_)
      , origNumCols_(source.origNumCols_) {
    }

    ~MultiVector() {
    }

    void
    initializeValues (size_t numRows,
                      size_t numCols,
                      const ArrayRCP<Scalar> &values,
                      size_t stride)
    {
      numRows_ = numRows;
      numCols_ = numCols;
      stride_ = stride;
      contigValues_ = values;
      origNumRows_ = numRows;
      origNumCols_ = numCols;
    }

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

    ArrayRCP<Scalar> getValuesNonConst () {
      return contigValues_;
    }

    ArrayRCP<const Scalar> getValues () const {
      return contigValues_;
    }

    ArrayRCP<Scalar>
    getValuesNonConst(size_t i) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValuesNonConst(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    ArrayRCP<const Scalar>
    getValues (size_t i) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      const bool inRange = i < numCols_; // i >= 0 since it's unsigned.
      TEUCHOS_TEST_FOR_EXCEPTION(
        contigValues_.is_null () || ! inRange,
        std::runtime_error,
        Teuchos::typeName(*this) << "::getValues(): index out of range or data structure not initialized.");
#endif
      return contigValues_.persistingView (stride_*i, numRows_);
    }

    void
    sumIntoLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] += newVal;
    }

    void
    replaceLocalValue (size_t i, size_t j, const Scalar& newVal) {
      contigValues_[i + j*getStride()] = newVal;
    }

    const MultiVector<Scalar,OpenMPNode>
    offsetView (size_t newNumRows,
                size_t newNumCols,
                size_t offsetRow,
                size_t offsetCol) const
    {
      MultiVector<Scalar,OpenMPNode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: offset row or column are "
        "out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetView: new dimensions are out "
        "of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    MultiVector<Scalar,OpenMPNode>
    offsetViewNonConst (size_t newNumRows,
                        size_t newNumCols,
                        size_t offsetRow,
                        size_t offsetCol)
    {
      MultiVector<Scalar,OpenMPNode> B (this->getNode ());

      TEUCHOS_TEST_FOR_EXCEPTION(
        offsetRow >= this->getOrigNumRows () || offsetCol >= this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: offset row or "
        "column are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested offset row and column are "
        << offsetRow << ", " << offsetCol << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows > this->getOrigNumRows () || newNumCols > this->getOrigNumCols (),
        std::invalid_argument,
        Teuchos::typeName (*this) << "::offsetViewNonConst: new dimensions "
        "are out of bounds.  The original multivector has dimensions "
        << this->getOrigNumRows () << " x " << this->getOrigNumCols ()
        << ", but your requested new dimensions are "
        << newNumRows << " x " << newNumCols << ".");

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

    RCP<OpenMPNode> getNode() const { return node_; }
    size_t getNumRows() const { return numRows_; }
    size_t getNumCols() const { return numCols_; }
    size_t getStride() const { return stride_; }
    size_t getOrigNumRows() const { return origNumRows_; }
    size_t getOrigNumCols() const { return origNumCols_; }

  private:
    RCP<OpenMPNode> node_;
    ArrayRCP<Scalar> contigValues_;
    size_t numRows_;
    size_t numCols_;
    size_t stride_;
    size_t origNumRows_;
    size_t origNumCols_;
  };
#endif // defined (HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPNODE)

} // namespace Kokkos

#endif /* KOKKOS_MULTIVECTOR_H */
