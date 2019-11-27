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
// ************************************************************************
//@HEADER

#ifndef __TSQR_Tsqr_Matrix_hpp
#define __TSQR_Tsqr_Matrix_hpp

#include "Tsqr_Util.hpp"
#include "Tsqr_MatView.hpp"
#include <stdexcept>
#include <sstream>
#include <limits>
#include <vector>

namespace TSQR {
  /// \class Matrix
  /// \brief A column-oriented dense matrix
  /// \author Mark Hoemmen
  ///
  /// A column-oriented dense matrix, with indices of type Ordinal and
  /// elements of type Scalar.
  ///
  /// \note This class resembles Teuchos::SerialDenseMatrix.  It
  ///   existed originally because there was a need for TSQR to build
  ///   independently of Teuchos.  That requirement no longer exists,
  ///   but for various reasons it has been helpful to keep Matrix
  ///   around.  In particular, I can change the interface of Matrix
  ///   without affecting other Teuchos users.
  template<class Ordinal, class Scalar>
  class Matrix {
  public:
    using non_const_value_type = typename std::remove_const<Scalar>::type;
    static_assert (std::is_same<non_const_value_type, Scalar>::value,
                   "Scalar must be nonconst.");
    using const_value_type = const non_const_value_type;
    using ordinal_type = Ordinal;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using reference = Scalar&;
    using const_reference = const Scalar&;

    using mat_view_type = MatView<ordinal_type, non_const_value_type>;
    using const_mat_view_type = MatView<ordinal_type, const_value_type>;

  private:
    static bool
    fits_in_size_t (const ordinal_type& ord)
    {
      const ordinal_type result = ordinal_type (size_t (ord));
      return (ord == result);
    }

    /// Check whether num_rows*num_cols makes sense as an amount of
    /// storage (for the num_rows by num_cols dense matrix).  Not
    /// making sense includes negative values for either parameter (if
    /// they are signed types), or overflow when computing their
    /// product.  Throw an exception of the appropriate type for any
    /// of these cases.  Otherwise, return num_rows*num_cols as a
    /// size_t.
    ///
    /// \param num_rows [in] Number of rows in the matrix
    /// \param num_cols [in] Number of columns in the matrix
    /// \return num_rows*num_cols
    size_t
    verified_alloc_size (const ordinal_type num_rows,
                         const ordinal_type num_cols) const
    {
      static_assert (std::numeric_limits<ordinal_type>::is_integer,
                     "ordinal_type must be an integer type.");
      // Quick exit also checks for zero num_cols (which prevents
      // division by zero in the tests below).
      if (num_rows == 0 || num_cols == 0) {
        return size_t(0);
      }

      // If ordinal_type is signed, make sure that num_rows and num_cols
      // are nonnegative.
      if (std::numeric_limits<ordinal_type>::is_signed) {
        if (num_rows < 0) {
          std::ostringstream os;
          os << "# rows (= " << num_rows << ") < 0";
          throw std::logic_error (os.str());
        }
        else if (num_cols < 0) {
          std::ostringstream os;
          os << "# columns (= " << num_cols << ") < 0";
          throw std::logic_error (os.str());
        }
      }

      // If ordinal_type is bigger than a size_t, do special range
      // checking.  The compiler warns (comparison of signed and
      // unsigned) if ordinal_type is a signed type and we try to do
      // "numeric_limits<size_t>::max() <
      // std::numeric_limits<ordinal_type>::max()", so instead we cast each
      // of num_rows and num_cols to size_t and back to ordinal_type again,
      // and see if we get the same result.  If not, then we
      // definitely can't return a size_t product of num_rows and
      // num_cols.
      if (! fits_in_size_t (num_rows)) {
        std::ostringstream os;
        os << "# rows (= " << num_rows << ") > max size_t value (= "
           << std::numeric_limits<size_t>::max() << ")";
        throw std::range_error (os.str());
      }
      else if (! fits_in_size_t (num_cols)) {
        std::ostringstream os;
        os << "# columns (= " << num_cols << ") > max size_t value (= "
           << std::numeric_limits<size_t>::max() << ")";
        throw std::range_error (os.str());
      }

      // Both num_rows and num_cols fit in a size_t, and are
      // nonnegative.  Now check whether their product also fits in a
      // size_t.
      //
      // Note: This may throw a SIGFPE (floating-point exception) if
      // num_cols is zero.  Be sure to check first (above).
      if (size_t (num_rows) >
          std::numeric_limits<size_t>::max() / size_t (num_cols)) {
        std::ostringstream os;
        os << "num_rows (= " << num_rows << ") * num_cols (= "
           << num_cols << ") > max size_t value (= "
           << std::numeric_limits<size_t>::max() << ")";
        throw std::range_error (os.str());
      }
      return size_t (num_rows) * size_t (num_cols);
    }

  public:
    //! Constructor with dimensions.
    Matrix (const ordinal_type num_rows,
            const ordinal_type num_cols) :
      nrows_ (num_rows),
      ncols_ (num_cols),
      A_ (verified_alloc_size (num_rows, num_cols))
    {}

    //! Constructor with dimensions and fill datum.
    Matrix (const ordinal_type num_rows,
            const ordinal_type num_cols,
            const non_const_value_type& value) :
      nrows_ (num_rows),
      ncols_ (num_cols),
      A_ (verified_alloc_size (num_rows, num_cols), value)
    {}

    /// \brief Copy constructor.
    ///
    /// We need an explicit copy constructor, because otherwise the
    /// default copy constructor would override the generic matrix
    /// view "copy constructor" below.
    Matrix (const Matrix& in) :
      nrows_ (in.extent(0)),
      ncols_ (in.extent(1)),
      A_ (verified_alloc_size (in.extent(0), in.extent(1)))
    {
      if (! in.empty()) {
        MatView<ordinal_type, non_const_value_type> this_view
          (extent(0), extent(1), data(), stride(1));
        MatView<ordinal_type, const_value_type> in_view
          (in.extent(0), in.extent(1), in.data(), in.stride(1));
        deep_copy (this_view, in_view);
      }
    }

    //! Default constructor (constructs an empty matrix).
    Matrix () = default;

    /// \brief "Copy constructor" from a matrix view type.
    ///
    /// This constructor allocates a new matrix and copies the
    /// elements of the input view into the resulting new matrix.
    /// MatrixViewType must have extent(0), extent(1), data(), and
    /// stride(1) methods that match MatView's methods.
    template<class MatrixViewType>
    Matrix (const MatrixViewType& in) :
      nrows_ (in.extent(0)),
      ncols_ (in.extent(1)),
      A_ (verified_alloc_size (in.extent(0), in.extent(1)))
    {
      if (A_.size() != 0) {
        MatView<ordinal_type, non_const_value_type> this_view
          (extent(0), extent(1), data(), stride(1));
        MatView<ordinal_type, const_value_type> in_view
          (in.extent(0), in.extent(1), in.data(), in.stride(1));
        deep_copy (this_view, in_view);
      }
    }

    /// \brief Non-const reference to element (i,j) of the matrix.
    ///
    /// \param i [in] Zero-based row index of the matrix.
    /// \param j [in] Zero-based column index of the matrix.
    reference operator() (const ordinal_type i,
                          const ordinal_type j) {
      return A_[i + j*stride(1)];
    }

    /// \brief Const reference to element (i,j) of the matrix.
    ///
    /// \param i [in] Zero-based row index of the matrix.
    /// \param j [in] Zero-based column index of the matrix.
    const_reference operator() (const ordinal_type i,
                                const ordinal_type j) const {
      return A_[i + j*stride(1)];
    }

    //! 1-D std::vector - style access.
    reference operator[] (const ordinal_type i) {
      return A_[i];
    }

    //! Equality: ONLY compares dimensions and pointers (shallow).
    template<class MatrixViewType>
    bool operator== (const MatrixViewType& B) const
    {
      if (data() != B.data() || extent(0) != B.extent(0) ||
          extent(1) != B.extent(1) || stride(1) != B.stride(1)) {
        return false;
      } else {
        return true;
      }
    }

    constexpr ordinal_type extent (const int r) const noexcept {
      return r == 0 ? nrows_ : (r == 1 ? ncols_ : ordinal_type(0));
    }

    constexpr ordinal_type stride(const int r) const noexcept {
      return r == 0 ? ordinal_type(1) : (r == 1 ? nrows_ : ordinal_type(0));
    }

    //! Whether the matrix is empty (has either zero rows or zero columns).
    bool empty() const { return extent(0) == 0 || extent(1) == 0; }

    //! A non-const pointer to the matrix data.
    pointer data()
    {
      return A_.size() != 0 ? A_.data () : nullptr;
    }

    //! A const pointer to the matrix data.
    const_pointer data() const
    {
      return A_.size() != 0 ? A_.data () : nullptr;
    }

    //! A non-const view of the matrix.
    mat_view_type view () {
      return mat_view_type (extent(0), extent(1), data(), stride(1));
    }

    //! A const view of the matrix.
    const_mat_view_type const_view () const {
      return const_mat_view_type (extent(0), extent(1),
                                  const_cast<const_pointer> (data()), stride(1));
    }

    /// Change the dimensions of the matrix.  Reallocate if necessary.
    /// Existing data in the matrix is invalidated.
    ///
    /// \param num_rows [in] New number of rows in the matrix
    /// \param num_cols [in] New number of columns in the matrix
    ///
    /// \warning This does <it>not</it> do the same thing as the
    ///   Matlab function of the same name.  In particular, it does
    ///   not reinterpret the existing matrix data using different
    ///   dimensions.
    void
    reshape (const ordinal_type num_rows, const ordinal_type num_cols)
    {
      if (num_rows == extent(0) && num_cols == extent(1))
        return; // no need to reallocate or do anything else

      const size_t alloc_size = verified_alloc_size (num_rows, num_cols);
      nrows_ = num_rows;
      ncols_ = num_cols;
      A_.resize (alloc_size);
    }

  private:
    //! Number of rows in the matrix.
    ordinal_type nrows_ = 0;
    //! Number of columns in the matrix.
    ordinal_type ncols_ = 0;
    /// \brief Where the entries of the matrix are stored.
    ///
    /// The matrix is stored using one-dimensional storage with
    /// column-major (Fortran-style) indexing.  This makes Matrix
    /// compatible with the BLAS and LAPACK.
    std::vector<non_const_value_type> A_;
  };

  template<class LO, class SC, class SourceScalar>
  void
  deep_copy (Matrix<LO, SC>& tgt, const SourceScalar& src)
  {
    deep_copy (tgt.view(), src);
  }

  template<class TargetOrdinal, class TargetScalar,
           class SourceOrdinal, class SourceScalar,
           template<class LO, class SC> class SourceMat>
  void
  deep_copy (Matrix<TargetOrdinal, TargetScalar>& tgt,
             const SourceMat<SourceOrdinal, SourceScalar>& src)
  {
    deep_copy (tgt.view(), src);
  }

  template<class LO, class SC>
  std::pair<MatView<LO, SC>, MatView<LO, SC>>
  partition_2x1 (Matrix<LO, SC>& A,
                 const typename Matrix<LO, SC>::ordinal_type nrows_top,
                 const bool b_contiguous_blocks = false)
  {
    return partition_2x1 (A.view(), nrows_top, b_contiguous_blocks);
  }

  template<class LO, class SC>
  std::pair<MatView<LO, const SC>, MatView<LO, const SC>>
  partition_2x1 (const Matrix<LO, SC>& A,
                 const typename Matrix<LO, SC>::ordinal_type nrows_top,
                 const bool b_contiguous_blocks = false)
  {
    return partition_2x1 (A.view(), nrows_top, b_contiguous_blocks);
  }
} // namespace TSQR

#endif // __TSQR_Tsqr_Matrix_hpp
