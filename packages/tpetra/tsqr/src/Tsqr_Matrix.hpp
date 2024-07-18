// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

    //! Constructor with dimensions.
    Matrix (const ordinal_type num_rows,
            const ordinal_type num_cols) :
      nrows_ (num_rows),
      ncols_ (num_cols),
      A_ (size_t (num_rows) * size_t (num_cols))
    {}

    //! Constructor with dimensions and fill datum.
    Matrix (const ordinal_type num_rows,
            const ordinal_type num_cols,
            const non_const_value_type& value) :
      nrows_ (num_rows),
      ncols_ (num_cols),
      A_ (size_t (num_rows) * size_t (num_cols), value)
    {}

    /// \brief Copy constructor.
    ///
    /// We need an explicit copy constructor, because otherwise the
    /// default copy constructor would override the generic matrix
    /// view "copy constructor" below.
    Matrix (const Matrix& in) :
      nrows_ (in.extent(0)),
      ncols_ (in.extent(1)),
      A_ (size_t (in.extent(0)) * size_t (in.extent(1)))
    {
      MatView<ordinal_type, const_value_type> in_view
        (in.extent(0), in.extent(1), in.data(), in.stride(1));
      deep_copy (*this, in_view);
    }

    //! Default constructor (constructs an empty matrix).
    Matrix () = default;

    /// \brief "Copy constructor" from a Matrix or MatrixView.
    ///
    /// This constructor allocates a new matrix and copies the
    /// elements of the input view into the resulting new matrix.
    /// MatrixViewType must have extent(0), extent(1), data(), and
    /// stride(1) methods that match MatView's methods.
    template<class MatrixViewType>
    Matrix (const MatrixViewType& in) :
      nrows_ (in.extent(0)),
      ncols_ (in.extent(1)),
      A_ (size_t (in.extent(0)) * size_t (in.extent(1)))
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

    constexpr ordinal_type extent (const int r) const noexcept {
      return r == 0 ? nrows_ : (r == 1 ? ncols_ : ordinal_type(0));
    }

    constexpr ordinal_type stride(const int r) const noexcept {
      return r == 0 ? ordinal_type(1) : (r == 1 ? nrows_ : ordinal_type(0));
    }

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

      const size_t alloc_size = size_t (num_rows) * size_t (num_cols);
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

  template<class LO, class SC>
  bool empty (const Matrix<LO, SC>& A) {
    return A.extent(0) == 0 || A.extent(1) == 0;
  }

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

  template<class LO, class TargetScalar, class SourceScalar>
  void
  copy_upper_triangle (Matrix<LO, TargetScalar>& R_out,
                       const MatView<LO, SourceScalar>& R_in)
  {
    copy_upper_triangle (R_out.view (), R_in);
  }

  template<class LO, class TargetScalar, class SourceScalar>
  void
  copy_upper_triangle (Matrix<LO, TargetScalar>& R_out,
                       const Matrix<LO, SourceScalar>& R_in)
  {
    auto R_out_view = R_out.view ();
    copy_upper_triangle (R_out_view, R_in.const_view ());
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
