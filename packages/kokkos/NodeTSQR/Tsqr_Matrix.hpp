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

#ifndef __TSQR_Tsqr_Matrix_hpp
#define __TSQR_Tsqr_Matrix_hpp

#include <Tsqr_Util.hpp>
#include <Tsqr_MatView.hpp>

#include <stdexcept>
#include <sstream>
#include <limits>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
  private:
    static bool 
    fits_in_size_t (const Ordinal& ord) 
    {
      const Ordinal result = static_cast< Ordinal > (static_cast< size_t > (ord));
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
    verified_alloc_size (const Ordinal num_rows,
			 const Ordinal num_cols) const 
    {
      if (! std::numeric_limits< Ordinal >::is_integer)
	throw std::logic_error("Ordinal must be an integer type");

      // Quick exit also checks for zero num_cols (which prevents
      // division by zero in the tests below).
      if (num_rows == 0 || num_cols == 0)
	return size_t(0);

      // If Ordinal is signed, make sure that num_rows and num_cols
      // are nonnegative.
      if (std::numeric_limits< Ordinal >::is_signed)
	{
	  if (num_rows < 0)
	    {
	      std::ostringstream os;
	      os << "# rows (= " << num_rows << ") < 0";
	      throw std::logic_error (os.str());
	    }
	  else if (num_cols < 0)
	    {
	      std::ostringstream os;
	      os << "# columns (= " << num_cols << ") < 0";
	      throw std::logic_error (os.str());
	    }
	}

      // If Ordinal is bigger than a size_t, do special range
      // checking.  The compiler warns (comparison of signed and
      // unsigned) if Ordinal is a signed type and we try to do
      // "numeric_limits<size_t>::max() <
      // std::numeric_limits<Ordinal>::max()", so instead we cast each
      // of num_rows and num_cols to size_t and back to Ordinal again,
      // and see if we get the same result.  If not, then we
      // definitely can't return a size_t product of num_rows and
      // num_cols.
      if (! fits_in_size_t (num_rows))
	{
	  std::ostringstream os;
	  os << "# rows (= " << num_rows << ") > max size_t value (= " 
	     << std::numeric_limits<size_t>::max() << ")";
	  throw std::range_error (os.str());
	}
      else if (! fits_in_size_t (num_cols))
	{
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
      if (static_cast<size_t>(num_rows) > 
	  std::numeric_limits<size_t>::max() / static_cast<size_t>(num_cols))
	{
	  std::ostringstream os;
	  os << "num_rows (= " << num_rows << ") * num_cols (= "
	     << num_cols << ") > max size_t value (= " 
	     << std::numeric_limits<size_t>::max() << ")";
	  throw std::range_error (os.str());
	}
      return static_cast<size_t>(num_rows) * static_cast<size_t>(num_cols);
    }
    
  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef Scalar* pointer_type;

    //! Constructor with dimensions.
    Matrix (const Ordinal num_rows, 
	    const Ordinal num_cols) :
      nrows_ (num_rows),
      ncols_ (num_cols),
      A_ (verified_alloc_size (num_rows, num_cols))
    {}

    //! Constructor with dimensions and fill datum.
    Matrix (const Ordinal num_rows,
	    const Ordinal num_cols,
	    const Scalar& value) :
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
      nrows_ (in.nrows()),
      ncols_ (in.ncols()),
      A_ (verified_alloc_size (in.nrows(), in.ncols()))
    {
      if (! in.empty())
	copy_matrix (nrows(), ncols(), get(), lda(), in.get(), in.lda());
    }

    //! Default constructor (constructs an empty matrix).
    Matrix () : nrows_(0), ncols_(0), A_(0) {}

    //! Trivial destructor.
    ~Matrix () {} 

    /// \brief "Copy constructor" from a matrix view type.
    ///
    /// This constructor allocates a new matrix and copies the
    /// elements of the input view into the resulting new matrix.
    /// MatrixViewType must have nrows(), ncols(), get(), and lda()
    /// methods that match MatView's methods.
    template<class MatrixViewType>
    Matrix (const MatrixViewType& in) :
      nrows_ (in.nrows()),
      ncols_ (in.ncols()),
      A_ (verified_alloc_size (in.nrows(), in.ncols()))
    {
      if (A_.size() != 0)
	copy_matrix (nrows(), ncols(), get(), lda(), in.get(), in.lda());
    }

    /// \brief Copy the matrix B into this matrix.
    ///
    /// \param B [in] A matrix with the same dimensions (numbers of
    ///   rows and columns) as this matrix.
    template<class MatrixViewType>
    void 
    copy (MatrixViewType& B)
    {
      const typename MatrixViewType::ordinal_type num_rows = B.nrows();
      const typename MatrixViewType::ordinal_type num_cols = B.ncols();
      if (num_rows != nrows() || num_cols != ncols())
	{
	  std::ostringstream os;
	  os << "Matrix::Copy: incompatible dimensions: attempt to assign " 
	     << num_rows << " x " << num_cols << " matrix to an " 
	     << (nrows()) << " x " << (ncols()) << "matrix";
	  throw std::logic_error (os.str());
	}
      copy_matrix (nrows(), ncols(), get(), lda(), B.get(), B.lda());
    }

    //! Fill all entries of the matrix with the given value.
    void 
    fill (const Scalar value)
    {
      fill_matrix (nrows(), ncols(), get(), lda(), value);
    }

    /// \brief Non-const reference to element (i,j) of the matrix.
    ///
    /// \param i [in] Zero-based row index of the matrix.
    /// \param j [in] Zero-based column index of the matrix.
    Scalar& operator() (const Ordinal i, const Ordinal j) {
      return A_[i + j*lda()];
    }

    /// \brief Const reference to element (i,j) of the matrix.
    ///
    /// \param i [in] Zero-based row index of the matrix.
    /// \param j [in] Zero-based column index of the matrix.
    const Scalar& operator() (const Ordinal i, const Ordinal j) const {
      return A_[i + j*lda()];
    }

    //! 1-D std::vector - style access.
    Scalar& operator[] (const Ordinal i) {
      return A_[i];
    }

    /// \brief Equality test: compares dimensions and entries.
    ///
    /// The test is templated so that B may have a different Ordinal
    /// type or even a different Scalar type than *this.
    template<class MatrixViewType>
    bool operator== (const MatrixViewType& B) const 
    {
      if (nrows() != B.nrows() || ncols() != B.ncols())
	return false;
    
      typedef typename MatrixViewType::ordinal_type second_ordinal_type;
      typedef typename MatrixViewType::scalar_type second_scalar_type;
      typedef typename MatrixViewType::pointer_type second_pointer_type;

      const ordinal_type A_nrows = nrows();
      const ordinal_type A_lda = lda();
      const ordinal_type A_ncols = ncols();
      const second_ordinal_type B_lda = B.lda();
      const scalar_type* A_j = get();
      const second_scalar_type* B_j = B.get();

      for (ordinal_type j = 0; j < A_ncols; ++j, A_j += A_lda, B_j += B_lda)
	for (ordinal_type i = 0; i < A_nrows; ++i)
	  if (A_j[i] != B_j[i])
	    return false;
      return true;
    }

    //! Number of rows in the matrix.
    Ordinal nrows() const { return nrows_; }

    //! Number of columns in the matrix.
    Ordinal ncols() const { return ncols_; }

    //! Leading dimension (a.k.a. stride) of the matrix.
    Ordinal lda() const { return nrows_; }

    //! Whether the matrix is empty (has either zero rows or zero columns).
    bool empty() const { return nrows() == 0 || ncols() == 0; }

    //! A non-const pointer to the matrix data.
    Scalar* 
    get() 
    { 
      if (A_.size() > 0)
	return &A_[0]; 
      else
	return static_cast<Scalar*> (NULL);
    }

    //! A const pointer to the matrix data.
    const Scalar* 
    get() const
    { 
      if (A_.size() > 0)
	return &A_[0]; 
      else
	return static_cast<const Scalar*> (NULL);
    }

    //! A non-const view of the matrix.
    MatView<Ordinal, Scalar> view () {
      return MatView<Ordinal, Scalar> (nrows(), ncols(), get(), lda());
    }

    //! A const view of the matrix.
    ConstMatView<Ordinal, Scalar> const_view () const {
      typedef ConstMatView<Ordinal, Scalar> const_view_type;
      return const_view_type (nrows(), ncols(), 
			      const_cast<const Scalar*> (get()), lda());
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
    reshape (const Ordinal num_rows, const Ordinal num_cols)
    {
      if (num_rows == nrows() && num_cols == ncols())
	return; // no need to reallocate or do anything else

      const size_t alloc_size = verified_alloc_size (num_rows, num_cols);
      nrows_ = num_rows;
      ncols_ = num_cols;
      A_.resize (alloc_size);
    }

  private:
    //! Number of rows in the matrix.
    Ordinal nrows_;
    //! Number of columns in the matrix.
    Ordinal ncols_;
    /// \brief Where the entries of the matrix are stored.
    ///
    /// The matrix is stored using one-dimensional storage with
    /// column-major (Fortran-style) indexing.  This makes Matrix
    /// compatible with the BLAS and LAPACK.
    std::vector<Scalar> A_;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_Matrix_hpp
