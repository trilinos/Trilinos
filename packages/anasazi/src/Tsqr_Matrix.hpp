#ifndef __TSQR_Tsqr_Matrix_hpp
#define __TSQR_Tsqr_Matrix_hpp

#include <Tsqr_Util.hpp>
#include <Tsqr_MatView.hpp>

#include <stdexcept>
#include <sstream>
#include <limits>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class Matrix
  /// \brief A column-oriented matrix.
  ///
  /// \note This class is a bit like Teuchos::SerialDenseMatrix.  It
  ///   exists because I didn't want TSQR to depend on Teuchos.  This
  ///   is not because I don't like Teuchos, but because I had/have
  ///   potential TSQR users who might not want to depend on Teuchos.
  template< class Ordinal, class Scalar >
  class Matrix {
  private:
    size_t
    verified_alloc_size (const Ordinal num_rows,
			 const Ordinal num_cols) const 
    {
      if (! std::numeric_limits< Ordinal >::is_integer)
	throw std::logic_error("Ordinal must be an integer type");

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

      // If Ordinal is bigger than a size_t, do special range checking.
      if (std::numeric_limits<size_t>::max() < std::numeric_limits< Ordinal >::max())
	{
	  if (num_rows > std::numeric_limits<size_t>::max())
	    {
	      std::ostringstream os;
	      os << "# rows (= " << num_rows << ") > max size_t value (= " 
		 << std::numeric_limits<size_t>::max() << ")";
	      throw std::range_error (os.str());
	    }
	  else if (num_cols > std::numeric_limits<size_t>::max())
	    {
	      std::ostringstream os;
	      os << "# columns (= " << num_cols << ") > max size_t value (= "
		 << std::numeric_limits<size_t>::max() << ")";
	      throw std::range_error (os.str());
	    }
	}

      // Both num_rows and num_cols <= max size_t value.
      // Now make sure that their product is also <= max size_t value.
      if (num_rows > std::numeric_limits<size_t>::max() / num_cols)
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
    typedef Scalar value_type; // deprecated
    typedef Ordinal ordinal_type;
    typedef Ordinal index_type; // deprecated
    typedef Scalar* pointer_type;

    Matrix (const Ordinal num_rows, 
	    const Ordinal num_cols)
    {
      const size_t alloc_size = verified_alloc_size (num_rows, num_cols);
      nrows_ = num_rows;
      ncols_ = num_cols;
      if (alloc_size == 0)
	A_ = NULL;
      else
	A_ = new double[alloc_size];
    }

    Matrix (const Ordinal num_rows,
	    const Ordinal num_cols,
	    const Scalar value) 
    {
      const size_t alloc_size = verified_alloc_size (num_rows, num_cols);
      nrows_ = num_rows;
      ncols_ = num_cols;
      if (alloc_size == 0)
	A_ = NULL;
      else
	{
	  A_ = new double[alloc_size];
	  fill_matrix (num_rows, num_cols, get(), lda(), value);
	}
    }

    // We need an explicit copy constructor, because for some reason
    // the default copy constructor (with shallow copies of pointers,
    // eeek! double free()s!!!) overrides the generic "copy
    // constructors" below.
    Matrix (const Matrix& in)
    {
      const size_t alloc_size = verified_alloc_size (in.nrows(), in.ncols());
      nrows_ = in.nrows();
      ncols_ = in.ncols();
      if (alloc_size == 0)
	A_ = NULL;
      else 
	{
	  A_ = new double[alloc_size];
	  copy_matrix (nrows(), ncols(), get(), lda(), in.get(), in.lda());
	}
    }

    Matrix () : nrows_(0), ncols_(0), A_(NULL) {}

    ~Matrix () { 
      if (A_ != NULL) 
	delete [] A_; 
      A_ = NULL;
    }

    template< class MatrixViewType >
    Matrix (const MatrixViewType& in)
    {
      const size_t alloc_size = verified_alloc_size (in.nrows(), in.ncols());
      nrows_ = in.nrows();
      ncols_ = in.ncols();
      if (alloc_size == 0)
	A_ = NULL;
      else 
	{
	  A_ = new double[alloc_size];
	  copy_matrix (nrows(), ncols(), get(), lda(), in.get(), in.lda());
	}
    }

    /// *this gets a deep copy of B.
    ///
    /// \note Assumes *this and B have the same dimensions (but not
    ///   necessarily the same strides).
    template< class MatrixViewType >
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

    void 
    fill (const Scalar value)
    {
      fill_matrix (nrows(), ncols(), get(), lda(), value);
    }

    /// 2-D Fortran array - style access (column-major order, though
    /// indexing is zero-based (C style) instead of one-based (Fortran
    /// style))
    Scalar& operator() (const Ordinal i, const Ordinal j) {
      return A_[i + j*lda()];
    }

    /// 1-D std::vector - style access
    Scalar& operator[] (const Ordinal i) {
      return A_[i];
    }

    template< class MatrixViewType >
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

    Ordinal nrows() const { return nrows_; }
    Ordinal ncols() const { return ncols_; }
    Ordinal lda() const { return nrows_; }
    bool empty() const { return nrows() == 0 || ncols() == 0; }
    Scalar* get() { return A_; }
    const Scalar* get() const { return A_; }

    MatView< Ordinal, Scalar > view () {
      return MatView< Ordinal, Scalar >(nrows(), ncols(), get(), lda());
    }

    /// Change the dimensions of the matrix.  Reallocate if necessary.
    /// Existing data in the matrix is invalidated.
    ///
    /// \param num_rows [in] New number of rows in the matrix
    /// \param num_cols [in] New number of columns in the matrix
    void
    reshape (const Ordinal num_rows, const Ordinal num_cols)
    {
      if (A_ != NULL)
	delete A_;

      const size_t alloc_size = verified_alloc_size (num_rows, num_cols);
      nrows_ = num_rows;
      ncols_ = num_cols;
      if (alloc_size == 0)
	A_ = NULL;
      else
	A_ = new double[alloc_size];
    }

  private:
    Ordinal nrows_, ncols_;
    Scalar* A_;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_Matrix_hpp
