#ifndef __TSQR_Tsqr_MatView_hpp
#define __TSQR_Tsqr_MatView_hpp

#include <cstring> // NULL

// Define for bounds checking and other safety features, undefine for speed.
// #define TSQR_MATVIEW_DEBUG 1

#ifdef TSQR_MATVIEW_DEBUG
#  include <limits>
#endif // TSQR_MATVIEW_DEBUG

#include <sstream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class FirstMatrixViewType, class SecondMatrixViewType >
  bool
  matrix_equal (FirstMatrixViewType& A, SecondMatrixViewType& B)
  {
    if (A.nrows() != B.nrows() || A.ncols() != B.ncols())
      return false;
    
    typedef typename FirstMatrixViewType::index_type first_index_type;
    typedef typename SecondMatrixViewType::index_type second_index_type;
    typedef typename FirstMatrixViewType::pointer_type first_pointer_type;
    typedef typename SecondMatrixViewType::pointer_type second_pointer_type;

    const first_index_type nrows = A.nrows();
    const first_index_type A_lda = A.lda();
    const first_index_type ncols = A.ncols();
    const second_index_type B_lda = B.lda();

    first_pointer_type A_j = A.get();
    second_pointer_type B_j = B.get();

    for (first_index_type j = 0; j < ncols; ++j, A_j += A_lda, B_j += B_lda)
      for (first_index_type i = 0; i < nrows; ++i)
	if (A_j[i] != B_j[i])
	  return false;

    return true;
  }

#ifdef TSQR_MATVIEW_DEBUG
  template< class Ordinal, class Scalar >
  class MatViewVerify {
  public:
    static void 
    verify (const Ordinal num_rows, 
	    const Ordinal num_cols, 
	    const Scalar* const A, 
	    const Ordinal leading_dim)
    {
      using std::endl;

      bool good = true;
      std::ostringstream os;
      if (! std::numeric_limits<Ordinal>::is_integer)
	{
	  good = false;
	  os << "Error: Ordinal type must be an integer." << endl;
	}
      if (std::numeric_limits<Ordinal>::is_signed)
	{
	  if (num_rows < 0)
	    {
	      good = false;
	      os << "Error: num_rows (= " << num_rows << ") < 0." << endl;
	    }
	  if (num_cols < 0)
	    {
	      good = false;
	      os << "Error: num_cols (= " << num_cols << ") < 0." << endl;
	    }
	  if (leading_dim < 0)
	    {
	      good = false;
	      os << "Error: leading_dim (= " << leading_dim << ") < 0." << endl;
	    }
	}
      if (leading_dim < num_rows)
	{
	  good = false;
	  os << "Error: leading_dim (= " << leading_dim << ") < num_rows (= " << num_rows << ")." << endl;
	}
      if (! good)
	throw std::invalid_argument (os.str());
    }
  };
#endif // TSQR_MATVIEW_DEBUG


  /// \class MatView
  /// 
  /// A read-and-write view of a column-oriented matrix.
  template< class Ordinal, class Scalar >
  class MatView {
  public:
    typedef Scalar value_type;
    typedef Ordinal index_type;
    typedef Scalar* pointer_type;

    MatView () : nrows_(0), ncols_(0), A_(NULL), lda_(0) {}

    MatView (const Ordinal num_rows, 
	     const Ordinal num_cols, 
	     Scalar* const A, 
	     const Ordinal leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      A_(A),
      lda_(leading_dim) 
    {
#ifdef TSQR_MATVIEW_DEBUG
      MatViewVerify< Ordinal, Scalar >::verify (num_rows, num_cols, A, leading_dim);
#endif // TSQR_MATVIEW_DEBUG
    }

    MatView (const MatView& view) :
      nrows_(view.nrows()),
      ncols_(view.ncols()),
      A_(view.get()),
      lda_(view.lda())
    {}

    MatView& operator= (const MatView& view) {
      if (this != &view)
	{
	  nrows_ = view.nrows();
	  ncols_ = view.ncols();
	  A_ = view.get();
	  lda_ = view.lda();
	}
      return *this;
    }

    /// \note The function is const, only because returning a
    /// reference to the matrix data doesn't change any members of
    /// *this.  Of course one may use the resulting reference to
    /// change an entry in the matrix, but that doesn't affect the
    /// MatView's properties.
    Scalar& operator() (const Ordinal i, const Ordinal j) const 
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed)
	{
	  if (i < 0 || i >= nrows())
	    throw std::invalid_argument("Row range invalid");
	  else if (j < 0 || j >= ncols())
	    throw std::invalid_argument("Column range invalid");
	}
      else
	{
	  if (i >= nrows())
	    throw std::invalid_argument("Row range invalid");
	  else if (j >= ncols())
	    throw std::invalid_argument("Column range invalid");
	}
      if (A_ == NULL)
	throw std::logic_error("Attempt to reference NULL data");
#endif // TSQR_MATVIEW_DEBUG
      return A_[i + j*lda()];
    }

    Ordinal nrows() const { return nrows_; }
    Ordinal ncols() const { return ncols_; }
    Ordinal lda() const { return lda_; }

    /// \note The function is const, only because returning A_ doesn't
    /// change any members of *this.  Of course one may use the
    /// resulting pointer to fiddle with entries in the matrix, but
    /// that doesn't affect the MatView's properties.
    pointer_type get() const { return A_; }
    bool empty() const { return nrows() == 0 || ncols() == 0; }

    /// Return a "row block" (submatrix of consecutive rows in the
    /// inclusive range [firstRow,lastRow]).
    MatView row_block (const Ordinal firstRow, const Ordinal lastRow) 
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed)
	{
	  if (firstRow < 0 || firstRow > lastRow || lastRow >= nrows())
	    throw std::invalid_argument ("Row range invalid");
	}
      else
	{
	  if (firstRow > lastRow || lastRow >= nrows())
	    throw std::invalid_argument ("Row range invalid");
	}
#endif // TSQR_MATVIEW_DEBUG
      return MatView (lastRow - firstRow + 1, ncols(), get() + firstRow, lda());
    }

    /// Split off and return the top cache block of nrows_top rows.
    /// Modify *this to be the "rest" of the matrix.
    ///
    /// \note Only use this method to split off a single cache block.
    ///   It breaks if you try to use it otherwise.
    ///
    /// \param nrows_top [in] Number of rows in the top block (which
    ///   this method returns)
    ///
    /// \param b_contiguous_blocks [in] Whether or not the entries of
    ///   the top block are stored contiguously in *this.  The default
    ///   is no (false).
    ///
    /// \return The top block of nrows_top rows.  Data is a shallow
    ///   copy of the data in *this.
    MatView split_top (const Ordinal nrows_top, 
		       const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed && nrows_top < 0)
	{
	  std::ostringstream os;
	  os << "nrows_top (= " << nrows_top << ") < 0";
	  throw std::invalid_argument (os.str());
	}
      else if (nrows_top > nrows())
	{
	  std::ostringstream os;
	  os << "nrows_top (= " << nrows_top << ") > nrows (= " << nrows() << ")";
	  throw std::invalid_argument (os.str());
	}
#endif // TSQR_MATVIEW_DEBUG

      Scalar* const A_top_ptr = get();
      Scalar* A_rest_ptr;
      const Ordinal nrows_rest = nrows() - nrows_top;
      Ordinal lda_top, lda_rest;
      if (b_contiguous_blocks)
	{
	  lda_top = nrows_top;
	  lda_rest = nrows_rest;
	  A_rest_ptr = A_top_ptr + nrows_top * ncols();
	}
      else
	{
	  lda_top = lda();
	  lda_rest = lda();
	  A_rest_ptr = A_top_ptr + nrows_top;
	}
      MatView A_top (nrows_top, ncols(), get(), lda_top);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_top;
    }

    /// Split off and return the bottom block.  Modify *this to be the
    /// "rest" of the matrix.
    MatView split_bottom (const Ordinal nrows_bottom, 
			  const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed && nrows_bottom < 0)
	throw std::invalid_argument ("nrows_bottom < 0");
      if (nrows_bottom > nrows())
	throw std::invalid_argument ("nrows_bottom > nrows");
#endif // TSQR_MATVIEW_DEBUG

      Scalar* const A_rest_ptr = get();
      Scalar* A_bottom_ptr;
      const Ordinal nrows_rest = nrows() - nrows_bottom;
      Ordinal lda_bottom, lda_rest;
      if (b_contiguous_blocks)
	{
	  lda_bottom = nrows_bottom;
	  lda_rest = nrows() - nrows_bottom;
	  A_bottom_ptr = A_rest_ptr + nrows_rest * ncols();
	}
      else
	{
	  lda_bottom = lda();
	  lda_rest = lda();
	  A_bottom_ptr = A_rest_ptr + nrows_rest;
	}
      MatView A_bottom (nrows_bottom, ncols(), A_bottom_ptr, lda_bottom);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_bottom;
    }

    void
    fill (const Scalar value) {
      const Ordinal num_rows = nrows();
      const Ordinal num_cols = ncols();
      const Ordinal stride = lda();

      Scalar* A_j = get();
      for (Ordinal j = 0; j < num_cols; ++j, A_j += stride)
	for (Ordinal i = 0; i < num_rows; ++i)
	  A_j[i] = value;
    }

    /// Deep copy (A := B)
    ///
    /// \note Assumes that B and *this have the same dimensions, that
    ///   their index types are compatible, etc.  See also the
    ///   "generic copy constructor" discussion in Matrix.hpp.
    template< class MatrixViewType >
    void
    copy (MatrixViewType& B) {
      const index_type num_rows = nrows();
      const index_type num_cols = ncols();
      const index_type A_lda = lda();
      const typename MatrixViewType::index_type B_lda = B.lda();

      if (B.nrows() != num_rows || B.ncols() != num_cols)
	throw std::invalid_argument("Dimensions of input matrix B "
				    "are not compatible with *this");
      value_type* A_j = get();
      typename MatrixViewType::pointer_type B_j = B.get();
      for (index_type j = 0; j < num_cols; ++j, A_j += A_lda, B_j += B_lda)
	for (index_type i = 0; i < num_rows; ++i)
	  A_j[i] = B_j[i];
    }

    bool operator== (const MatView& rhs) const {
      return nrows() == rhs.nrows() && ncols() == rhs.ncols() && 
	lda() == rhs.lda() && get() == rhs.get();
    }

    bool operator!= (const MatView& rhs) const {
      return nrows() != rhs.nrows() || ncols() != rhs.ncols() || 
	lda() != rhs.lda() || get() != rhs.get();
    }

  private:
    index_type nrows_, ncols_, lda_;
    value_type* A_;
  };


  /// \class ConstMatView
  /// 
  /// A read-only view of a column-oriented matrix.
  ///
  /// \note Implicit promotion of a MatView to a ConstMatView is
  /// forbidden, because it violates the expectation that ConstMatView
  /// points to a matrix that doesn't change during the computation.
  template< class Ordinal, class Scalar >
  class ConstMatView {
  public:
    typedef Scalar value_type;
    typedef Ordinal index_type;
    typedef const Scalar* pointer_type;

    ConstMatView () : nrows_(0), ncols_(0), A_(NULL), lda_(0) {}

    ConstMatView (const Ordinal num_rows, 
		  const Ordinal num_cols, 
		  const Scalar* const A, 
		  const Ordinal leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      A_(A),
      lda_(leading_dim) 
    {
#ifdef TSQR_MATVIEW_DEBUG
      MatViewVerify< Ordinal, Scalar >::verify (num_rows, num_cols, A, leading_dim);
#endif // TSQR_MATVIEW_DEBUG
    }

    ConstMatView (const ConstMatView& view) :
      nrows_(view.nrows()),
      ncols_(view.ncols()),
      A_(view.get()),
      lda_(view.lda())
    {}

    ConstMatView& operator= (const ConstMatView& view) {
      if (this == &view)
	return *this;
      else
	{
	  nrows_ = view.nrows();
	  ncols_ = view.ncols();
	  A_ = view.get();
	  lda_ = view.lda();
	}
    }

    const Scalar& operator() (const Ordinal i, const Ordinal j) const 
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed)
	{
	  if (i < 0 || i >= nrows())
	    throw std::invalid_argument("Row range invalid");
	  else if (j < 0 || j >= ncols())
	    throw std::invalid_argument("Column range invalid");
	}
      else
	{
	  if (i >= nrows())
	    throw std::invalid_argument("Row range invalid");
	  else if (j >= ncols())
	    throw std::invalid_argument("Column range invalid");
	}
      if (A_ == NULL)
	throw std::logic_error("Attempt to reference NULL data");
#endif // TSQR_MATVIEW_DEBUG
      return A_[i + j*lda()];
    }

    Ordinal nrows() const { return nrows_; }
    Ordinal ncols() const { return ncols_; }
    Ordinal lda() const { return lda_; }
    pointer_type get() const { return A_; }
    bool empty() const { return nrows() == 0 || ncols() == 0; }

    /// Return a "row block" (submatrix of consecutive rows in the
    /// inclusive range [firstRow,lastRow]).
    ConstMatView rowBlock (const Ordinal firstRow, 
			   const Ordinal lastRow) const
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (firstRow < 0 || lastRow >= nrows())
	throw std::invalid_argument ("Row range invalid");
#endif // TSQR_MATVIEW_DEBUG
      return ConstMatView (lastRow - firstRow + 1, ncols(), get() + firstRow, lda());
    }


    /// Split off and return the top block.  Modify *this to be the
    /// "rest" of the matrix.
    ///
    /// \note Only use this method to split off a single cache block.
    ///   It breaks if you try to use it otherwise.
    ///
    /// \param nrows_top [in] Number of rows in the top block (which
    ///   this method returns)
    ///
    /// \param b_contiguous_blocks [in] Whether or not the entries of
    ///   the top block are stored contiguously in *this.  The default
    ///   is no (false).
    ///
    /// \return The top block of nrows_top rows.  Data is a shallow
    ///   copy of the data in *this.
    ConstMatView split_top (const Ordinal nrows_top, 
			    const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed && nrows_top < 0)
	throw std::invalid_argument ("nrows_top < 0");
      if (nrows_top > nrows())
	throw std::invalid_argument ("nrows_top > nrows");
#endif // TSQR_MATVIEW_DEBUG

      pointer_type const A_top_ptr = get();
      pointer_type A_rest_ptr;
      const Ordinal nrows_rest = nrows() - nrows_top;
      Ordinal lda_top, lda_rest;
      if (b_contiguous_blocks)
	{
	  lda_top = nrows_top;
	  lda_rest = nrows_rest;
	  A_rest_ptr = A_top_ptr + nrows_top * ncols();
	}
      else
	{
	  lda_top = lda();
	  lda_rest = lda();
	  A_rest_ptr = A_top_ptr + nrows_top;
	}
      ConstMatView A_top (nrows_top, ncols(), get(), lda_top);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_top;
    }


    /// Split off and return the bottom block.  Modify *this to be the
    /// "rest" of the matrix.
    ConstMatView split_bottom (const Ordinal nrows_bottom, 
			       const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits< Ordinal >::is_signed && nrows_bottom < 0)
	throw std::invalid_argument ("nrows_bottom < 0");
      if (nrows_bottom > nrows())
	throw std::invalid_argument ("nrows_bottom > nrows");
#endif // TSQR_MATVIEW_DEBUG

      pointer_type const A_rest_ptr = get();
      pointer_type A_bottom_ptr;
      const index_type nrows_rest = nrows() - nrows_bottom;
      index_type lda_bottom, lda_rest;
      if (b_contiguous_blocks)
	{
	  lda_bottom = nrows_bottom;
	  lda_rest = nrows() - nrows_bottom;
	  A_bottom_ptr = A_rest_ptr + nrows_rest * ncols();
	}
      else
	{
	  lda_bottom = lda();
	  lda_rest = lda();
	  A_bottom_ptr = A_rest_ptr + nrows_rest;
	}
      ConstMatView A_bottom (nrows_bottom, ncols(), A_bottom_ptr, lda_bottom);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_bottom;
    }

    bool operator== (const ConstMatView& rhs) const {
      return nrows() == rhs.nrows() && ncols() == rhs.ncols() && 
	lda() == rhs.lda() && get() == rhs.get();
    }

    bool operator!= (const ConstMatView& rhs) const {
      return nrows() != rhs.nrows() || ncols() != rhs.ncols() || 
	lda() != rhs.lda() || get() != rhs.get();
    }


  private:
    index_type nrows_, ncols_, lda_;
    pointer_type A_;
  };

} // namespace TSQR


#endif // __TSQR_Tsqr_MatView_hpp
