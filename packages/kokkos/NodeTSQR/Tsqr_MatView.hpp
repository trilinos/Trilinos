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

  template< class MatrixViewType1, class MatrixViewType2 >
  void
  matrixCopy (MatrixViewType1& A, const MatrixViewType2& B)
  {
    const typename MatrixViewType1::ordinal_type A_nrows = A.nrows();
    const typename MatrixViewType1::ordinal_type A_ncols = A.ncols();
    if (A_nrows != B.nrows() || A_ncols != B.ncols())
      {
	using std::endl;
	std::ostringstream os;
	os << "matrixCopy: dimensions of A (output matrix) "
	  "and B (input matrix) are not compatible." 
	   << endl
	   << "A is " << A.nrows() << " x " << A.ncols() 
	   << ", and B is " << B.nrows() << " x " << B.ncols() << ".";
	throw std::invalid_argument(os.str());
      }
    for (typename MatrixViewType1::ordinal_type j = 0; j < A_ncols; ++j)
      {
	typename MatrixViewType1::scalar_type* const A_j = &A(0,j);
	const typename MatrixViewType2::scalar_type* const B_j = &B(0,j);
	for (typename MatrixViewType1::ordinal_type i = 0; i < A_nrows; ++i)
	  A_j[i] = B_j[i];
      }
  }

  template< class FirstMatrixViewType, class SecondMatrixViewType >
  bool
  matrix_equal (FirstMatrixViewType& A, SecondMatrixViewType& B)
  {
    if (A.nrows() != B.nrows() || A.ncols() != B.ncols())
      return false;
    
    typedef typename FirstMatrixViewType::ordinal_type first_ordinal_type;
    typedef typename SecondMatrixViewType::ordinal_type second_ordinal_type;
    typedef typename FirstMatrixViewType::pointer_type first_pointer_type;
    typedef typename SecondMatrixViewType::pointer_type second_pointer_type;

    const first_ordinal_type nrows = A.nrows();
    const first_ordinal_type A_lda = A.lda();
    const first_ordinal_type ncols = A.ncols();
    const second_ordinal_type B_lda = B.lda();

    first_pointer_type A_j = A.get();
    second_pointer_type B_j = B.get();

    for (first_ordinal_type j = 0; j < ncols; ++j, A_j += A_lda, B_j += B_lda)
      for (first_ordinal_type i = 0; i < nrows; ++i)
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


  // Forward declaration
  template< class Ordinal, class Scalar >
  class ConstMatView;

  // Forward declaration
  template< class Ordinal, class Scalar >
  class Matrix;

  /// \class MatView
  /// 
  /// A read-and-write view of a column-oriented matrix.
  template< class Ordinal, class Scalar >
  class MatView {
  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef Scalar* pointer_type;

    /// \note g++ with -Wall wants A_ to be initialized after lda_,
    /// otherwise it emits a compiler warning.
    MatView () : nrows_(0), ncols_(0), lda_(0), A_(NULL) {}

    MatView (const Ordinal num_rows, 
	     const Ordinal num_cols, 
	     Scalar* const A, 
	     const Ordinal leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      lda_(leading_dim),
      A_(A)
    {
#ifdef TSQR_MATVIEW_DEBUG
      MatViewVerify< Ordinal, Scalar >::verify (num_rows, num_cols, A, leading_dim);
#endif // TSQR_MATVIEW_DEBUG
    }

    MatView (const MatView& view) :
      nrows_(view.nrows()),
      ncols_(view.ncols()),
      lda_(view.lda()),
      A_(view.get())
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
    fill (const scalar_type& value) 
    {
      const ordinal_type num_rows = nrows();
      const ordinal_type num_cols = ncols();
      const ordinal_type stride = lda();

      scalar_type* A_j = get();
      for (ordinal_type j = 0; j < num_cols; ++j, A_j += stride)
	for (ordinal_type i = 0; i < num_rows; ++i)
	  A_j[i] = value;
    }

    /// Deep copy (A := B)
    ///
    /// \note Assumes that B and *this have the same dimensions
    void
    copy (const MatView< ordinal_type, scalar_type >& B) {
      matrixCopy (*this, B);
    }
    void
    copy (const ConstMatView< ordinal_type, scalar_type >& B) {
      matrixCopy (*this, B);
    }
    void
    copy (const Matrix< ordinal_type, scalar_type >& B) {
      matrixCopy (*this, B);
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
    ordinal_type nrows_, ncols_, lda_;
    scalar_type* A_;
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
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef const Scalar* pointer_type;

    ConstMatView () : nrows_(0), ncols_(0), lda_(0), A_(NULL) {}

    /// \note g++ with -Wall wants A_ to be initialized after lda_,
    /// otherwise it emits a compiler warning.
    ConstMatView (const Ordinal num_rows, 
		  const Ordinal num_cols, 
		  const Scalar* const A, 
		  const Ordinal leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      lda_(leading_dim),
      A_(A)
    {
#ifdef TSQR_MATVIEW_DEBUG
      MatViewVerify< Ordinal, Scalar >::verify (num_rows, num_cols, A, leading_dim);
#endif // TSQR_MATVIEW_DEBUG
    }

    ConstMatView (const ConstMatView& view) :
      nrows_(view.nrows()),
      ncols_(view.ncols()),
      lda_(view.lda()),
      A_(view.get())
    {}

    ConstMatView& operator= (const ConstMatView& view) {
      if (this != &view)
	{
	  nrows_ = view.nrows();
	  ncols_ = view.ncols();
	  lda_ = view.lda();
	  A_ = view.get();
	}
      return *this;
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
      const ordinal_type nrows_rest = nrows() - nrows_bottom;
      ordinal_type lda_bottom, lda_rest;
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
    ordinal_type nrows_, ncols_, lda_;
    pointer_type A_;
  };

} // namespace TSQR


#endif // __TSQR_Tsqr_MatView_hpp
