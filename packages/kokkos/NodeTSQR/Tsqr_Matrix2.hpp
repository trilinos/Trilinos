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

#ifndef __TSQR_Tsqr_Matrix2_hpp
#define __TSQR_Tsqr_Matrix2_hpp

#include <Teuchos_ConstTypeTraits.hpp>
#include <cstring> // NULL

// Define for bounds checking and other safety features, undefine for speed.
// #define TSQR_MATVIEW_DEBUG 1

#ifdef TSQR_MATVIEW_DEBUG
#  include <limits>
#endif // TSQR_MATVIEW_DEBUG

#include <sstream>
#include <stdexcept>

namespace TSQR {

  /// \class Matrix
  /// \brief A column-oriented matrix, or a view of one.
  /// \author Mark Hoemmen
  ///
  template<class Ordinal, class Scalar>
  class Matrix {
  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef Scalar* pointer_type;

    //! Default constructor (empty, nonowning).
    Matrix () : nrows_(0), ncols_(0), lda_(0), A_(NULL), owns_ (false) {}

    //! Owning, allocating constructor.
    Matrix (const Ordinal numRows, const Ordinal numCols) : 
      nrows_(numRows), 
      ncols_(numCols), 
      lda_(numRows), 
      A_ (new scalar_type [numRows * numCols]), 
      owns_(true)
    {}

    //! Destructor.
    ~Matrix () 
    {
      if (owns_ && A_ != NULL)
	delete [] A_;
      A_ = NULL;
    }

    // Nonallocating constructor (nonowning by default).
    Matrix (const Ordinal numRows, 
	    const Ordinal numCols, 
	    Scalar* const A, 
	    const Ordinal stride,
	    const bool owns=false)
      nrows_(numRows),
      ncols_(numCols),
      lda_(stride),
      A_(A),
      owns_(owns)
    {}

    //! A view of the matrix.  Const since the pointer address is.
    Matrix 
    view () const 
    {
      return Matrix (nrows(), ncols(), get(), lda(), false);
    }

    void
    assign (const Matrix& rhs)
    {
      copy_matrix (nrows(), ncols(), get(), lda(), rhs.get(), rhs.ldb());
    }

    /// \brief Copy constructor.  
    ///
    /// Nonowning matrices are copied by view (shallowly).  Owning
    /// matrices are copied deeply.
    Matrix (const Matrix& rhs)
    {
      nrows_ = rhs.nrows();
      ncols_ = rhs.ncols();
      lda_ = rhs.lda();
      owns_ = rhs.owns();
      if (rhs.owns())
	{
	  A_ = new scalar_type [nrows_ * ncols_];
	  assign (rhs);
	}
      else
	A_ = rhs.get();
    }

    Matrix& operator= (const MatView& rhs) 
    {
      if (this != &view)
	{
	  nrows_ = rhs.nrows();
	  ncols_ = rhs.ncols();
	  lda_ = rhs.lda();
	  owns_ = rhs.owns();
	  if (rhs.owns())
	    {
	      A_ = new scalar_type [nrows_ * ncols_];
	      assign (rhs);
	    }
	  else
	    A_ = rhs.get();
	}
      return *this;
    }

    /// \brief Return a reference to element (i,j) (zero-indexed).
    ///
    /// \note The function is const, only because returning a
    /// reference to the matrix data doesn't change any members of
    /// *this.  Of course one may use the resulting reference to
    /// change an entry in the matrix.
    Scalar& operator() (const Ordinal i, const Ordinal j) const 
    {
      return (get())[i + j*lda()];
    }

    Ordinal nrows() const { return nrows_; }
    Ordinal ncols() const { return ncols_; }
    Ordinal lda() const { return lda_; }
    bool owns() const { return owns_; }

    /// \brief Return a pointer to the matrix data.
    ///
    /// \note The function is const because the pointer address is
    ///   const.  Of course one may use the resulting pointer to
    ///   change entries in the matrix.
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

#endif // __TSQR_Tsqr_Matrix2_hpp
