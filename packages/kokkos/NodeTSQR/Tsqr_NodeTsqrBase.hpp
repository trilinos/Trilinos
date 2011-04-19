// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __TSQR_NodeTsqrBase_hpp
#define __TSQR_NodeTsqrBase_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_Matrix.hpp>

#include <Teuchos_Describable.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <cstring> // size_t
#include <limits>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class NodeTsqrBase
  /// \brief Common interface for intranode TSQR.
  ///
  /// NodeTsqrBase provides a generic interface for TSQR operations
  /// within a node ("intranode").
  template<class LocalOrdinal, class Scalar>
  class NodeTsqrBase : public Teuchos::Describable {
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef LocalOrdinal ordinal_type;

    //! Constructor
    NodeTsqrBase() {}

    //! Virtual destructor ensures safe polymorphic destruction.
    virtual ~NodeTsqrBase() {}

    /// \brief One-line description of this object.
    ///
    /// This implements Teuchos::Describable::description().
    /// Subclasses should override this to provide a more specific
    /// description of their implementation.  Subclasses may also
    /// implement Teuchos::Describable::describe(), which for this
    /// class has a simple default implementation that calls
    /// description() with appropriate indenting.
    virtual std::string description () const {
      return std::string("Intranode Tall Skinny QR (TSQR) interface");
    }

    //! Does the factorization make an R factor with nonnegative diagonal?
    virtual bool 
    QR_produces_R_factor_with_nonnegative_diagonal () const = 0;

    /// \brief Cache block A_in into A_out.
    ///
    /// \param nrows [in] Number of rows in A_in and A_out.
    /// \param ncols [in] Number of columns in A_in and A_out.
    /// \param A_out [out] Result of cache-blocking A_in.
    /// \param A_in [in] Matrix to cache block, stored in column-major
    ///   order with leading dimension lda_in.
    /// \param lda_in [in] Leading dimension of A_in.  (See the LAPACK
    ///   documentation for a definition of "leading dimension.")
    ///   lda_in >= nrows.
    virtual void
    cache_block (const LocalOrdinal nrows,
		 const LocalOrdinal ncols, 
		 Scalar A_out[],
		 const Scalar A_in[],
		 const LocalOrdinal lda_in) const = 0;

    /// \brief Un - cache block A_in into A_out.
    ///
    /// A_in is a matrix produced by \c cache_block().  It is
    /// organized as contiguously stored cache blocks.  This method
    /// reorganizes A_in into A_out as an ordinary matrix stored in
    /// column-major order with leading dimension lda_out.
    ///
    /// \param nrows [in] Number of rows in A_in and A_out.
    /// \param ncols [in] Number of columns in A_in and A_out.
    /// \param A_out [out] Result of un-cache-blocking A_in.
    ///   Matrix stored in column-major order with leading
    ///   dimension lda_out.
    /// \param lda_out [in] Leading dimension of A_out.  (See the
    ///   LAPACK documentation for a definition of "leading
    ///   dimension.")  lda_out >= nrows.
    /// \param A_in [in] Matrix to un-cache-block.
    virtual void
    un_cache_block (const LocalOrdinal nrows,
		    const LocalOrdinal ncols,
		    Scalar A_out[],
		    const LocalOrdinal lda_out,		    
		    const Scalar A_in[]) const = 0;

    //! Compute QR factorization (implicitly stored Q factor) of A.
    virtual factor_output_type
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols, 
	    Scalar A[],
	    const LocalOrdinal lda,
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguousCacheBlocks) = 0;

    //! Apply implicit Q factor stored in Q and factorOutput to C.
    virtual void
    apply (const ApplyType& applyType,
	   const LocalOrdinal nrows,
	   const LocalOrdinal ncols_Q,
	   const Scalar Q[],
	   const LocalOrdinal ldq,
	   const FactorOutput& factorOutput,
	   const LocalOrdinal ncols_C,
	   Scalar C[],
	   const LocalOrdinal ldc,
	   const bool contiguousCacheBlocks) = 0;

    //! Compute explicit Q factor from result of factor().
    virtual void
    explicit_Q (const LocalOrdinal nrows,
		const LocalOrdinal ncols_Q,
		const Scalar Q[],
		const LocalOrdinal ldq,
		const factor_output_type& factorOutput,
		const LocalOrdinal ncols_C,
		Scalar C[],
		const LocalOrdinal ldc,
		const bool contiguousCacheBlocks) = 0;

    /// \brief Compute Q*B
    ///
    /// Compute matrix-matrix product Q*B, where Q is nrows by ncols
    /// and B is ncols by ncols.  Respect cache blocks of Q.
    virtual void
    Q_times_B (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       Scalar Q[],
	       const LocalOrdinal ldq,
	       const Scalar B[],
	       const LocalOrdinal ldb,
	       const bool contiguousCacheBlocks) const = 0;

    /// \brief Fill the nrows by ncols matrix A with zeros.
    /// 
    /// Fill the matrix A with zeros, in a way that respects the cache
    /// blocking scheme.
    ///
    /// \param nrows [in] Number of rows in A
    /// \param ncols [in] Number of columns in A 
    /// \param A [out] nrows by ncols column-major-order dense matrix 
    ///   with leading dimension lda
    /// \param lda [in] Leading dimension of A: lda >= nrows
    /// \param contiguousCacheBlocks [in] Whether the cache blocks
    ///   in A are stored contiguously.
    virtual void
    fill_with_zeros (const LocalOrdinal nrows,
		     const LocalOrdinal ncols,
		     Scalar A[],
		     const LocalOrdinal lda, 
		     const bool contiguousCacheBlocks) = 0;

    /// \brief Return topmost cache block of C
    ///
    /// \param C [in] Matrix (view), supporting the usual nrows(),
    ///   ncols(), get(), lda() interface.
    /// \param contiguousCacheBlocks [in] Whether the cache blocks
    ///   in C are stored contiguously.
    ///
    /// Return a view of the topmost cache block (on this node) of the
    /// given matrix C.  This is not necessarily square, though it
    /// must have at least as many rows as columns.  For a square
    /// ncols by ncols block, as needed by Tsqr::apply(), do as 
    /// follows:
    /// \code 
    /// MatrixViewType top = this->top_block (C, contig);
    /// MatView< LocalOrdinal, Scalar > square (ncols, ncols, top.get(), top.lda());
    /// \endcode
    template< class MatrixViewType >
    virtual MatrixViewType
    top_block (const MatrixViewType& C, 
	       const bool contiguousCacheBlocks) const = 0;
  };


} // namespace TSQR

#endif // __TSQR_NodeTsqrBase_hpp
