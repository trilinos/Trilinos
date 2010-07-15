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

#ifndef __TSQR_Tsqr_NodeTsqr_hpp
#define __TSQR_Tsqr_NodeTsqr_hpp

#include <Tsqr_MatView.hpp>
#include <cstring> // size_t

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class LocalOrdinal, class Scalar >
  class NodeTsqrFactorOutput {
  public:
    NodeTsqrFactorOutput () {}
    virtual ~NodeTsqrFactorOutput() = 0;
  };

  template< class LocalOrdinal, class Scalar >
  class NodeTsqr {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef NodeTsqrFactorOutput factor_output_type;

    NodeTsqr(const size_t cache_block_size = 0) {}
    virtual ~NodeTsqr() {}

    /// Whether or not the R factor from the QR factorization has a
    /// nonnegative diagonal.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const;

    virtual size_t
    cache_block_size() const = 0;

    virtual void
    cache_block (const LocalOrdinal nrows,
		 const LocalOrdinal ncols, 
		 Scalar A_out[],
		 const Scalar A_in[],
		 const LocalOrdinal lda_in) const = 0;

    virtual void
    un_cache_block (const LocalOrdinal nrows,
		    const LocalOrdinal ncols,
		    Scalar A_out[],
		    const LocalOrdinal lda_out,		    
		    const Scalar A_in[]) const = 0;

    virtual factor_output_type
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols, 
	    Scalar A[],
	    const LocalOrdinal lda,
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguous_cache_blocks) = 0;

    virtual void
    apply (const std::string& op,
	   const LocalOrdinal nrows,
	   const LocalOrdinal ncols_Q,
	   const Scalar* const Q,
	   const LocalOrdinal ldq,
	   const FactorOutput& factor_output,
	   const LocalOrdinal ncols_C,
	   Scalar* const C,
	   const LocalOrdinal ldc,
	   const bool contiguous_cache_blocks) = 0;

    virtual void
    explicit_Q (const LocalOrdinal nrows,
		const LocalOrdinal ncols_Q,
		const Scalar Q[],
		const LocalOrdinal ldq,
		const factor_output_type* const factor_output,
		const LocalOrdinal ncols_C,
		Scalar C[],
		const LocalOrdinal ldc,
		const bool contiguous_cache_blocks) = 0;

    // Return a view of the topmost cache block (on this node) of the
    // given matrix C.  TSQR::Tsqr needs this for apply().
    template< class MatrixViewType >
    virtual MatrixViewType
    top_block (const MatrixViewType& C, 
	       const bool contiguous_cache_blocks = false) const = 0;

    virtual void
    fill_with_zeros (const LocalOrdinal nrows,
		     const LocalOrdinal ncols,
		     Scalar C[],
		     const LocalOrdinal ldc, 
		     const bool contiguous_cache_blocks = false) = 0;
  };
} // namespace TSQR


#endif __TSQR_Tsqr_NodeTsqr_hpp
