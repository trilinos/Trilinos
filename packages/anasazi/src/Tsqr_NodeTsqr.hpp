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
