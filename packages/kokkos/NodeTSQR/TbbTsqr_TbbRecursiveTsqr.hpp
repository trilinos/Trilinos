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

#ifndef __TSQR_TbbRecursiveTsqr_hpp
#define __TSQR_TbbRecursiveTsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_CacheBlocker.hpp>
#include <Tsqr_SequentialTsqr.hpp>
#include <TbbTsqr_Partitioner.hpp>

#include <stdexcept>
#include <string>
#include <utility> // std::pair
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    /// \class TbbRecursiveTsqr
    /// \brief Non-parallel "functioning stub" implementation of \c TbbTsqr.
    ///
    template< class LocalOrdinal, class Scalar >
    class TbbRecursiveTsqr {
    public:
      /// \brief Constructor.
      ///
      /// \param num_cores [in] Maximum parallelism to use (i.e.,
      ///   maximum number of partitions into which to divide the
      ///   matrix to factor).
      ///
      /// \param cache_size_hint [in] Approximate cache size in bytes
      ///   per CPU core.  A hint, not a command.  If zero, set to a
      ///   reasonable default.
      TbbRecursiveTsqr (const size_t num_cores = 1,
			const size_t cache_size_hint = 0);

      /// Number of cores to use to solve the problem (i.e., number of
      /// subproblems into which to divide the main problem, to solve
      /// it in parallel).
      size_t ncores() const { return ncores_; }

      /// \brief Cache size hint (in bytes) used for the factorization.
      ///
      /// This method is deprecated, because the name is misleading.
      /// Please call \c cache_size_hint() instead.
      size_t cache_block_size() const { return seq_.cache_size_hint(); }

      //! Cache size hint (in bytes) used for the factorization.
      size_t cache_size_hint() const { return seq_.cache_size_hint(); }

      //! Results of SequentialTsqr for each core.
      typedef typename SequentialTsqr<LocalOrdinal, Scalar>::FactorOutput SeqOutput;

      /// \typedef ParOutput
      /// \brief Array of ncores "local tau arrays" from parallel TSQR.
      ///
      /// Local Q factors are stored in place.
      typedef std::vector<std::vector<Scalar> > ParOutput;

      /// \typedef FactorOutput
      /// \brief Return type of factor().
      ///
      /// factor() returns a pair: the results of SequentialTsqr for
      /// data on each core, and the results of combining the data on
      /// the cores.
      typedef typename std::pair<std::vector<SeqOutput>, ParOutput> FactorOutput;

      /// Copy the nrows by ncols matrix A_in (with leading dimension
      /// lda_in >= nrows) into A_out, such that cache blocks are
      /// arranged contiguously in memory.
      void
      cache_block (const LocalOrdinal nrows,
		   const LocalOrdinal ncols, 
		   Scalar A_out[],
		   const Scalar A_in[],
		   const LocalOrdinal lda_in) const;

      /// Copy the nrows by ncols matrix A_in, whose cache blocks are
      /// arranged contiguously in memory, into A_out (with leading
      /// dimension lda_out >= nrows), which is in standard
      /// column-major order.
      void
      un_cache_block (const LocalOrdinal nrows,
		      const LocalOrdinal ncols,
		      Scalar A_out[],
		      const LocalOrdinal lda_out,		    
		      const Scalar A_in[]) const;

      /// Compute the QR factorization of the nrows by ncols matrix A
      /// (with leading dimension lda >= nrows), returning a
      /// representation of the Q factor (which includes data stored
      /// in-place in A), and overwriting R (an ncols by ncols matrix
      /// in column-major order with leading dimension ldr >= ncols)
      /// with the R factor.
      FactorOutput 
      factor (const LocalOrdinal nrows,
	      const LocalOrdinal ncols, 
	      Scalar A[],
	      const LocalOrdinal lda,
	      Scalar R[],
	      const LocalOrdinal ldr,
	      const bool contiguous_cache_blocks) const;

      /// Apply the Q factor computed by factor() (which see) to the
      /// nrows by ncols_C matrix C, with leading dimension ldc >=
      /// nrows.
      void
      apply (const std::string& op,
	     const LocalOrdinal nrows,
	     const LocalOrdinal ncols_C,
	     Scalar C[],
	     const LocalOrdinal ldc,
	     const LocalOrdinal ncols_Q,
	     const Scalar Q[],
	     const LocalOrdinal ldq,
	     const FactorOutput& factor_output,
	     const bool contiguous_cache_blocks) const;

      /// Compute the explicit representation of the Q factor computed
      /// by factor().
      void 
      explicit_Q (const LocalOrdinal nrows,
		  const LocalOrdinal ncols_Q_in,
		  const Scalar Q_in[],
		  const LocalOrdinal ldq_in,
		  const LocalOrdinal ncols_Q_out,
		  Scalar Q_out[],
		  const LocalOrdinal ldq_out,
		  const FactorOutput& factor_output,
		  const bool contiguous_cache_blocks) const;

    private:
      size_t ncores_;
      TSQR::SequentialTsqr<LocalOrdinal, Scalar> seq_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;

      typedef MatView<LocalOrdinal, Scalar> mat_view;
      typedef ConstMatView<LocalOrdinal, Scalar> const_mat_view;
      typedef std::pair<const_mat_view, const_mat_view> const_split_t;
      typedef std::pair<mat_view, mat_view> split_t;
      typedef std::pair<const_mat_view, mat_view> top_blocks_t;
      typedef std::vector<top_blocks_t> array_top_blocks_t;

      void
      explicit_Q_helper (const size_t P_first, 
			 const size_t P_last,
			 MatView< LocalOrdinal, Scalar >& Q_out,
			 const bool contiguous_cache_blocks) const;

      /// \return MatView of the topmost block (good for combining the
      ///   R factors and extracting the final R factor result).
      MatView<LocalOrdinal, Scalar>
      factor_helper (const size_t P_first, 
		     const size_t P_last,
		     const size_t depth,
		     MatView< LocalOrdinal, Scalar > A,
		     std::vector< SeqOutput >& seq_outputs,
		     ParOutput& par_outputs,
		     Scalar R[],
		     const LocalOrdinal ldr,
		     const bool contiguous_cache_blocks) const;

      bool
      apply_helper_empty (const size_t P_first,
			  const size_t P_last,
			  const_mat_view &Q,
			  mat_view& C) const;

      /// Build array of ncores() blocks, one for each partition.
      /// Each block is the topmost block in that partition.  This is
      /// useful for apply_helper.
      void
      build_partition_array (const size_t P_first,
			     const size_t P_last,
			     array_top_blocks_t& top_blocks,
			     const_mat_view& Q,
			     mat_view& C,
			     const bool contiguous_cache_blocks) const;

      /// Apply Q (not Q^T or Q^H, which is why we don't ask for "op")
      /// to C.
      void
      apply_helper (const size_t P_first, 
		    const size_t P_last,
		    const_mat_view Q,
		    mat_view C,
		    array_top_blocks_t& top_blocks, 
		    const FactorOutput& factor_output,
		    const bool contiguous_cache_blocks) const;

      /// Apply Q^T or Q^H to C.
      ///
      /// \return Views of the topmost partitions of Q resp. C.
      std::pair< ConstMatView< LocalOrdinal, Scalar >, MatView< LocalOrdinal, Scalar > >
      apply_transpose_helper (const std::string& op,
			      const size_t P_first, 
			      const size_t P_last,
			      const_mat_view Q,
			      mat_view C,
			      const FactorOutput& factor_output,
			      const bool contiguous_cache_blocks) const;

      void 
      factor_pair (const size_t P_top,
		   const size_t P_bot,
		   mat_view& A_top,
		   mat_view& A_bot,
		   std::vector< std::vector< Scalar > >& par_outputs,
		   const bool contiguous_cache_blocks) const;

      void
      apply_pair (const std::string& trans,
		  const size_t P_top,
		  const size_t P_bot,
		  const_mat_view& Q_bot,
		  const std::vector< std::vector< Scalar > >& tau_arrays,
		  mat_view& C_top,
		  mat_view& C_bot,
		  const bool contiguous_cache_blocks) const;

      void 
      cache_block_helper (MatView< LocalOrdinal, Scalar >& A_out,
			  ConstMatView< LocalOrdinal, Scalar >& A_in,
			  const size_t P_first,
			  const size_t P_last) const;

      void 
      un_cache_block_helper (MatView< LocalOrdinal, Scalar >& A_out,
			     const ConstMatView< LocalOrdinal, Scalar >& A_in,
			     const size_t P_first,
			     const size_t P_last) const;

    }; // class TbbRecursiveTsqr
  } // namespace TBB
} // namespace TSQR

#include <TSQR/TBB/TbbRecursiveTsqr_Def.hpp>

#endif // __TSQR_TbbRecursiveTsqr_hpp
