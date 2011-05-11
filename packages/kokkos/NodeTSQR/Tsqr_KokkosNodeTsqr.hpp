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

#ifndef __TSQR_KokkosNodeTsqr_hpp
#define __TSQR_KokkosNodeTsqr_hpp

#include <Tsqr_CacheBlocker.hpp>
#include <Tsqr_Combine.hpp>
#include <Tsqr_NodeTsqr.hpp>

#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_ScalarTraits.hpp>

//#define KNR_DEBUG 1
#ifdef KNR_DEBUG
#  include <iostream>
#endif // KNR_DEBUG

namespace TSQR {

  namespace details {

    /// \brief Half-exclusive range of my partition's cache block indices.
    ///
    /// \c FactorFirstPass (used by the factor() method of \c
    /// KokkosNodeTsqr) breaks up the matrix into contiguous
    /// partitions of row blocks.  The index argument of Kokkos'
    /// parallel_for is the (zero-based) partition index.  This
    /// function returns the half-exclusive range of the cache block
    /// indices belonging to the partition partitionIndex.
    ///
    /// \param numRows [in] Number of rows in the matrix.
    /// \param numCols [in] Number of columns in the matrix.
    /// \param partitionIndex [in] Zero-based index of the partition.
    ///   This is specifically an int and not a LocalOrdinal, because
    ///   partition indices are arguments to Kokkos Node API methods
    ///   parallel_for and parallel_reduce.  Cache block indices are
    ///   of LocalOrdinal type and should not be mixed with partition
    ///   indices, even though in most cases LocalOrdinal == int.
    /// \param numPartitions [in] Total number of partitions; a
    ///   positive integer.
    /// \param strategy [in] The cache blocking strategy to use.
    ///   
    /// \return (start cache block index, end cache block index).
    ///   This is a half-exclusive range: it does not include the end
    ///   point.  Thus, if the two indices are equal, the range is
    ///   empty.
    template<class LocalOrdinal, class Scalar>
    std::pair<LocalOrdinal, LocalOrdinal>
    cacheBlockIndexRange (const LocalOrdinal numRows,
			  const LocalOrdinal numCols,
			  const int partitionIndex,
			  const int numPartitions,
			  const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy)
    {
#ifdef KNR_DEBUG
      using std::cerr;
      using std::endl;
      // cerr << "cacheBlockIndexRange(numRows=" << numRows 
      // 	   << ", numCols=" << numCols
      // 	   << ", partitionIndex=" << partitionIndex 
      // 	   << ", numPartitions=" << numPartitions
      // 	   << ", strategy)" << endl;
#endif // KNR_DEBUG

      // The input index is a zero-based index of the current
      // partition (not the "current cache block" -- a partition
      // contains zero or more cache blocks).  If the input index is
      // out of range, then return, since there is nothing to do.
      //
      // The nice thing about partitioning over cache blocks is that
      // the cache blocking strategy guarantees that exactly one of
      // the following is true:
      // 
      // 1. The partition is empty (contains zero cache blocks)
      // 2. All cache blocks in the partition are valid (none 
      //    contains more columns than rows)

      // Return an empty partition (an empty cache block range) if
      // the partition index is out of range.
      if (partitionIndex >= numPartitions)
	return std::make_pair (LocalOrdinal(0), LocalOrdinal(0));

      const LocalOrdinal numRowsCacheBlock = 
	strategy.cache_block_num_rows (numCols);
      const LocalOrdinal numCacheBlocks = 
	strategy.num_cache_blocks (numRows, numCols, numRowsCacheBlock);

#ifdef KNR_DEBUG
      // cerr << "numRowsCacheBlock=" << numRowsCacheBlock
      // 	   << ", numCacheBlocks=" << numCacheBlocks
      // 	   << endl;
#endif // KNR_DEBUG

      // Figure out how many cache blocks my partition contains.  If
      // the number of partitions doesn't evenly divide the number
      // of cache blocks, we spread out the remainder among the
      // first few threads.
      const LocalOrdinal quotient = numCacheBlocks / numPartitions;
      const LocalOrdinal remainder = numCacheBlocks - quotient * numPartitions;
      const LocalOrdinal myNumCacheBlocks = 
	(partitionIndex < remainder) ? (quotient + 1) : quotient;

#ifdef KNR_DEBUG
      // cerr << "Partition " << partitionIndex << ": quotient=" << quotient 
      // 	   << ", remainder=" << remainder << ", myNumCacheBlocks=" 
      // 	   << myNumCacheBlocks << endl;
#endif // KNR_DEBUG

      // If there are no cache blocks, there is nothing to factor.
      // Return an empty cache block range to indicate this.
      if (myNumCacheBlocks == 0)
	return std::make_pair (LocalOrdinal(0), LocalOrdinal(0));

      // Index of my first cache block (inclusive).
      const LocalOrdinal myFirstCacheBlockIndex = 
	(partitionIndex < remainder) ? 
	partitionIndex * (quotient+1) :
	remainder * (quotient+1) + (partitionIndex - remainder) * quotient;
      // Index of my last cache block (exclusive).
      const LocalOrdinal myLastCacheBlockIndex = 
	(partitionIndex+1 < remainder) ? 
	(partitionIndex+1) * (quotient+1) :
	remainder * (quotient+1) + (partitionIndex+1 - remainder) * quotient;
      // Sanity check.
      if (myLastCacheBlockIndex <= myFirstCacheBlockIndex)
	{
	  std::ostringstream os;
	  os << "Partition " << (partitionIndex+1) << " of " 
	     << numPartitions << ":  My range of cache block indices [" 
	     << myFirstCacheBlockIndex << ", " << myLastCacheBlockIndex 
	     << ") is empty.";
	  throw std::logic_error(os.str());
	}
      return std::make_pair (myFirstCacheBlockIndex, myLastCacheBlockIndex);
    }


    /// \class FactorFirstPass
    /// \brief First pass of KokkosNodeTsqr's factorization.
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class FactorFirstPass {
    private:
      MatView<LocalOrdinal, Scalar> A_;
      // While tauArrays_ is shared among tasks (i.e., partitions),
      // there are no race conditions among entries, since each
      // partition writes its own entry.  Ditto for topBlocks_.
      std::vector<std::vector<Scalar> >& tauArrays_;
      std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      std::vector<Scalar>
      factorFirstCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
			     const MatView<LocalOrdinal, Scalar>& A_top,
			     std::vector<Scalar>& work)
      {
	std::vector<Scalar> tau (A_top.ncols());

	// We should only call this if A_top.ncols() > 0 and therefore
	// work.size() > 0, but we've already checked for that, so we
	// don't have to check again.
	combine.factor_first (A_top.nrows(), A_top.ncols(), A_top.get(), 
			      A_top.lda(), &tau[0], &work[0]);
	return tau;
      }

      std::vector<Scalar>
      factorCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
			const MatView<LocalOrdinal, Scalar>& A_top, 
			const MatView<LocalOrdinal, Scalar>& A_cur,
			std::vector<Scalar>& work)
      {
	std::vector<Scalar> tau (A_top.ncols());

	// We should only call this if A_top.ncols() > 0 and therefore
	// tau.size() > 0 and work.size() > 0, but we've already
	// checked for that, so we don't have to check again.
	combine.factor_inner (A_cur.nrows(), A_top.ncols(), 
			      A_top.get(), A_top.lda(),
			      A_cur.get(), A_cur.lda(),
			      &tau[0], &work[0]);
	return tau;
      }

      /// \brief Factor the given cache block range using sequential TSQR.
      ///
      /// \param cbIndices [in] Half-exclusive range of cache block indices.
      /// \param partitionIndex [in] Zero-based index of my partition.
      ///
      /// \return A view of the top block of the cache block range.
      MatView<LocalOrdinal, Scalar>
      factor (const std::pair<LocalOrdinal, LocalOrdinal> cbIndices,
	      const int partitionIndex)
      {
#ifdef KNR_DEBUG
	using std::cerr;
	using std::endl;
#endif // KNR_DEBUG

	typedef MatView<LocalOrdinal, Scalar> view_type;
	typedef CacheBlockRange<view_type> range_type;

	// Workspace is created here, because it must not be shared
	// among threads.
	std::vector<Scalar> work (A_.ncols());

	// Range of cache blocks to factor.
	range_type cbRange (A_, strategy_, 
			    cbIndices.first, 
			    cbIndices.second, 
			    contiguousCacheBlocks_);
	// Iterator in the forward direction over the range of cache
	// blocks to factor.
	typedef typename CacheBlockRange<view_type>::iterator range_iter_type;
	range_iter_type cbIter = cbRange.begin();

	// Remember the top (first) block.
	view_type A_top = *cbIter;
	if (A_top.empty())
	  return A_top;
	TEST_FOR_EXCEPTION(cbIndices.first >= cbIndices.second,
			   std::logic_error,
			   "FactorFirstPass::factor: A_top is not empty, but "
			   "the cache block index range " << cbIndices.first 
			   << "," << cbIndices.second << " is empty.  Please "
			   "report this bug to the Kokkos developers.");

	// Current cache block index.
	LocalOrdinal curTauIdx = cbIndices.first;

	// Factor the first cache block.
	Combine<LocalOrdinal, Scalar> combine;
	tauArrays_[curTauIdx++] = factorFirstCacheBlock (combine, A_top, work);

	// Move past the first cache block.
	++cbIter;

	// Number of cache block(s) we have factored thus far.
	LocalOrdinal count = 1;

	// Factor the remaining cache block(s).
	range_iter_type cbEnd = cbRange.end();
	while (cbIter != cbEnd)
	  {
	    view_type A_cur = *cbIter;
	    // Iteration over cache blocks of a partition should
	    // always result in nonempty cache blocks.
	    TEST_FOR_EXCEPTION(A_cur.empty(), std::logic_error,
			       "FactorFirstPass::factor: The current cache bloc"
			       "k (the " << count << "-th to factor in the rang"
			       "e [" << cbIndices.first << ","
			       << cbIndices.second << ") of cache block indices"
			       ") in partition " << (partitionIndex+1) << " (ou"
			       "t of " << numPartitions_ << " partitions) is em"
			       "pty.  Please report this bug to the Kokkos deve"
			       "lopers.");
	    TEST_FOR_EXCEPTION(static_cast<size_t>(curTauIdx) >= tauArrays_.size(),
			       std::logic_error,
			       "FactorFirstPass::factor: curTauIdx (= " 
			       << curTauIdx << ") >= tauArrays_.size() (= " 
			       << tauArrays_.size() << ").  Please report this "
			       "bug to the Kokkos developers.");
	    tauArrays_[curTauIdx++] = 
	      factorCacheBlock (combine, A_top, A_cur, work);
	    ++count;
	    ++cbIter;
	  }
#ifdef KNR_DEBUG
	cerr << "Factored " << count << " cache blocks" << endl;
#endif // KNR_DEBUG
	return A_top;
      }

    public:
      /// \brief Constructor
      ///
      /// \param A [in/out] On input: View of the matrix to factor.  
      ///   On output: (Part of) the implicitly stored Q factor.
      ///   (The other part is tauArrays.)
      /// \param tauArrays [out] Where to write the "TAU" arrays
      ///   (implicit factorization results) for each cache block.
      ///   (TAU is what LAPACK's QR factorization routines call this
      ///   array; see the LAPACK documentation for an explanation.)
      ///   Indexed by the cache block index; one TAU array per cache
      ///   block.
      /// \param strategy [in] Cache blocking strategy to use.
      /// \param numPartitions [in] Number of partitions (positive
      ///   integer), and therefore the maximum parallelism available
      ///   to the algorithm.  Oversubscribing processors is OK, but
      ///   should not be done to excess.  This is an int, and not a
      ///   LocalOrdinal, because it is the argument to Kokkos' 
      ///   parallel_for.
      /// \param contiguousCacheBlocks [in] Whether the cache blocks
      ///   of A are stored contiguously.
      FactorFirstPass (const MatView<LocalOrdinal,Scalar>& A, 
		       std::vector<std::vector<Scalar> >& tauArrays,
		       std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks,
		       const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
		       const int numPartitions,
		       const bool contiguousCacheBlocks = false) :
	A_ (A), 
	tauArrays_ (tauArrays), 
	topBlocks_ (topBlocks),
	strategy_ (strategy), 
	numPartitions_ (numPartitions),
	contiguousCacheBlocks_ (contiguousCacheBlocks)
      {
	TEST_FOR_EXCEPTION(A_.empty(), std::logic_error,
			   "TSQR::FactorFirstPass constructor: A is empty.  "
			   "Please report this bug to the Kokkos developers.");
	TEST_FOR_EXCEPTION(numPartitions < 1, std::logic_error,
			   "TSQR::FactorFirstPass constructor: numPartitions "
			   "must be positive, but numPartitions = " 
			   << numPartitions << ".  Please report this bug to "
			   "the Kokkos developers.");
      }

      /// \brief First pass of intranode TSQR factorization.
      ///
      /// Invoked by Kokkos' parallel_for template method.  This
      /// routine parallelizes over contiguous partitions of the
      /// matrix.  Each partition in turn contains cache blocks.
      /// Partitions do not break up cache blocks.  (This ensures that
      /// the cache blocking scheme is the same as that used by
      /// SequentialTsqr, as long as the cache blocking strategies are
      /// the same.  However, the implicit Q factor is not compatible
      /// with that of SequentialTsqr.)
      ///
      /// This method also saves a view of the top block of the
      /// partition in the topBlocks_ array.  This is useful for the
      /// next factorization pass.
      ///
      /// \param partitionIndex [in] Zero-based index of the
      ///   partition.  If greater than or equal to the number of
      ///   partitions, this routine does nothing.
      ///
      /// \warning This routine almost certainly won't work in CUDA.
      ///   If it does, it won't be efficient.  If you are interested
      ///   in a GPU TSQR routine, please contact the author (Mark
      ///   Hoemmen <mhoemme@sandia.gov>) of this code to discuss the
      ///   possibilities.  For this reason, we have not added the
      ///   KERNEL_PREFIX method prefix.
      ///
      /// \note Unlike typical Kokkos work-data pairs (WDPs) passed
      ///   into parallel_for, this one is not declared inline.  This
      ///   method is heavyweight enough that an inline declaration is
      ///   unlikely to improve performance.
      void 
      execute (const int partitionIndex) 
      {
#ifdef KNR_DEBUG
	using std::cerr;
	using std::endl;
	// cerr << "FactorFirstPass::execute (" << partitionIndex << ")" << endl;
#endif // KNR_DEBUG

	if (partitionIndex < 0 || partitionIndex >= numPartitions_)
	  return;
	else if (A_.empty())
	  return;
	else 
	  {
	    const std::pair<LocalOrdinal, LocalOrdinal> cbIndices = 
	      cacheBlockIndexRange (A_.nrows(), A_.ncols(), partitionIndex, 
				    numPartitions_, strategy_);
#ifdef KNR_DEBUG
	    cerr << "Partition " << partitionIndex 
		 << ": Factoring cache block indices ["
		 << cbIndices.first << ", " << cbIndices.second << ")"
		 << endl;
#endif // KNR_DEBUG
	    // It's legitimate, though suboptimal, for some partitions
	    // not to get any work to do (in this case, not to get any
	    // cache blocks to factor).
	    if (cbIndices.second <= cbIndices.first)
	      return;
	    else
	      topBlocks_[partitionIndex] = factor (cbIndices, partitionIndex);
	  }
      }
    };


    /// \class ApplyFirstPass
    /// \brief "First" pass of applying KokkosNodeTsqr's implicit Q factor.
    /// \author Mark Hoemmen
    ///
    /// We call this ApplyFirstPass as a reminder that this algorithm
    /// has the same form as FactorFirstPass and uses the results of
    /// the latter, even though ApplyFirstPass is really the last pass
    /// of applying the implicit Q factor.
    template<class LocalOrdinal, class Scalar>
    class ApplyFirstPass {
    private:
      ApplyType applyType_;
      ConstMatView<LocalOrdinal, Scalar> Q_;
      const std::vector<std::vector<Scalar> >& tauArrays_;
      const std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks_;
      MatView<LocalOrdinal, Scalar> C_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool explicitQ_, contiguousCacheBlocks_;

      void
      applyFirstCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
			    const ApplyType& applyType,
			    const ConstMatView<LocalOrdinal, Scalar>& Q_top,
			    const std::vector<Scalar>& tau,
			    const MatView<LocalOrdinal, Scalar>& C_top,
			    std::vector<Scalar>& work)
      {
	TEST_FOR_EXCEPTION(tau.size() < static_cast<size_t> (Q_top.ncols()), 
			   std::logic_error,
			   "ApplyFirstPass::applyFirstCacheBlock: tau.size() "
			   "(= " << tau.size() << ") < number of columns " 
			   << Q_top.ncols() << " in the Q factor.  Please "
			   "report this bug to the Kokkos developers.");

	// If we get this far, it's fair to assume that we have
	// checked whether tau and work have nonzero lengths.
	combine.apply_first (applyType, C_top.nrows(), C_top.ncols(),
			     Q_top.ncols(), Q_top.get(), Q_top.lda(),
			     &tau[0], C_top.get(), C_top.lda(), &work[0]);
      }

      void
      applyCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
		       const ApplyType& applyType,
		       const ConstMatView<LocalOrdinal, Scalar>& Q_cur,
		       const std::vector<Scalar>& tau,
		       const MatView<LocalOrdinal, Scalar>& C_top,
		       const MatView<LocalOrdinal, Scalar>& C_cur,
		       std::vector<Scalar>& work)
      {
	TEST_FOR_EXCEPTION(tau.size() < static_cast<size_t> (Q_cur.ncols()), 
			   std::logic_error,
			   "ApplyFirstPass::applyCacheBlock: tau.size() "
			   "(= " << tau.size() << ") < number of columns " 
			   << Q_cur.ncols() << " in the Q factor.  Please "
			   "report this bug to the Kokkos developers.");

	// If we get this far, it's fair to assume that we have
	// checked whether tau and work have nonzero lengths.
	combine.apply_inner (applyType, C_cur.nrows(), C_cur.ncols(), 
			     Q_cur.ncols(), Q_cur.get(), Q_cur.lda(), 
			     &tau[0],
			     C_top.get(), C_top.lda(),
			     C_cur.get(), C_cur.lda(),
			     &work[0]);
      }

      /// \fn apply
      /// \brief Apply the sequential part of the implicit Q factor to C.
      ///
      /// \param applyType [in] Whether we are applying Q, Q^T, or Q^H.
      /// \param cbIndices [in] Half-exclusive range of cache block
      ///   indices.
      /// \param partitionIndex [in] The argument to \c execute(); the
      ///   index of the partition which instance of ApplyFirstPass
      ///   is currently processing.
      void
      apply (const ApplyType& applyType,
	     const std::pair<LocalOrdinal, LocalOrdinal> cbIndices,
	     const int partitionIndex)
      {
#ifdef KNR_DEBUG
	using std::cerr;
	using std::endl;
#endif // KNR_DEBUG
	typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
	typedef MatView<LocalOrdinal, Scalar> view_type;
	typedef CacheBlockRange<const_view_type> const_range_type;
	typedef CacheBlockRange<view_type> range_type;
	typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;

	if (cbIndices.first >= cbIndices.second)
	  return; // My range of cache blocks is empty; nothing to do

	// Q_range: Range of cache blocks in the Q factor.
	// C_range: Range of cache blocks in the matrix C.
	const_range_type Q_range (Q_, strategy_, 
				  cbIndices.first, cbIndices.second,
				  contiguousCacheBlocks_);
	range_type C_range (C_, strategy_, 
			    cbIndices.first, cbIndices.second,
			    contiguousCacheBlocks_);
	TEST_FOR_EXCEPTION(Q_range.empty(), std::logic_error,
			   "Q_range is empty, but the range of cache block "
			   "indices [" << cbIndices.first << ", " 
			   << cbIndices.second << ") is not empty.  Please "
			   "report this bug to the Kokkos developers.");
	TEST_FOR_EXCEPTION(C_range.empty(), std::logic_error,
			   "C_range is empty, but the range of cache block "
			   "indices [" << cbIndices.first << ", " 
			   << cbIndices.second << ") is not empty.  Please "
			   "report this bug to the Kokkos developers.");

	// Task-local workspace array of length C_.ncols().  Workspace
	// must be per task, else there will be race conditions as
	// different tasks attempt to write to and read from the same
	// workspace simultaneously.
	std::vector<Scalar> work (C_.ncols());

	Combine<LocalOrdinal, Scalar> combine;
	if (applyType.transposed())
	  {
	    typename const_range_type::iterator Q_rangeIter = Q_range.begin();
	    typename range_type::iterator C_rangeIter = C_range.begin();
	    TEST_FOR_EXCEPTION(Q_rangeIter == Q_range.end(), std::logic_error,
			       "The Q cache block range claims to be nonempty, "
			       "but the iterator range is empty.  Please report"
			       " this bug to the Kokkos developers.");
	    TEST_FOR_EXCEPTION(C_rangeIter == C_range.end(), std::logic_error,
			       "The C cache block range claims to be nonempty, "
			       "but the iterator range is empty.  Please report"
			       " this bug to the Kokkos developers.");

	    // Q_top: Topmost cache block in the cache block range of Q.
	    // C_top: Topmost cache block in the cache block range of C.
	    const_view_type Q_top = *Q_rangeIter;
	    view_type C_top = *C_rangeIter;
	    if (explicitQ_)
	      { 
		C_top.fill (Teuchos::ScalarTraits<Scalar>::zero());
		if (partitionIndex == 0)
		  for (LocalOrdinal j = 0; j < C_top.ncols(); ++j)
		    C_top(j,j) = Teuchos::ScalarTraits<Scalar>::one();
	      }
	    LocalOrdinal curTauIndex = cbIndices.first;

	    // Apply the first block.
	    applyFirstCacheBlock (combine, applyType, Q_top, 
				  tauArrays_[curTauIndex++], C_top, work);

	    // Apply the rest of the blocks, if any.
	    ++Q_rangeIter;
	    ++C_rangeIter;
	    while (Q_rangeIter != Q_range.end())
	      {
		TEST_FOR_EXCEPTION(C_rangeIter == C_range.end(),
				   std::logic_error,
				   "When applying Q^T or Q^H to C: The Q cache "
				   "block iterator is not yet at the end, but "
				   "the C cache block iterator is.  Please "
				   "report this bug to the Kokkos developers.");
		const_view_type Q_cur = *Q_rangeIter;
		view_type C_cur = *C_rangeIter;
		++Q_rangeIter;
		++C_rangeIter;
		if (explicitQ_)
		  C_cur.fill (Teuchos::ScalarTraits<Scalar>::zero());
		applyCacheBlock (combine, applyType, Q_cur, 
				 tauArrays_[curTauIndex++], 
				 C_top, C_cur, work);
	      }
	  }
	else
	  {
	    // Q_top: Topmost cache block in the cache block range of Q.
	    // C_top: Topmost cache block in the cache block range of C.
	    const_view_type Q_top = *(Q_range.begin());
	    view_type C_top = *(C_range.begin());

	    if (explicitQ_)
	      {
		// We've already filled the top ncols x ncols block of
		// C_top with data (that's the result of applying the
		// internode part of the Q factor via DistTsqr).
		// However, we still need to fill the rest of C_top
		// (everything but the top ncols rows of C_top) with
		// zeros.
		view_type C_top_rest (C_top.nrows() - C_top.ncols(), C_top.ncols(),
				      C_top.get() + C_top.ncols(), C_top.lda());
		C_top_rest.fill (Teuchos::ScalarTraits<Scalar>::zero());
	      }
	    LocalOrdinal curTauIndex = cbIndices.second-1;

	    // When applying Q (rather than Q^T or Q^H), we apply the
	    // cache blocks in reverse order.
	    typename const_range_type::iterator Q_rangeIter = Q_range.rbegin();
	    typename range_type::iterator C_rangeIter = C_range.rbegin();
	    TEST_FOR_EXCEPTION(Q_rangeIter == Q_range.rend(), std::logic_error,
			       "The Q cache block range claims to be nonempty, "
			       "but the iterator range is empty.  Please report"
			       " this bug to the Kokkos developers.");
	    TEST_FOR_EXCEPTION(C_rangeIter == C_range.rend(), std::logic_error,
			       "The C cache block range claims to be nonempty, "
			       "but the iterator range is empty.  Please report"
			       " this bug to the Kokkos developers.");

	    // Equality of cache block range iterators only tests the
	    // cache block index, not reverse-ness.  This means we can
	    // compare a reverse-direction iterator (Q_rangeIter) with
	    // a forward-direction iterator (Q_range.begin()).
	    //
	    // We do this because we need to handle the topmost block
	    // of Q_range separately (applyFirstCacheBlock(), rather
	    // than applyCacheBlock()).
	    while (Q_rangeIter != Q_range.begin())
	      {
		const_view_type Q_cur = *Q_rangeIter;
		view_type C_cur = *C_rangeIter;

		if (explicitQ_)
		  C_cur.fill (Teuchos::ScalarTraits<Scalar>::zero());
#ifdef KNR_DEBUG
		cerr << "tauArrays_[curTauIndex=" << curTauIndex << "].size() = " 
		     << tauArrays_[curTauIndex].size() << endl;
#endif // KNR_DEBUG
		TEST_FOR_EXCEPTION(curTauIndex < cbIndices.first, std::logic_error,
				   "curTauIndex=" << curTauIndex << " out of valid "
				   "range [" << cbIndices.first << "," 
				   << cbIndices.second << ").  Please report this "
				   "bug to the Kokkos developers.");
		applyCacheBlock (combine, applyType, Q_cur, 
				 tauArrays_[curTauIndex--], 
				 C_top, C_cur, work);
		++Q_rangeIter;
		++C_rangeIter;
	      }
	    TEST_FOR_EXCEPTION(curTauIndex < cbIndices.first, std::logic_error,
			       "curTauIndex=" << curTauIndex << " out of valid "
			       "range [" << cbIndices.first << "," 
			       << cbIndices.second << ").  Please report this "
			       "bug to the Kokkos developers.");
#ifdef KNR_DEBUG
	    cerr << "tauArrays_[curTauIndex=" << curTauIndex << "].size() = " 
		 << tauArrays_[curTauIndex].size() << endl;
#endif // KNR_DEBUG
	    // Apply the first block.
	    applyFirstCacheBlock (combine, applyType, Q_top, 
				  tauArrays_[curTauIndex--], C_top, work);
	  }
      }

    public:
      /// \brief Constructor
      ///
      /// \param applyType [in] Whether we are applying Q, Q^T, or Q^H.
      /// \param A [in/out] On input: View of the matrix to factor.  
      ///   On output: (Part of) the implicitly stored Q factor.
      ///   (The other part is tauArrays.)
      /// \param tauArrays [in] Where to write the "TAU" arrays
      ///   (implicit factorization results) for each cache block.
      ///   (TAU is what LAPACK's QR factorization routines call this
      ///   array; see the LAPACK documentation for an explanation.)
      ///   Indexed by the cache block index; one TAU array per cache
      ///   block.
      /// \param strategy [in] Cache blocking strategy to use.
      /// \param numPartitions [in] Number of partitions (positive
      ///   integer), and therefore the maximum parallelism available
      ///   to the algorithm.  Oversubscribing processors is OK, but
      ///   should not be done to excess.  This is an int, and not a
      ///   LocalOrdinal, because it is the argument to Kokkos' 
      ///   parallel_for.
      /// \param contiguousCacheBlocks [in] Whether the cache blocks
      ///   of A are stored contiguously.
      ApplyFirstPass (const ApplyType& applyType,
		      const ConstMatView<LocalOrdinal,Scalar>& Q, 
		      const std::vector<std::vector<Scalar> >& tauArrays,
		      const std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks,
		      const MatView<LocalOrdinal,Scalar>& C,
		      const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
		      const int numPartitions,
		      const bool explicitQ = false,
		      const bool contiguousCacheBlocks = false) :
	applyType_ (applyType),
	Q_ (Q),
	tauArrays_ (tauArrays), 
	topBlocks_ (topBlocks),
	C_ (C),
	strategy_ (strategy), 
	numPartitions_ (numPartitions),
	explicitQ_ (explicitQ),
	contiguousCacheBlocks_ (contiguousCacheBlocks)
      {}

      /// \brief First pass of applying intranode TSQR's implicit Q factor.
      ///
      /// Invoked by Kokkos' parallel_for template method.  This
      /// routine parallelizes over contiguous partitions of the C
      /// matrix.  Each partition in turn contains cache blocks.  We
      /// take care not to break up the cache blocks among partitions;
      /// this ensures that the cache blocking scheme is the same as
      /// SequentialTsqr uses.  (However, the implicit Q factor is not
      /// compatible with that of SequentialTsqr.)
      ///
      /// \param partitionIndex [in] Zero-based index of the partition
      ///   which this instance of ApplyFirstPass is currently
      ///   processing.  If greater than or equal to the number of
      ///   partitions, this routine does nothing.
      ///
      /// \warning This routine almost certainly won't work in CUDA.
      ///   If it does, it won't be efficient.  If you are interested
      ///   in a GPU TSQR routine, please contact the author (Mark
      ///   Hoemmen <mhoemme@sandia.gov>) of this code to discuss the
      ///   possibilities.
      ///
      /// \note Unlike typical Kokkos work-data pairs (WDPs) passed
      ///   into parallel_for, this one is not declared inline.  This
      ///   method is heavyweight enough that an inline declaration is
      ///   unlikely to improve performance.
      void 
      execute (const int partitionIndex) 
      {
	if (partitionIndex < 0 || partitionIndex >= numPartitions_)
	  return;
	else if (Q_.empty())
	  return;
	else if (C_.empty())
	  return;

	// We use the same cache block indices for Q and for C.
	std::pair<LocalOrdinal, LocalOrdinal> cbIndices = 
	  cacheBlockIndexRange (Q_.nrows(), Q_.ncols(), partitionIndex, 
				numPartitions_, strategy_);
	if (cbIndices.second <= cbIndices.first)
	  return;
	{
	  std::pair<size_t, size_t> cbInds (static_cast<size_t> (cbIndices.first),
					    static_cast<size_t> (cbIndices.second));
	  TEST_FOR_EXCEPTION(cbIndices.first < static_cast<LocalOrdinal>(0), 
			     std::logic_error,
			     "TSQR::ApplyFirstPass::execute: cacheBlockIndexRa"
			     "nge(" << Q_.nrows() << ", " << Q_.ncols() << ", "
			     << partitionIndex << ", " << numPartitions_ << ","
			     " strategy) returned a cache block range " 
			     << cbIndices.first << "," << cbIndices.second 
			     << " with negative starting index.  Please report"
			     " this bug to the Kokkos developers.");
	  TEST_FOR_EXCEPTION(cbInds.second > tauArrays_.size(),
			     std::logic_error,
			     "TSQR::ApplyFirstPass::execute: cacheBlockIndexRa"
			     "nge(" << Q_.nrows() << ", " << Q_.ncols() << ", " 
			     << partitionIndex << ", " << numPartitions_ << ","
			     " strategy) returned a cache block range "
			     << cbIndices.first << "," << cbIndices.second 
			     << " with starting index larger than the number o"
			     "f tau arrays " << tauArrays_.size() << ".  Pleas"
			     "e report this bug to the Kokkos developers.");
	}

	apply (applyType_, cbIndices, partitionIndex);
      }

    };


    /// \class CacheBlockWDP
    /// \brief Kokkos work-data pair (WDP) for KokkosNodeTsqr's (un_)cache_block() methods.
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class CacheBlockWDP {
    private:
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef CacheBlockRange<const_view_type> const_range_type;
      typedef CacheBlockRange<view_type> range_type;
      
      const_view_type A_in_;
      view_type A_out_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool unblock_;

      /// \brief Copy one range of cache blocks into another.
      ///
      /// \param cbInputRange [in] Range of input cache blocks.
      /// \param cbOutputRange [out] Range of output cache blocks.
      void
      copyRange (const_range_type& cbInputRange, range_type& cbOutputRange)
      {
	typedef typename const_range_type::iterator input_iter_type;
	typedef typename range_type::iterator output_iter_type;
	
	input_iter_type inputIter = cbInputRange.begin();
	output_iter_type outputIter = cbOutputRange.begin();

	input_iter_type inputEnd = cbInputRange.end();
	output_iter_type outputEnd = cbOutputRange.end();
	
	while (inputIter != inputEnd)
	  {
	    const_view_type A_in_cur = *inputIter;
	    view_type A_out_cur = *outputIter;
	    A_out_cur.copy (A_in_cur);
	    ++inputIter;
	    ++outputIter;
	  }
      }

    public:
      /// \brief Constructor
      ///
      /// \param A_in [in] The matrix to (un-)cache-block.
      /// \param A_out [in/out] Result of (un-)cache-blocking the
      ///   matrix A_in.
      /// \param strategy [in] Cache blocking strategy.
      /// \param numPartitions [in] Number of partitions; maximum
      ///   available parallelism.
      /// \param unblock [in] If false, cache-block A_in (a matrix in
      ///   column-major order) into A_out.  If true, un-cache-block
      ///   A_in into A_out (a matrix in column-major order).
      CacheBlockWDP (const ConstMatView<LocalOrdinal, Scalar> A_in,
		     const MatView<LocalOrdinal, Scalar> A_out,
		     const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
		     const int numPartitions,
		     const bool unblock) :
	A_in_ (A_in), 
	A_out_ (A_out), 
	strategy_ (strategy), 
	numPartitions_ (numPartitions),
	unblock_ (unblock)
      {
	TEST_FOR_EXCEPTION(A_in_.nrows() != A_out_.nrows() || 
			   A_in_.ncols() != A_out_.ncols(), 
			   std::invalid_argument,
			   "A_in and A_out do not have the same dimensions: "
			   "A_in is " << A_in_.nrows() << " by " 
			   << A_in_.ncols() << ", but A_out is " 
			   << A_out_.nrows() << " by " 
			   << A_out_.ncols() << ".");
	TEST_FOR_EXCEPTION(numPartitions_ < 1, 
			   std::invalid_argument,
			   "The number of partitions " << numPartitions_ 
			   << " is not a positive integer.");
      }
      
      /// \brief Method called by Kokkos' parallel_for.
      ///
      /// \param partitionIndex [in] Zero-based index of the partition
      ///   of the matrix.  We parallelize over partitions.
      ///   Partitions respect cache blocks.
      void 
      execute (const int partitionIndex)
      {
	if (partitionIndex < 0 || partitionIndex >= numPartitions_ || 
	    A_in_.empty())
	  return;
	else
	  {
	    typedef std::pair<LocalOrdinal, LocalOrdinal> index_range_type;
	    const index_range_type cbIndices = 
	      cacheBlockIndexRange (A_in_.nrows(), A_in_.ncols(), 
				    partitionIndex, numPartitions_, strategy_);
	    // It's perfectly legal for a partitioning to assign zero
	    // cache block indices to a particular partition.  In that
	    // case, this task has nothing to do.
	    if (cbIndices.first >= cbIndices.second)
	      return;
	    else
	      {
		// If unblock_ is false, then A_in_ is in column-major
		// order, and we want to cache-block it into A_out_.  If
		// unblock_ is true, then A_in_ is cache-blocked, and we
		// want to un-cache-block it into A_out_ (a matrix in
		// column-major order).
		const_range_type inputRange (A_in_, strategy_, cbIndices.first,
					     cbIndices.second, unblock_);
		range_type outputRange (A_out_, strategy_, cbIndices.first, 
					cbIndices.second, ! unblock_);
		copyRange (inputRange, outputRange);
	      }
	  }
      }
    };

    /// \class MultWDP
    /// \brief Kokkos work-data pair (WDP) for \c KokkosNodeTsqr::Q_times_B().
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class MultWDP {
    private:
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef CacheBlockRange<view_type> range_type;

      view_type Q_;
      const_view_type B_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      void
      multBlock (BLAS<LocalOrdinal, Scalar>& blas,
		 const view_type& Q_cur,
		 Matrix<LocalOrdinal, Scalar>& Q_temp)
      {
	const LocalOrdinal numCols = Q_cur.ncols();

	// GEMM doesn't like aliased arguments, so we use a copy.  We
	// only copy the current cache block, rather than all of Q;
	// this saves memory.
	Q_temp.reshape (Q_cur.nrows(), numCols);
	Q_temp.copy (Q_cur);

	// Q_cur := Q_temp * B.
	blas.GEMM ("N", "N", Q_cur.nrows(), numCols, numCols, 
		   Teuchos::ScalarTraits<Scalar>::one(),
		   Q_temp.get(), Q_temp.lda(), B_.get(), B_.lda(),
		   Scalar(0), Q_cur.get(), Q_cur.lda());
      }

      /// \brief Multiply (in place) each cache block in the range by B_.
      ///
      /// \param cbRange [in/out] Range of cache blocks.
      void
      multRange (range_type& cbRange)
      {
	typedef typename range_type::iterator iter_type;
	iter_type iter = cbRange.begin();
	iter_type end = cbRange.end();

	// Temporary storage for the BLAS' matrix-matrix multiply
	// routine (which forbids aliasing of any input argument and
	// the output argument).
	Matrix<LocalOrdinal, Scalar> Q_temp;
	BLAS<LocalOrdinal, Scalar> blas;
	while (iter != end)
	  {
	    view_type Q_cur = *iter;
	    multBlock (blas, Q_cur, Q_temp);
	    ++iter;
	  }
      }

    public:
      /// \brief Constructor
      /// 
      /// \param Q [in/out] Matrix to multiply in place by B.
      /// \param B [in] \f$Q := Q * B\f$.
      /// \param strategy [in] Cache-blocking strategy.
      /// \param numPartitions [in] Number of partitions of the matrix
      ///   Q; maximum available parallelism.
      /// \param contiguousCacheBlocks [in] Whether the cache blocks
      ///   of Q are stored contiguously.
      MultWDP (const MatView<LocalOrdinal, Scalar> Q,
	       const ConstMatView<LocalOrdinal, Scalar> B,
	       const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
	       const int numPartitions,
	       const bool contiguousCacheBlocks) :
	Q_ (Q), 
	B_ (B),
	strategy_ (strategy), 
	numPartitions_ (numPartitions),
	contiguousCacheBlocks_ (contiguousCacheBlocks)
      {}

      /// \brief Method called by Kokkos' parallel_for.
      ///
      /// \param partitionIndex [in] Zero-based index of the partition
      ///   of the matrix.  We parallelize over partitions.
      ///   Partitions respect cache blocks.
      void
      execute (const int partitionIndex) 
      {
	if (partitionIndex < 0 || partitionIndex >= numPartitions_ || 
	    Q_.empty())
	  return;
	else
	  {
	    typedef std::pair<LocalOrdinal, LocalOrdinal> index_range_type;
	    const index_range_type cbIndices = 
	      cacheBlockIndexRange (Q_.nrows(), Q_.ncols(), partitionIndex, 
				    numPartitions_, strategy_);
	    if (cbIndices.first >= cbIndices.second)
	      return;
	    else
	      {
		range_type range (Q_, strategy_, cbIndices.first, 
				  cbIndices.second, contiguousCacheBlocks_);
		multRange (range);
	      }
	  }
      }
    };


    /// \class FillWDP
    /// \brief Kokkos work-data pair (WDP) for \c KokkosNodeTsqr::fill_with_zeros().
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class FillWDP {
    private:
      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef CacheBlockRange<view_type> range_type;

      view_type A_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      const Scalar value_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      //! Fill (in place) each cache block in the range with value.
      void
      fillRange (range_type& cbRange, const Scalar value)
      {
	typedef typename range_type::iterator iter_type;
	iter_type iter = cbRange.begin();
	iter_type end = cbRange.end();
	while (iter != end)
	  {
	    view_type A_cur = *iter;
	    A_cur.fill (value);
	    ++iter;
	  }
      }

    public:
      /// \brief Constructor
      /// 
      /// \param A [in/out] Matrix to fill with the value.
      /// \param strategy [in] Cache-blocking strategy.
      /// \param value [in] The value with which to fill A.
      /// \param numPartitions [in] Number of partitions of 
      ///   the matrix A; maximum available parallelism.
      /// \param contiguousCacheBlocks [in] Whether the cache 
      ///   blocks of A are stored contiguously.
      FillWDP (const MatView<LocalOrdinal, Scalar> A,
	       const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
	       const Scalar value,
	       const int numPartitions,
	       const bool contiguousCacheBlocks) :
	A_ (A),
	strategy_ (strategy), 
	value_ (value),
	numPartitions_ (numPartitions),
	contiguousCacheBlocks_ (contiguousCacheBlocks)
      {}

      /// \brief Method called by Kokkos' parallel_for.
      ///
      /// \param partitionIndex [in] Zero-based index of the partition
      ///   of the matrix.  We parallelize over partitions.
      ///   Partitions respect cache blocks.
      void
      execute (const int partitionIndex) 
      {
	if (partitionIndex < 0 || partitionIndex >= numPartitions_ || 
	    A_.empty())
	  return;
	else
	  {
	    typedef std::pair<LocalOrdinal, LocalOrdinal> index_range_type;
	    const index_range_type cbIndices = 
	      cacheBlockIndexRange (A_.nrows(), A_.ncols(), partitionIndex, 
				    numPartitions_, strategy_);
	    if (cbIndices.first >= cbIndices.second)
	      return;
	    else 
	      {
		range_type range (A_, strategy_, cbIndices.first, 
				  cbIndices.second, contiguousCacheBlocks_);
		fillRange (range, value_);
	      }
	  }
      }
    };
  } // namespace details

  /// \class KokkosNodeTsqrFactorOutput
  /// \brief Part of KokkosNodeTsqr's implicit Q representation.
  /// \author Mark Hoemmen
  ///
  /// The \c KokkoNodeTsqr::factor() method represents the Q factor of
  /// the matrix A implicitly.  Part of that representation is in the
  /// A matrix on output, and the other part is returned as an object
  /// of this type.  The apply() and explicit_Q() methods need both
  /// parts of the implicit Q representation in order to do their
  /// work.
  template<class LocalOrdinal, class Scalar>
  struct KokkosNodeTsqrFactorOutput {
    /// \brief Constructor
    ///
    /// \param theNumCacheBlocks [in] Total number of cache blocks
    ///   (over all partitions).
    /// \param theNumPartitions [in] Number of partitions.  This is
    ///   an int because partition indices are ints, and the latter
    ///   are ints because they end up as range arguments to Kokkos'
    ///   parallel_for.
    KokkosNodeTsqrFactorOutput (const size_t theNumCacheBlocks, 
				const int theNumPartitions) :
      firstPassTauArrays (theNumCacheBlocks)
    {
      // Protect the cast to size_t from a negative number of
      // partitions.
      TEST_FOR_EXCEPTION(theNumPartitions < 1, std::invalid_argument,
			 "TSQR::KokkosNodeTsqrFactorOutput: Invalid number of "
			 "partitions " << theNumPartitions << "; number of "
			 "partitions must be a positive integer.");
      // If there's only one partition, we don't even need a second
      // pass (it's just sequential TSQR), and we don't need a TAU
      // array for the top partition.
      secondPassTauArrays.resize (static_cast<size_t> (theNumPartitions-1));
      topBlocks.resize (static_cast<size_t> (theNumPartitions));
    }

    //! Total number of cache blocks in the matrix (over all partitions).
    int numCacheBlocks() const { return firstPassTauArrays.size(); }

    //! Number of partitions of the matrix; max available parallelism.
    int numPartitions() const { return topBlocks.size(); }

    //! TAU arrays from the first pass; one per cache block.
    std::vector<std::vector<Scalar> > firstPassTauArrays;

    /// \brief TAU arrays from the second pass.
    ///
    /// There is one TAU array per partition, except for the topmost
    /// partition.
    ///
    /// For now, KokkosNodeTsqr::factor() uses only two passes over
    /// the matrix.  firstPassTauArrays contains the result of the
    /// pass over cache blocks, and secondPassTauArrays contains the
    /// result of combining the upper triangular R factors from the
    /// first pass.  Later, we may add more passes, in which case we
    /// will likely combine firstPassTauArrays and secondPassTauArrays
    /// into a single std::vector (variable number of passes) or
    /// Teuchos::Tuple (fixed number of passes).
    std::vector<std::vector<Scalar> > secondPassTauArrays;

    /// \brief Views of the topmost cache blocks in each partition.
    ///
    /// One entry for each partition.
    std::vector<MatView<LocalOrdinal, Scalar> > topBlocks;
  };

  /// \class KokkosNodeTsqr
  /// \brief Intranode TSQR using the Kokkos Node API.
  /// \author Mark Hoemmen
  ///
  /// This version of Intranode TSQR factors the matrix in two passes.
  /// The first pass parallelizes over partitions, doing Sequential
  /// TSQR over each partition.  The second pass combines the R
  /// factors from the partitions, and is not currently parallel.
  /// Thus, the overall algorithm is similar to that of \c TbbTsqr,
  /// except that:
  /// - TbbTsqr partitions differently; KokkosNodeTsqr's partitions
  ///   use the same layout of cache blocks as SequentialTsqr, whereas
  ///   TbbTsqr uses a different layout.
  /// - TbbTsqr reduces the R factors in parallel; it only needs one
  ///   "pass."
  ///
  template<class LocalOrdinal, class Scalar, class NodeType>    
  class KokkosNodeTsqr : 
    public NodeTsqr<LocalOrdinal, Scalar, KokkosNodeTsqrFactorOutput<LocalOrdinal, Scalar> >,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef Scalar scalar_type;
    typedef NodeType node_type;

    /// \typedef FactorOutput
    /// \brief Part of the implicit Q representation returned by factor().
    typedef typename NodeTsqr<LocalOrdinal, Scalar, KokkosNodeTsqrFactorOutput<LocalOrdinal, Scalar> >::factor_output_type FactorOutput;

  private:
    //! Implementation of fundamental TSQR kernels.
    Combine<LocalOrdinal, Scalar> combine_;

    //! Workspace for Combine operations.
    mutable std::vector<Scalar> work_;

    //! Pointer to the Kokkos Node object.
    Teuchos::RCP<const node_type> node_;

    //! Cache blocking strategy.
    CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;

    /// \brief Number of partitions; max available parallelism.
    ///
    /// The number of partitions is an int rather than a LocalOrdinal,
    /// to ensure that it is always stored in the ParameterList with
    /// the same type, despite the type of LocalOrdinal.  Besides,
    /// Kokkos wants an int anyway.
    int numPartitions_;

    /// \brief Default number of partitions.
    ///
    /// \param node [in] Kokkos Node object to be used for intranode TSQR.
    ///
    /// This method may in the future try to "learn" the optimal
    /// number of partitions.  For now, it's a constant.  Later, we
    /// may even try to "learn" the best value, perhaps even at
    /// runtime.  As a result, this method may not necessarily return
    /// the same value each time it is called.
    ///
    /// \warning We may change this method to take an RCP to a const
    ///   Kokkos node_type instance, if the Kokkos Node API later
    ///   supports queries for available computational resources
    ///   (e.g., number of CPU cores per node).
    int 
    defaultNumPartitions () const
    {
      // Currently the Kokkos Node API does not give us access to the
      // amount of available parallelism, so we return a constant.
      return 8;
    }

  public:
    /// \brief Constructor
    /// 
    /// \param node [in] Pointer to a Kokkos Node instance.
    /// \param params [in/out] List of parameters.  Missing parameters
    ///   will be filled in with default values.
    KokkosNodeTsqr (const Teuchos::RCP<const node_type>& node,
		    const Teuchos::RCP<Teuchos::ParameterList>& params) :
      node_ (node)
    {
      setParameterList (params);
    }

    /// \brief One-line description of this object.
    ///
    /// This implements Teuchos::Describable::description().
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "KokkosNodeTsqr<LocalOrdinal=" 
	 << TypeNameTraits<LocalOrdinal>::name()
	 << ", Scalar=" 
	 << TypeNameTraits<Scalar>::name()
	 << ", NodeType="
	 << TypeNameTraits<node_type>::name()
	 << ">: \"Cache Size Hint\"=" << strategy_.cache_size_hint()
	 << ", \"Size of Scalar\"=" << strategy_.size_of_scalar()
	 << ", \"Num Partitions\"=" << numPartitions_;
      return os.str();
    }

    /// \brief Validate and read in parameters.
    ///
    /// \param paramList [in/out] On input: non-null parameter list
    ///   containing zero or more of the parameters in \c
    ///   getValidParameters().  On output: missing parameters (i.e.,
    ///   parameters in \c getValidParameters() but not in the input
    ///   list) are filled in with default values.
    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& paramList) 
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;

      RCP<ParameterList> plist;
      if (paramList.is_null())
	plist = rcp (new ParameterList (*getValidParameters ()));
      else
	{
	  plist = paramList;
	  plist->validateParametersAndSetDefaults (*getValidParameters ());
	}
      // Get values of parameters.  We do this "transactionally" so
      // that (except for validation and filling in defaults above)
      // this method has the strong exception guarantee (it either
      // returns, or throws an exception with no externally visible
      // side effects).
      size_t cacheSizeHint, sizeOfScalar;
      int numPartitions;
      try {
	cacheSizeHint = plist->get<size_t> ("Cache Size Hint");
	sizeOfScalar = plist->get<size_t> ("Size of Scalar");
	numPartitions = plist->get<int> ("Num Partitions");
      } catch (Teuchos::Exceptions::InvalidParameter& e) {
	std::ostringstream os;
	os << "Failed to read default parameters after setting defaults.  Pleas"
	  "e report this bug to the Kokkos developers.  Original exception mess"
	  "age: " << e.what();
	throw std::logic_error(os.str());
      }
      numPartitions_ = numPartitions;

      // Recreate the cache blocking strategy.
      typedef CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_type;
      strategy_ = strategy_type (cacheSizeHint, sizeOfScalar);

      // Save the input parameter list.
      setMyParamList (plist);
    }

    /// \brief Default valid parameter list.
    ///
    /// The returned list contains all parameters accepted by \c
    /// KokkosNodeTsqr, with their default values and documentation.
    /// This method is reentrant and should be thread safe.  
    ///
    /// \note This method creates a new parameter list each time it is
    ///   called.  This saves storage for the common case of \c
    ///   setParameterList() being called infrequently (since \c
    ///   setParameterList() calls this method once).  If you find
    ///   yourself calling \c setParameterList() often, you might want
    ///   to change the implementation of getValidParameters() to
    ///   store the valid parameter list as member data.  Calling \c
    ///   setParameterList() often would be unusual for a class like
    ///   this one, whose configuration options are parameters related
    ///   to hardware that are unlikely to change at run time.
    Teuchos::RCP<const Teuchos::ParameterList> 
    getValidParameters() const
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = parameterList ("Intranode TSQR");
      params->set ("Cache Size Hint", 
		   static_cast<size_t>(0), 
		   std::string("Cache size in bytes; a hint for TSQR.  Set to t"
			       "he size of the largest private cache per CPU co"
			       "re, or the fraction of shared cache per core.  "
			       "If zero, we pick a reasonable default."));
      params->set ("Size of Scalar", 
		   sizeof(Scalar),
		   std::string ("Size in bytes of the Scalar type.  In most "
				"cases, the default sizeof(Scalar) is fine.  "
				"Set a non-default value only when Scalar's "
				"data is dynamically allocated (such as for a "
				"type with precision variable at run time)."));

      // The number of partitions is an int rather than a
      // LocalOrdinal, to ensure that it is always stored with the
      // same type, despite the type of LocalOrdinal.  Besides, Kokkos
      // wants an int anyway.
      params->set ("Num Partitions",
		   defaultNumPartitions (),
		   std::string ("Number of partitions; the maximum available pa"
				"rallelelism in intranode TSQR.  Slight oversub"
				"scription is OK; undersubscription may have a "
				"performance cost."));
      return params;
    }

    //! Factor the matrix A (see \c NodeTsqr documentation).
    FactorOutput
    factor (const LocalOrdinal numRows, 
	    const LocalOrdinal numCols,
	    Scalar A[],
	    const LocalOrdinal lda,
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguousCacheBlocks) const
    {
      MatView<LocalOrdinal, Scalar> A_view (numRows, numCols, A, lda);
      MatView<LocalOrdinal, Scalar> R_view (numCols, numCols, R, ldr);
      return factorImpl (A_view, R_view, contiguousCacheBlocks);
    }
    

    //! Apply the Q factor to C (see \c NodeTsqr documentation).
    void
    apply (const ApplyType& applyType,
	   const LocalOrdinal nrows,
	   const LocalOrdinal ncols_Q,
	   const Scalar Q[],
	   const LocalOrdinal ldq,
	   const FactorOutput& factorOutput,
	   const LocalOrdinal ncols_C,
	   Scalar C[],
	   const LocalOrdinal ldc,
	   const bool contiguousCacheBlocks) const
    {
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;

      const_view_type Q_view (nrows, ncols_Q, Q, ldq);
      view_type C_view (nrows, ncols_C, C, ldc);
      applyImpl (applyType, Q_view, factorOutput, C_view,
		 false, contiguousCacheBlocks);
    }

    //! Compute the explicit Q factor (see \c NodeTsqr documentation).
    void
    explicit_Q (const LocalOrdinal nrows,
		const LocalOrdinal ncols_Q,
		const Scalar Q[],
		const LocalOrdinal ldq,
		const FactorOutput& factorOutput,
		const LocalOrdinal ncols_C,
		Scalar C[],
		const LocalOrdinal ldc,
		const bool contiguousCacheBlocks) const
    {
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;

      const_view_type Q_view (nrows, ncols_Q, Q, ldq);
      view_type C_view (nrows, ncols_C, C, ldc);
      applyImpl (ApplyType::NoTranspose, Q_view, factorOutput,
		 C_view, true, contiguousCacheBlocks);
    }

    /// \brief Whether the R factor always has a nonnegative diagonal.
    /// 
    /// See the \c NodeTsqr documentation.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return combine_.QR_produces_R_factor_with_nonnegative_diagonal ();
    }

    /// \brief Cache size hint in bytes.
    ///
    /// This method is deprecated, because the name is misleading (the
    /// return value is not the size of one cache block, even though
    /// it is used to pick the "typical" cache block size).  Use \c
    /// cache_size_hint() instead (which returns the same value).
    /// 
    /// See the \c NodeTsqr documentation for details.
    size_t cache_block_size() const {
      return strategy_.cache_size_hint();
    }

    /// \brief Cache size hint in bytes.
    /// 
    /// See the \c NodeTsqr documentation for details.
    size_t cache_size_hint() const {
      return strategy_.cache_size_hint();
    }


    //! Fill A with zeros (see \c NodeTsqr documentation).
    void
    fill_with_zeros (const LocalOrdinal nrows,
		     const LocalOrdinal ncols,
		     Scalar A[],
		     const LocalOrdinal lda, 
		     const bool contiguousCacheBlocks) const
    {
      typedef MatView<LocalOrdinal, Scalar> view_type;
      view_type A_view (nrows, ncols, A, lda);

      typedef details::FillWDP<LocalOrdinal, Scalar> fill_wdp_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      fill_wdp_type filler (A_view, strategy_, STS::zero(), 
			    numPartitions_, contiguousCacheBlocks);
      node_->parallel_for (0, numPartitions_, filler);
    }

    //! Cache block A (see \c NodeTsqr documentation).
    void
    cache_block (const LocalOrdinal nrows,
		 const LocalOrdinal ncols, 
		 Scalar A_out[],
		 const Scalar A_in[],
		 const LocalOrdinal lda_in) const
    {
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      const_view_type A_in_view (nrows, ncols, A_in, lda_in);

      // The leading dimension of A_out doesn't matter here, since its
      // cache blocks are to be stored contiguously.  We set it
      // arbitrarily to a sensible value.
      typedef MatView<LocalOrdinal, Scalar> view_type;
      view_type A_out_view (nrows, ncols, A_out, nrows);

      typedef details::CacheBlockWDP<LocalOrdinal, Scalar> cb_wdp_type;
      cb_wdp_type cacheBlocker (A_in_view, A_out_view, strategy_, 
				numPartitions_, false);
      node_->parallel_for (0, numPartitions_, cacheBlocker);
    }

    //! Un - cache block A (see \c NodeTsqr documentation).
    void
    un_cache_block (const LocalOrdinal nrows,
		    const LocalOrdinal ncols,
		    Scalar A_out[],
		    const LocalOrdinal lda_out,		    
		    const Scalar A_in[]) const
    {
      // The leading dimension of A_in doesn't matter here, since its
      // cache blocks are contiguously stored.  We set it arbitrarily
      // to a sensible value.
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      const_view_type A_in_view (nrows, ncols, A_in, nrows);

      typedef MatView<LocalOrdinal, Scalar> view_type;
      view_type A_out_view (nrows, ncols, A_out, lda_out);

      typedef details::CacheBlockWDP<LocalOrdinal, Scalar> cb_wdp_type;
      cb_wdp_type cacheBlocker (A_in_view, A_out_view, strategy_, 
				numPartitions_, true);
      node_->parallel_for (0, numPartitions_, cacheBlocker);
    }

    //! Compute Q := Q*B in place (see \c NodeTsqr documentation).
    void
    Q_times_B (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       Scalar Q[],
	       const LocalOrdinal ldq,
	       const Scalar B[],
	       const LocalOrdinal ldb,
	       const bool contiguousCacheBlocks) const
    {
      typedef MatView<LocalOrdinal, Scalar> view_type;
      view_type Q_view (nrows, ncols, Q, ldq);

      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      const_view_type B_view (ncols, ncols, B, ldb);

      typedef details::MultWDP<LocalOrdinal, Scalar> mult_wdp_type;
      mult_wdp_type mult (Q_view, B_view, strategy_, numPartitions_, 
			  contiguousCacheBlocks);
      node_->parallel_for (0, numPartitions_, mult);
    }

  private:

    FactorOutput
    factorImpl (MatView<LocalOrdinal, Scalar> A,
		MatView<LocalOrdinal, Scalar> R,
		const bool contiguousCacheBlocks) const
    {
      if (A.empty())
	{
	  TEST_FOR_EXCEPTION(! R.empty(), std::logic_error,
			     "KokkosNodeTsqr::factorImpl: A is empty, but R "
			     "is not.  Please report this bug to the Kokkos "
			     "developers.");
	  return FactorOutput (0, 0);
	}
      const LocalOrdinal numRowsPerCacheBlock = 
	strategy_.cache_block_num_rows (A.ncols());
      const LocalOrdinal numCacheBlocks = 
	strategy_.num_cache_blocks (A.nrows(), A.ncols(), numRowsPerCacheBlock);
      //
      // Compute the first factorization pass (over partitions).
      //
      FactorOutput result (numCacheBlocks, numPartitions_);
      typedef details::FactorFirstPass<LocalOrdinal, Scalar> first_pass_type;
      first_pass_type firstPass (A, result.firstPassTauArrays, 
				 result.topBlocks, strategy_, 
				 numPartitions_, contiguousCacheBlocks);
      // parallel_for wants an exclusive range.
      node_->parallel_for (0, numPartitions_, firstPass);

      // Each partition collected a view of its top block, where that
      // partition's R factor is stored.  The second pass reduces
      // those R factors.  We do this on one thread to avoid the
      // overhead of parallelizing it.  If the typical use case is
      // oversubscription, you should parallelize this step with
      // multiple passes.  Note that we can't use parallel_reduce,
      // because the tree topology matters.
      factorSecondPass (result.topBlocks, result.secondPassTauArrays, 
			numPartitions_);

      // The "topmost top block" contains the resulting R factor.
      typedef MatView<LocalOrdinal, Scalar> view_type;
      const view_type& R_top = result.topBlocks[0];
      TEST_FOR_EXCEPTION(R_top.empty(), std::logic_error,
			 "After factorSecondPass: result.topBlocks[0] is an "
			 "empty MatView.  Please report this bug to the Kokkos "
			 "developers.");
      view_type R_top_square (R_top.ncols(), R_top.ncols(), 
			      R_top.get(), R_top.lda());
      R.fill (Teuchos::ScalarTraits<Scalar>::zero());
      // Only copy the upper triangle of R_top into R.
      copy_upper_triangle (R.ncols(), R.ncols(), R.get(), R.lda(),
			   R_top.get(), R_top.lda());
      return result;
    }

    void
    applyImpl (const ApplyType& applyType, 
	       const ConstMatView<LocalOrdinal, Scalar>& Q, 
	       const FactorOutput& factorOutput,
	       const MatView<LocalOrdinal, Scalar>& C, 
	       const bool explicitQ,
	       const bool contiguousCacheBlocks) const
    {
      using details::cacheBlockIndexRange;
      typedef details::ApplyFirstPass<LocalOrdinal, Scalar> first_pass_type;
      typedef CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;

      TEST_FOR_EXCEPTION(numPartitions_ != factorOutput.numPartitions(),
			 std::invalid_argument,
			 "applyImpl: KokkosNodeTsqr's number of partitions " 
			 << numPartitions_ << " does not match the given "
			 "factorOutput's number of partitions " 
			 << factorOutput.numPartitions() << ".  This likely "
			 "means that the given factorOutput object comes from "
			 "a different instance of KokkosNodeTsqr.  Please "
			 "report this bug to the Kokkos developers.");
      const int numParts = numPartitions_;
      first_pass_type firstPass (applyType, Q, factorOutput.firstPassTauArrays,
				 factorOutput.topBlocks, C, strategy_,
				 numParts, explicitQ, contiguousCacheBlocks);
      // Get a view of each partition's top block of the C matrix.
      std::vector<view_type> topBlocksOfC (numParts);
      {
	typedef std::pair<LocalOrdinal, LocalOrdinal> index_range_type;
	typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;
	blocker_type C_blocker (C.nrows(), C.ncols(), strategy_);

	// For each partition, collect its top block of C.
	for (int partIdx = 0; partIdx < numParts; ++partIdx)
	  {
	    const index_range_type cbIndices = 
	      cacheBlockIndexRange (C.nrows(), C.ncols(), partIdx, 
				    numParts, strategy_);
	    if (cbIndices.first >= cbIndices.second)
	      topBlocksOfC[partIdx] = view_type (0, 0, NULL, 0);
	    else
	      topBlocksOfC[partIdx] = 
		C_blocker.get_cache_block (C, cbIndices.first, 
					   contiguousCacheBlocks);
	  }
      }

      if (applyType.transposed())
	{
	  // parallel_for wants an exclusive range.
	  node_->parallel_for (0, numPartitions_, firstPass);
	  applySecondPass (applyType, factorOutput, topBlocksOfC, 
			   strategy_, explicitQ);
	}
      else
	{
	  applySecondPass (applyType, factorOutput, topBlocksOfC, 
			   strategy_, explicitQ);
	  // parallel_for wants an exclusive range.
	  node_->parallel_for (0, numPartitions_, firstPass);
	}
    }

    std::vector<Scalar>
    factorPair (const MatView<LocalOrdinal, Scalar>& R_top,
		const MatView<LocalOrdinal, Scalar>& R_bot) const
    {
      TEST_FOR_EXCEPTION(R_top.empty(), std::logic_error,
			 "R_top is empty!");
      TEST_FOR_EXCEPTION(R_bot.empty(), std::logic_error,
			 "R_bot is empty!");
      TEST_FOR_EXCEPTION(work_.size() == 0, std::logic_error,
			 "Workspace array work_ has length zero.");
      TEST_FOR_EXCEPTION(work_.size() < static_cast<size_t> (R_top.ncols()),
			 std::logic_error,
			 "Workspace array work_ has length = " 
			 << work_.size() << " < R_top.ncols() = " 
			 << R_top.ncols() << ".");

      std::vector<Scalar> tau (R_top.ncols());

      // Our convention for such helper methods is for the immediate
      // parent to allocate workspace (the work_ array in this case).
      //
      // The statement below only works if R_top and R_bot have a
      // nonzero (and the same) number of columns, but we have already
      // checked that above.
      combine_.factor_pair (R_top.ncols(), R_top.get(), R_top.lda(),
			    R_bot.get(), R_bot.lda(), &tau[0], &work_[0]);
      return tau;
    }

    void
    factorSecondPass (std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks,
		      std::vector<std::vector<Scalar> >& tauArrays,
		      const int numPartitions) const
    {
      if (numPartitions <= 1)
	return; // Done!
      TEST_FOR_EXCEPTION (topBlocks.size() < static_cast<size_t>(numPartitions), 
			  std::logic_error,
			  "KokkosNodeTsqr::factorSecondPass: topBlocks.size() "
			  "(= " << topBlocks.size() << ") < numPartitions (= " 
			  << numPartitions << ").  Please report this bug to "
			  "the Kokkos developers.");
      TEST_FOR_EXCEPTION (tauArrays.size() < static_cast<size_t>(numPartitions-1), 
			  std::logic_error,
			  "KokkosNodeTsqr::factorSecondPass: topBlocks.size() "
			  "(= " << topBlocks.size() << ") < numPartitions-1 (= " 
			  << (numPartitions-1) << ").  Please report this bug "
			  "to the Kokkos developers.");
      // The top partition (partition index zero) should always be
      // nonempty if we get this far, so its top block should also be
      // nonempty.
      TEST_FOR_EXCEPTION(topBlocks[0].empty(), std::logic_error,
			 "KokkosNodeTsqr::factorSecondPass: topBlocks[0] is "
			 "empty.  Please report this bug to the Kokkos "
			 "developers.");
      // However, other partitions besides the top one might be empty,
      // in which case their top blocks will be empty.  We skip over
      // the empty partitions in the loop below.
      work_.resize (static_cast<size_t> (topBlocks[0].ncols()));
      for (int partIdx = 1; partIdx < numPartitions; ++partIdx)
	if (! topBlocks[partIdx].empty())
	  tauArrays[partIdx-1] = factorPair (topBlocks[0], topBlocks[partIdx]);
    }

    void
    applyPair (const ApplyType& applyType,
	       const MatView<LocalOrdinal, Scalar>& R_bot,
	       const std::vector<Scalar>& tau,
	       const MatView<LocalOrdinal, Scalar>& C_top,
	       const MatView<LocalOrdinal, Scalar>& C_bot) const
    {
      // Our convention for such helper methods is for the immediate
      // parent to allocate workspace (the work_ array in this case).
      //
      // The statement below only works if C_top, R_bot, and C_bot
      // have a nonzero (and the same) number of columns, but we have
      // already checked that above.
      combine_.apply_pair (applyType, C_top.ncols(), R_bot.ncols(), 
			   R_bot.get(), R_bot.lda(), &tau[0], 
			   C_top.get(), C_top.lda(),
			   C_bot.get(), C_bot.lda(), &work_[0]);
    }

    void
    applySecondPass (const ApplyType& applyType,
		     const FactorOutput& factorOutput,
		     std::vector<MatView<LocalOrdinal, Scalar> >& topBlocksOfC,
		     const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
		     const bool explicitQ) const
    {
      const int numParts = factorOutput.numPartitions();
      if (numParts <= 1)
	return; // Done!
      TEST_FOR_EXCEPTION(topBlocksOfC.size() != static_cast<size_t>(numParts),
			 std::logic_error,
			 "KokkosNodeTsqr:applySecondPass: topBlocksOfC.size() ("
			 "= " << topBlocksOfC.size() << ") != number of partiti"
			 "ons (= " << numParts << ").  Please report this bug t"
			 "o the Kokkos developers.");
      TEST_FOR_EXCEPTION(factorOutput.secondPassTauArrays.size() != static_cast<size_t>(numParts-1),
			 std::logic_error,
			 "KokkosNodeTsqr:applySecondPass: factorOutput.secondPassTauArrays.size() ("
			 "= " << factorOutput.secondPassTauArrays.size() << ") != number of partiti"
			 "ons minus 1 (= " << (numParts-1) << ").  Please report this bug t"
			 "o the Kokkos developers.");

      const LocalOrdinal numCols = topBlocksOfC[0].ncols();
      work_.resize (static_cast<size_t> (numCols));

      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef typename std::vector<view_type>::size_type size_type;

      // Top blocks of C are the whole cache blocks.  We only want to
      // affect the top ncols x ncols part of each of those blocks in
      // this method.
      view_type C_top_square (numCols, numCols, topBlocksOfC[0].get(), 
			      topBlocksOfC[0].lda());
      if (applyType.transposed())
	{
	  // Don't include the topmost (index 0) partition in the
	  // iteration; that corresponds to C_top_square.
	  for (int partIdx = 1; partIdx < numParts; ++partIdx)
	    {
	      // It's legitimate for some partitions not to have any
	      // cache blocks.  In that case, their top block will be
	      // empty, and we can skip over them.
	      const view_type& C_cur = topBlocksOfC[partIdx];
	      if (! C_cur.empty())
		{
		  view_type C_cur_square (numCols, numCols, C_cur.get(), 
					  C_cur.lda());
		  // If explicitQ: We've already done the first pass and
		  // filled the top blocks of C.
		  applyPair (applyType, factorOutput.topBlocks[partIdx],
			     factorOutput.secondPassTauArrays[partIdx-1], 
			     C_top_square, C_cur_square);
		}
	    }
	}
      else
	{
	  // In non-transposed mode, when computing the first
	  // C.ncols() columns of the explicit Q factor, intranode
	  // TSQR would run after internode TSQR (i.e., DistTsqr)
	  // (even if only running on a single node in non-MPI mode).
	  // Therefore, internode TSQR is responsible for filling the
	  // top block of this node's part of the C matrix.
	  //
	  // Don't include the topmost partition in the iteration;
	  // that corresponds to C_top_square.
	  for (int partIdx = numParts - 1; partIdx > 0; --partIdx)
	    {
	      // It's legitimate for some partitions not to have any
	      // cache blocks.  In that case, their top block will be
	      // empty, and we can skip over them.
	      const view_type& C_cur = topBlocksOfC[partIdx];
	      if (! C_cur.empty())
		{
		  view_type C_cur_square (numCols, numCols, 
					  C_cur.get(), C_cur.lda());
		  // The "first" pass (actually the last, only named
		  // "first" by analogy with factorFirstPass()) will
		  // fill the rest of these top blocks.  For now, we
		  // just fill the top n x n part of the top blocks
		  // with zeros.
		  if (explicitQ)
		    C_cur_square.fill (Teuchos::ScalarTraits<Scalar>::zero());
		  applyPair (applyType, factorOutput.topBlocks[partIdx],
			     factorOutput.secondPassTauArrays[partIdx-1],
			     C_top_square, C_cur_square);
		}
	    }
	}
    }

  protected:

    /// \brief Return the topmost cache block of the matrix C.
    ///
    /// NodeTsqr's top_block() method must be implemented using its
    /// subclasses' const_top_block() method.  This is because
    /// top_block() is a template method, and template methods cannot
    /// be virtual.
    ///
    /// \param C [in] View of a matrix, with at least as many rows as
    ///   columns.
    /// \param contiguous_cache_blocks [in] Whether the cache blocks
    ///   of C are stored contiguously.
    ///
    /// \return View of the topmost cache block of the matrix C.
    ConstMatView<LocalOrdinal, Scalar>
    const_top_block (const ConstMatView<LocalOrdinal, Scalar>& C, 
		     const bool contiguous_cache_blocks) const 
    {
      typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;
      blocker_type blocker (C.nrows(), C.ncols(), strategy_);

      // C_top_block is a view of the topmost cache block of C.
      // C_top_block should have >= ncols rows, otherwise either cache
      // blocking is broken or the input matrix C itself had fewer
      // rows than columns.
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      const_view_type C_top = blocker.top_block (C, contiguous_cache_blocks);
      return C_top;
    }
  };
} // namespace TSQR

#endif // __TSQR_KokkosNodeTsqr_hpp
