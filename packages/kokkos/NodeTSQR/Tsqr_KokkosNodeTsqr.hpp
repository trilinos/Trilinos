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

#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace TSQR {

  namespace details {

    /// \brief Half-exclusive range of my cache block indices.
    ///
    /// FactorFirstPass breaks up the matrix into contiguous
    /// partitions of row blocks.  The index argument of Kokkos'
    /// parallel_for is the (zero-based) partition index.
    ///
    /// \param partitionIndex [in] Zero-based index of the
    ///   partition.  This is an int and not a LocalOrdinal, because
    ///   partition indices are arguments to Kokkos Node API
    ///   methods.  That is OK because the number of partitions
    ///   shouldn't be too much more than the amount of available
    ///   parallelism on this node.  Cache block indices are of
    ///   LocalOrdinal type, so don't mix them up with partition
    ///   indices.
    ///
    /// \return (start cache block index, end cache block index).
    ///   This is a half-exclusive range not including the end
    ///   point.  Thus, if the two indices are equal, the range is
    ///   empty.
    template<class LocalOrdinal, Scalar>
    std::pair<LocalOrdinal, LocalOrdinal>
    cacheBlockIndexRange (const LocalOrdinal numRows,
			  const LocalOrdinal numCols,
			  const int partitionIndex,
			  const int numPartitions,
			  const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy)
    {
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

      // Figure out how many cache blocks my partition contains.  If
      // the number of partitions doesn't evenly divide the number
      // of cache blocks, we spread out the remainder among the
      // first few threads.
      const LocalOrdinal quotient = numCacheBlocks / numPartitions;
      const LocalOrdinal remainder = numCacheBlocks % numPartitions;
      const LocalOrdinal myNumCacheBlocks = 
	(partitionIndex < remainder) ? (quotient + 1) : quotient;

      // If there are no cache blocks, there is nothing to factor.
      // Return an empty cache block range to indicate this.
      if (myNumCacheBlocks == 0)
	return std::make_pair (LocalOrdinal(0), LocalOrdinal(0));

      // Index of my first cache block (inclusive).
      myFirstCacheBlockIndex = 
	(partitionIndex < remainder) ? 
	partitionIndex * (quotient+1) :
	remainder * (quotient+1) + (partitionIndex - remainder) * quotient;
      // Index of my last cache block (exclusive).
      myLastCacheBlockIndex = 
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

    /// \class CacheBlockRange
    /// \brief "Iterator" over continuous range of cache blocks.
    /// \author Mark Hoemmen
    ///
    /// It's a pain to implement an actual C++ iterator (though Boost
    /// can help with that a little).  Thus, I've implemented a simple
    /// "sequence" over a contiguous range of cache blocks.  This is
    /// useful for \c KokkosNodeTsqr, in particular for \c
    /// FactorFirstPass.  Sequential TSQR's factorization is forward
    /// iteration over a sequence cache blocks, and applying the Q
    /// factor or computing the explicit Q factor is reverse
    /// iteration.
    ///
    /// This class is templated so that it works with a MatView or a
    /// ConstMatView.
    template<class MatrixViewType>
    class CacheBlockRange {
    public:
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      typedef typename MatrixViewType::scalar_type scalar_type;

      /// \brief Constructor
      /// 
      /// \param A [in] View of the matrix to factor.
      /// \param strategy [in] Cache blocking strategy to use (copied
      ///   on input).
      /// \param startIndex [in] Starting index of the cache block
      ///   sequence.
      /// \param endIndex [in] Ending index (exclusive) of the cache
      ///   block sequence.  Precondition: startIndex <= endIndex.  If
      ///   equal, the sequence is empty.
      /// \param reverse [in] false for forward iteration, true for
      ///   reverse iteration.
      /// \param contiguousCacheBlocks [in] Whether cache blocks in A 
      ///   are stored contiguously.
      CacheBlockRange (MatrixViewType A,
		       const CacheBlockingStrategy<ordinal_type, scalar_type>& strategy,
		       const ordinal_type startIndex,
		       const ordinal_type endIndex,
		       const bool reverse,
		       const bool contiguousCacheBlocks=false) :
	A_ (A), 
	curIndex_ (reverse ? (endIndex-1) : startIndex),
	endIndex_ (endIndex),
	blocker_ (A.nrows(), A.ncols(), strategy), 
	reverse_ (reverse),
	contiguousCacheBlocks_ (contiguousCacheBlocks)
      {}

      //! View the current cache block in the sequence.
      MatrixViewType 
      operator*() const 
      {
	if (reverse_)
	  {
	    if (curIndex_ < startIndex_)
	      return MatrixViewType (0, 0, NULL, 0); // empty
	    else
	      return blocker_.cache_block (A_, curIndex_, contiguousCacheBlocks_);
	  }
	else
	  {
	    if (curIndex_ >= endIndex_)
	      return MatrixViewType (0, 0, NULL, 0); // empty
	    else
	      return blocker_.cache_block (A_, curIndex_, contiguousCacheBlocks_);
	  }
      }

      //! Advance the sequence, if possible.  Otherwise do nothing.
      void
      advance ()
      {
	if (reverse_)
	  {
	    if (curIndex_ >= startIndex_)
	      curIndex_--;
	  }
	else
	  {
	    if (curIndex_ < endIndex_)
	      curIndex_++;
	  }
      }

      //! Prefix operator; advance the sequence, if possible.
      CacheBlockRange<MatrixViewType>& 
      operator++()
      {
	advance();
	return *this;
      }

      /// \brief Postfix operator; advance the sequence, if possible.
      ///
      /// \note The prefix operator++ is likely more efficient, since
      ///   it does not copy any data.
      CacheBlockRange<MatrixViewType>
      operator++(int)
      {
	CacheBlockRange<MatrixViewType> retval = *this;
	++(*this);
	return *retval;
      }

    private:
      MatrixViewType A_;
      ordinal_type curIndex_, endIndex_;
      CacheBlocker<ordinal_type, scalar_type> blocker_;
      bool reverse_, contiguousCacheBlocks_;
    };

    template<class LocalOrdinal, class Scalar>
    class FactorFirstPass {
    private:
      MatView<LocalOrdinal, Scalar> A_;
      std::vector<Scalar> work_;
      std::vector<std::vector<Scalar> >& tauArrays_;
      std::vector<MatView<LocalOrdinal, Scalar> > topBlocks_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      std::vector<Scalar>
      factorFirstCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
			     const MatView<LocalOrdinal, Scalar>& A_top)
      {
	std::vector<Scalar> tau (A_top.ncols());
	if (A_top.ncols() > 0)
	  combine.factor_first (A_top.nrows(), A_top.ncols(), A_top.get(), 
				A_top.lda(), &tau[0], &work_[0]);
	return tau;
      }

      std::vector<Scalar>
      factorCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
			const MatView<LocalOrdinal, Scalar>& A_top, 
			const MatView<LocalOrdinal, Scalar>& A_cur)
      {
	std::vector<Scalar> tau (A_top.ncols());
	if (A_top.ncols() > 0)
	  combine.factor_inner (A_cur.nrows(), A_top.ncols(), 
				A_top.get(), A_top.lda(),
				A_cur.get(), A_cur.lda(),
				&tau[0], &work_[0]);
	return tau;
      }

      /// \fn factor
      /// \brief Factor the given cache block range using SequentialTsqr.
      ///
      /// \param cbRange [in/out] Range of cache blocks to factor.
      /// \param tauIter [out] Output iterator of std::vector<Scalar>.
      ///   The "TAU" arrays from each cache block factored go here,
      ///   so this output sequence should be at least as long as the
      ///   number of cache blocks in the range.
      ///
      /// \return A view of the top block of the cache block range.
      template<class TauArraysOuterIterType>
      MatrixViewType
      factor (CacheBlockRange<MatView<LocalOrdinal, Scalar> >& cbRange,
	      TauArraysOuterIterType tauIter)
      {
	typedef MatView<LocalOrdinal, Scalar> view_type;
	Combine<LocalOrdinal, Scalar> combine;

	// Remember the top (first) block.
	view_type A_top = *cbRange;
	if (A_top.empty())
	  return A_top;

	// Factor the first block.
	*tauIter++ = factorFirstCacheBlock (combine, A_top);

	// Factor the remaining cache block(s).  CacheBlockRange's
	// prefix operator++ is more efficient than its postfix
	// operator++, since the latter has to make a copy of the
	// CacheBlockRange.
	++cbRange;
	view_type A_cur = *cbRange;
	while (! A_cur.empty())
	  {
	    *tauIter++ = factorCacheBlock (combine, A_top, A_cur);
	    ++cbRange;
	    A_cur = *cbRange;
	  }
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
	work_ (A.ncols()),
	tauArrays_ (tauArrays), 
	topBlocks_ (topBlocks),
	strategy_ (strategy), 
	numPartitions_ (numPartitions),
	contiguousCacheBlocks_ (contiguousCacheBlocks)
      {}

      /// \brief First pass of intranode TSQR factorization.
      ///
      /// Invoked by Kokkos' parallel_for template method.  This
      /// routine parallelizes over contiguous partitions of the
      /// matrix.  Each partition in turn contains cache blocks;
      /// partitions do not break up cache blocks.  (This ensures that
      /// the cache blocking scheme is the same as SequentialTsqr
      /// uses.  However, the implicit Q factor is not compatible with
      /// that of SequentialTsqr.)
      ///
      /// This method also saves a view of the top block of the
      /// partition.  This is useful for the next factorization pass.
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
      operator execute (const int partitionIndex) 
      {
	using Teuchos::ArrayView;
	using Teuchos::arrayViewFromVector;
	using std::vector;

	std::pair<LocalOrdinal, LocalOrdinal> cbIndices = 
	  cacheBlockIndexRange (A_.nrows(), A_.ncols(), partitionIndex, 
				numPartitions_, strategy_);
	if (A_.empty() || cbIndices.second <= cbIndices.first)
	  return;

	typedef MatView<LocalOrdinal, Scalar> view_type;
	typedef CacheBlockRange<view_type> range_type;
	range_type range (A_, strategy_, cbIndices.first, cbIndices.second,
			  false, contiguousCacheBlocks_);
	// We could just pass in &tauArrays_[cbIndices.first] (a raw
	// pointer) for the output iterator over the TAU arrays.
	// However, we would like a little bit more safety than that.
	// Teuchos::ArrayView uses ArrayRCP as the iterator when
	// HAVE_TEUCHOS_ARRAY_BOUNDSCHECK is defined, and uses a raw
	// pointer otherwise.  ArrayRCP is not thread safe, but that
	// is OK here, since the ArrayRCP is not being shared among
	// threads (we are inside a single thread's "execute"
	// routine).
	ArrayView<vector<Scalar> > tauView = 
	  arrayViewFromVector (tauArrays_).view (cbIndices.first, 
						 cbIndices.second - cbIndices.first);
	topBlocks_[partitionIndex] = factor (range, tauView.begin());
      }
    };


    /// \class ApplyFirstPass
    /// \brief "First" (actually last) pass of applying implicit Q factor.
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class ApplyFirstPass {
    private:
      ApplyType applyType_;
      ConstMatView<LocalOrdinal, Scalar> Q_;
      const std::vector<std::vector<Scalar> >& tauArrays_;
      const std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks_;
      MatView<LocalOrdinal, Scalar> C_;
      std::vector<Scalar> work_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool explicitQ_, contiguousCacheBlocks_;

      std::vector<Scalar>
      applyFirstCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
			    const ApplyType& applyType,
			    const ConstMatView<LocalOrdinal, Scalar>& Q_top,
			    const std::vector<Scalar>& tau,
			    const MatView<LocalOrdinal, Scalar>& C_top) const
      {
	if (Q_top.ncols() > 0 && C_top.ncols() > 0)
	combine.apply_first (applyType, C_top.nrows(), C_top.ncols(),
			     Q_top.ncols(), Q_top.get(), Q_top.lda(),
			     &tau[0], C_top.get(), C_top.lda(), &work_[0]);
      }

      void
      applyCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
		       const ApplyType& applyType,
		       const ConstMatView<LocalOrdinal, Scalar>& Q_cur,
		       const std::vector<Scalar>& tau,
		       const MatView<LocalOrdinal, Scalar>& C_top,
		       const MatView<LocalOrdinal, Scalar>& C_cur)
      {
	if (Q_cur.ncols() > 0 && C_cur.ncols() > 0)
	  combine.apply_inner (applyType, C_cur.nrows(), C_cur.ncols(), 
			       Q_cur.ncols(), Q_cur.get(), Q_cur.lda(), 
			       &tau[0],
			       C_top.get(), C_top.lda(),
			       C_cur.get(), C_cur.lda(),
			       &work_[0]);
      }

      /// \fn apply
      /// \brief Apply the sequential part of the implicit Q factor to C.
      ///
      /// \param applyType [in] Whether we are applying Q, Q^T, or Q^H.
      /// \param Q_range [in] Range of cache blocks in Q factor; reverse
      ///   range if not applyType.transposed(), else forward range.
      /// \param C_range [in/out] Range of cache blocks in C; reverse
      ///   range if not applyType.transposed(), else forward range.
      /// \param C_top [in/out] Topmost (last) block in C_range.
      /// \param tauIter [in] Forward iterator of std::vector<Scalar>;
      ///   the "TAU" arrays corresponding to the cache blocks in
      ///   Q_range.
      template<class TauArraysFwdIterType>
      void
      apply (const ApplyType& applyType,
	     CacheBlockRange<ConstMatView<LocalOrdinal, Scalar> >& Q_range,
	     CacheBlockRange<MatView<LocalOrdinal, Scalar> >& C_range,
	     const MatView<LocalOrdinal, Scalar>& C_top,
	     TauArraysFwdIterType tauIter)
      {
	typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
	typedef MatView<LocalOrdinal, Scalar> view_type;

	if (C_top.empty())
	  return;

	Combine<LocalOrdinal, Scalar> combine;
	if (applyType.transposed())
	  {
	    if (explicitQ_)
	      { 
		C_top.fill (Scalar(0));
		if (partitionIndex == 0)
		  for (LocalOrdinal j = 0; j < C_top.ncols(); ++j)
		    C_top(j,j) = Scalar(1);
	      }
	    // Apply the first block.
	    const_view_type Q_cur = *Q_range;
	    view_type C_cur = *C_range;
	    applyFirstCacheBlock (combine, applyType, Q_cur, *tauIter++, C_cur);

	    // Apply the rest of the blocks, if any.
	    ++Q_range;
	    ++C_range;
	    Q_cur = *Q_range;
	    C_cur = *C_range;
	    while (! C_cur.empty())
	      {
		if (explicitQ_)
		  C_cur.fill (Scalar(0));
		applyCacheBlock (combine, applyType, Q_cur, *tauIter++, 
				 C_top, C_cur);
		++Q_range;
		++C_range;
		Q_cur = *Q_range;
		C_cur = *C_range;
	      }
	  }
	else
	  {
	    if (explicitQ_)
	      {
		// We've already filled the top ncols x ncols block of C_top.
		view_type C_top_rest (C_top.nrows() - C_top.ncols(), C_top.ncols(),
				      C_top.get() + C_top.ncols(), C_top.lda());
		C_top_rest.fill (Scalar(0));
	      }
	    // Fetch the current element, advance to the next.
	    const_view_type Q_cur = *Q_range;
	    view_type C_cur = *C_range;
	    ++Q_range;
	    ++C_range;
	    // If the range is empty, C_top should be an empty block.
	    while (C_cur != C_top)
	      {
		if (explicitQ_)
		  C_cur.fill (Scalar(0));

		applyCacheBlock (combine, applyType, Q_cur, *tauIter++, 
				 C_top, C_cur);
		Q_cur = *Q_range;
		C_cur = *C_range;
		++Q_range;
		++C_range;
	      }
	    // Apply the first block.
	    applyFirstCacheBlock (combine, applyType, Q_cur, *tauIter++, C_cur);
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
		      std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks,
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
	work_ (Q_.ncols()),
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
      /// \param partitionIndex [in] Zero-based index of the
      ///   partition.  If greater than or equal to the number of
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
      operator execute (const int partitionIndex) 
      {
	// We use the same cache block indices for Q and for C.
	std::pair<LocalOrdinal, LocalOrdinal> cbIndices = 
	  cacheBlockIndexRange (Q_.nrows(), Q_.ncols(), partitionIndex, 
				numPartitions_, strategy_);
	if (Q_.empty() || cbIndices.second <= cbIndices.first)
	  return;

	typedef CacheBlockRange<ConstMatView<LocalOrdinal, Scalar> > const_range_type;
	typedef CacheBlockRange<MatView<LocalOrdinal, Scalar> > range_type;
	const bool reversed = applyType_.transposed() ? false : true;
	const_range_type Q_range (Q_, strategy_, 
				  cbIndices.first, cbIndices.second,
				  reversed, contiguousCacheBlocks_);
	range_type C_range (C_, strategy_, 
			    cbIndices.first, cbIndices.second,
			    reversed, contiguousCacheBlocks_);
	apply (applyType, Q_range, C_range, C_top, 
	       reversed ? &tauArrays_[cbIndices.second-1] : &tauArrays_[cbIndices.first]);
      }
    };


    /// \class CacheBlockWDP
    /// \brief Work-data pair (WDP) for (un-)cache-blocking a matrix.
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class CacheBlockWDP {
    private:
      ConstMatView<LocalOrdinal, Scalar> A_in_;
      MatView<LocalOrdinal, Scalar> A_out_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool unblock_;

      /// \fn copyRange
      /// \brief Copy one block range into another.
      ///
      /// \param cbInputRange [in] Range of input cache blocks.
      /// \param cbOutputRange [out] Range of output cache blocks.
      void
      copyRange (CacheBlockRange<ConstMatView<LocalOrdinal, Scalar> >& cbInputRange,
		 CacheBlockRange<MatView<LocalOrdinal, Scalar> >& cbOutputRange)
      {
	typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
	typedef MatView<LocalOrdinal, Scalar> view_type;

	const_view_type A_in_cur = *cbInputRange;
	view_type A_out_cur = *cbOutputRange;
	while (! A_in_cur.empty() && ! A_out_cur.empty())
	  {
	    A_out_cur.copy (A_in_cur);
	    ++cbInputRange;
	    ++cbOutputRange;
	    A_in_cur = *cbInputRange;
	    A_out_cur = *cbOutputRange;
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
      {}

      /// \brief Method called by Kokkos' parallel_for.
      ///
      /// \param partitionIndex [in] Zero-based index of the partition
      ///   of the matrix.  We parallelize over partitions.
      ///   Partitions respect cache blocks.
      void 
      operator execute (const int partitionIndex) 
      {
	typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
	typedef MatView<LocalOrdinal, Scalar> view_type;
	typedef CacheBlockRange<const_view_type> const_range_type;
	typedef CacheBlockRange<view_type> range_type;
	
	std::pair<LocalOrdinal, LocalOrdinal> cbIndices = 
	  cacheBlockIndexRange (A_.nrows(), A_.ncols(), partitionIndex, 
				numPartitions_, strategy_);
	if (A_in_.empty() || A_out_.empty() || 
	    cbIndices.second <= cbIndices.first)
	  return;

	// If unblock_ is false, then A_in_ is in column-major order,
	// and we want to cache-block it into A_out_.  If unblock_ is
	// true, then A_in_ is cache-blocked, and we want to
	// un-cache-block it into A_out_ (a matrix in column-major
	// order).
	const_range_type inputRange (A_in_, strategy_, cbIndices.first, 
				     cbIndices.second, false, unblock_);
	range_type outputRange (A_out_, strategy_, cbIndices.first, 
				cbIndices.second, false, ! unblock_);
	copyRange (inputRange, outputRange);
      }
    };

    /// \class MultWDP
    /// \brief Kokkos work-data pair for \c KokkosNodeTsqr::Q_times_B().
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class MultWDP {
    private:
      MatView<LocalOrdinal, Scalar> Q_;
      ConstMatView<LocalOrdinal, Scalar> B_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      void
      multBlock (BLAS<LocalOrdinal, Scalar>& blas,
		 const MatView<LocalOrdinal, Scalar>& Q_cur,
		 const Matrix<LocalOrdinal, Scalar>& Q_temp)
      {
	// GEMM doesn't like aliased arguments, so we use a copy.  We
	// only copy the current cache block, rather than all of Q;
	// this saves memory.
	Q_temp.reshape (Q_cur.nrows(), Q_cur.ncols());
	Q_temp.copy (Q_cur);

	// Q_cur := Q_temp * B.
	blas.GEMM ("N", "N", Q_cur.nrows(), ncols, ncols, Scalar(1),
		   Q_temp.get(), Q_temp.lda(), B_.get(), B_.lda(),
		   Scalar(0), Q_cur.get(), Q_cur.lda());
      }

      /// \brief Multiply (in place) each cache block in the range by B_.
      ///
      /// \param cbRange [in/out] Range of cache blocks.
      void
      multRange (CacheBlockRange<MatView<LocalOrdinal, Scalar> >& cbRange)
      {
	typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
	typedef MatView<LocalOrdinal, Scalar> view_type;

	view_type Q_cur = *cbRange;
	if (Q_cur.empty())
	  return;
	else {
	  const LocalOrdinal numCols = Q_cur.ncols();
	  Matrix<LocalOrdinal, Scalar> Q_temp;
	  BLAS<LocalOrdinal, Scalar> blas;
	  do {
	    multBlock (blas, Q_cur, Q_temp);
	    ++cbRange;
	    Q_cur = *cbRange;
	  } while (! Q_cur.empty());
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
      operator execute (const int partitionIndex) 
      {
	typedef MatView<LocalOrdinal, Scalar> view_type;
	typedef CacheBlockRange<view_type> range_type;
	
	std::pair<LocalOrdinal, LocalOrdinal> cbIndices = 
	  cacheBlockIndexRange (A_.nrows(), A_.ncols(), partitionIndex, 
				numPartitions_, strategy_);
	if (Q_.empty() || cbIndices.second <= cbIndices.first)
	  return;

	range_type range (Q_, strategy_, cbIndices.first, 
			  cbIndices.second, false, contiguousCacheBlocks_);
	multRange (range);
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
      firstPassTauArrays (numCacheBlocks),
      secondPassTauArrays (numPartitions),
      topBlocks (numPartitions)
    {}

    //! Total number of cache blocks in the matrix (over all partitions).
    int numCacheBlocks() const { return firstPassTauArrays.size(); }

    //! Number of partitions of the matrix; max available parallelism.
    int numPartitions() const { return secondPassTauArrays.size(); }

    //! TAU arrays from the first pass; one per cache block.
    std::vector<std::vector<Scalar> > firstPassTauArrays;

    /// \brief TAU arrays from the second pass; one per partition.
    ///
    /// For now, KokkosNodeTsqr::factor() uses only two passes over
    /// the matrix.  firstPassTauArrays contains the result of the
    /// pass over cache blocks, and secondPassTauArrays contains the
    /// result of combining the upper triangular R factors from the
    /// first pass.  Later, we may add more passes, in which case we
    /// will likely combine firstPassTauArrays and secondPassTauArrays
    /// into a single std::vector or Teuchos::Tuple.
    std::vector<std::vector<Scalar> > secondPassTauArrays;

    //! Views of the topmost cache blocks in each partition.
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
    //! Size in bytes of a cache block.
    size_t cacheBlockSize_;

    //! Size in bytes of the Scalar type.
    size_t sizeOfScalar_;

    /// \brief Number of partitions; max available parallelism.
    ///
    /// The number of partitions is an int rather than a LocalOrdinal,
    /// to ensure that it is always stored in the ParameterList with
    /// the same type, despite the type of LocalOrdinal.  Besides,
    /// Kokkos wants an int anyway.
    int numPartitions_;

    /// \brief Default number of partitions.
    ///
    /// \warning Do not expect the returned value to be the same, even
    ///   at run time on subsequent invocations.  This method may in
    ///   the future try to "learn" the optimal number of partitions.
    int 
    defaultNumPartitions (const Teuchos::RCP<const node_type>& node) const 
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
	 << ">: cacheBlockSize=" << cacheBlockSize_
	 << ", sizeOfScalar=" << sizeOfScalar_
	 << ", numPartitions=" << numPartitions_;
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
      using Teuchos::RCP;
      RCP<ParameterList> plist;
      if (paramList.is_null())
	plist = getValidParameters ();
      else
	{
	  plist = paramList;
	  plist->validateParametersAndSetDefaults (getValidParameters());
	}
      // Get values of parameters.  We do this "transactionally" so
      // that (except for validation and filling in defaults above)
      // this method has the strong exception guarantee (it either
      // returns, or throws an exception with no externally visible
      // side effects).
      size_t cacheBlockSize, sizeOfScalar;
      int numPartitions;
      try {
	cacheBlockSize = plist->get<size_t> ("Cache Block Size");
	sizeOfScalar = plist->get<size_t> ("Size of Scalar");
	numPartitions = plist.get<int> ("Num Partitions");
      } catch (Teuchos::Exceptions::InvalidParameter& e) {
	std::ostringstream os;
	os << "Failed to read default parameters after setting defaults.  Pleas"
	  "e report this bug to the Kokkos developers.  Original exception mess"
	  "age: " << e.what();
	throw std::logic_error(os.str());
      }
      cacheBlockSize_ = cacheBlockSize;
      sizeOfScalar_ = sizeOfScalar;
      numPartitions_ = numPartitions;

      // Save the input parameter list.
      //
      // ??? Do we need this cast when the base class is not a
      // template class?
      typedef ParameterListAcceptorDefaultBase pla_type;
      static_cast<pla_type*>(this)->paramList_ = plist;
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

      RCP<ParameterList> params = parameterList ("Node TSQR");
      params->set ("Cache Block Size", 
		   static_cast<size_t>(0), 
		   std::string("Size in bytes of the cache blocks used by TSQR."
			       "  Set to the size of the largest private cache "
			       "per CPU core, or the fraction of shared cache p"
			       "er core.  Slight underestimation is OK."));
      params->set ("Size of Scalar", 
		   sizeof(Scalar),
		   std::string ("Size of the Scalar type.  In most cases, the d"
				"efault sizeof(Scalar) is fine.  Set a non-defa"
				"ult value only when Scalar is an variable-prec"
				"ision type whose data is dynamically allocated"
				"."));
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
      return factorImpl (MatView<LocalOrdinal, Scalar> (numRows, numCols, A, lda),
			 MatView<LocalOrdinal, Scalar> (numCols, numCols, R, ldr),
			 contiguousCacheBlocks);
    }
    

    //! Apply the Q factor to C (see \c NodeTsqr documentation).
    void
    apply (const ApplyType& applyType,
	   const Ordinal nrows,
	   const Ordinal ncols_Q,
	   const Scalar Q[],
	   const Ordinal ldq,
	   const FactorOutput& factorOutput,
	   const Ordinal ncols_C,
	   Scalar C[],
	   const Ordinal ldc,
	   const bool contiguousCacheBlocks) const
    {
      ConstMatView<LocalOrdinal, Scalar> Q_view (nrows, ncols_Q, Q, ldq);
      MatView<LocalOrdinal, Scalar> C_view (nrows, ncols_C, C, ldc);
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

      applyImpl (ApplyType::NoTranspose, 
		 const_view_type (nrows, ncols_Q, Q, ldq), 
		 factorOutput,
		 view_type (nrows, ncols_C, C, ldc), 
		 true, 
		 contiguousCacheBlocks);
    }

    /// \brief Whether the R factor always has a nonnegative diagonal.
    /// 
    /// See the \c NodeTsqr documentation.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return combine_.QR_produces_R_factor_with_nonnegative_diagonal ();
    }

    //! Fill A with zeros (see \c NodeTsqr documentation).
    void
    fill_with_zeros (const Ordinal nrows,
		     const Ordinal ncols,
		     Scalar A[],
		     const Ordinal lda, 
		     const bool contiguousCacheBlocks) const;

    //! Cache block A (see \c NodeTsqr documentation).
    void
    cache_block (const Ordinal nrows,
		 const Ordinal ncols, 
		 Scalar A_out[],
		 const Scalar A_in[],
		 const Ordinal lda_in) const
    {
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;

      const_view_type A_in_view (nrows, ncols, A_in, lda_in);
      view_type A_out_view (nrows, ncols, A_out, lda_out);
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy (cacheBlockSize_, sizeOfScalar_);
      CacheBlockWDP cacheBlocker (A_in_view, A_out_view, strategy, numPartitions_, false);

      pNode_->parallel_for (0, numPartitions_, cacheBlocker);
    }

    //! Un - cache block A (see \c NodeTsqr documentation).
    void
    un_cache_block (const Ordinal nrows,
		    const Ordinal ncols,
		    Scalar A_out[],
		    const Ordinal lda_out,		    
		    const Scalar A_in[]) const
    {
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_type;
      typedef CacheBlockWDP<LocalOrdinal, Scalar> cb_wdp_type;

      const_view_type A_in_view (nrows, ncols, A_in, lda_in);
      view_type A_out_view (nrows, ncols, A_out, lda_out);
      strategy_type strategy (cacheBlockSize_, sizeOfScalar_);
      cb_wdp_type cacheBlocker (A_in_view, A_out_view, strategy, 
				numPartitions_, true);
      pNode_->parallel_for (0, numPartitions_, cacheBlocker);
    }

    //! Compute Q := Q*B in place (see \c NodeTsqr documentation).
    void
    Q_times_B (const Ordinal nrows,
	       const Ordinal ncols,
	       Scalar Q[],
	       const Ordinal ldq,
	       const Scalar B[],
	       const Ordinal ldb,
	       const bool contiguousCacheBlocks) const
    {
      typedef ConstMatView<LocalOrdinal, Scalar> const_view_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_type;
      typedef MultWDP<LocalOrdinal, Scalar> mult_type;

      view_type Q_view (nrows, ncols, Q, ldq);
      const_view_type B_view (ncols, ncols, B, ldb);
      strategy_type strategy (cacheBlockSize_, sizeOfScalar_);
      mult_type mult (Q_view, B_view, strategy, numPartitions_, 
		      contiguousCacheBlocks);
      pNode_->parallel_for (0, numPartitions_, mult);
    }

    //! Return the topmost cache block of C (see \c NodeTsqr documentation).
    MatView<Ordinal, Scalar>
    top_block (const MatView<Ordinal, Scalar>& C,
	       const bool contiguousCacheBlocks) const
    {
      typedef CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;
      typedef const std::pair<LocalOrdinal, LocalOrdinal> range_type;
      typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;

      strategy_type strategy (cacheBlockSize_, sizeOfScalar_);
      const range_type range = 
	cacheBlockIndexRange (C.nrows(), C.ncols(), int(0),
			      numPartitions_, strategy);
      blocker_type C_blocker (C.nrows(), C.ncols(), strategy);
      if (range.first == range.second)
	return view_type (0, 0, NULL, 0); // empty cache block
      else
	return C_blocker.get_cache_block (C, range.first, 
					  contiguousCacheBlocks);
    }

  private:

    FactorOutput
    factorImpl (const MatView<LocalOrdinal, Scalar>& A,
		const MatView<LocalOrdinal, Scalar>& R,
		const bool contiguousCacheBlocks) const
    {
      using details::FactorFirstPass;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy (cacheBlockSize_, 
							    sizeOfScalar_);
      const LocalOrdinal numRowsPerCacheBlock = 
	strategy.cache_block_num_rows (A.ncols());
      const LocalOrdinal numCacheBlocks = 
	strategy.num_cache_blocks (A.nrows(), A.ncols(), numRowsPerCacheBlock);
      //
      // Compute the first factorization pass (over partitions).
      //
      FactorOutput result (numCacheBlocks, numPartitions);
      FactorFirstPass firstPass (A, result.firstPassTauArrays, 
				 result.topBlocks, strategy, 
				 numPartitions_, contiguousCacheBlocks);
      // parallel_for wants an exclusive range.
      pNode_->parallel_for (0, numPartitions_, firstPass);

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
      R.copy (result.topBlocks[0]);
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
      using details::ApplyFirstPass;
      typedef ApplyFirstPass<LocalOrdinal, Scalar> first_pass_type;
      typedef CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;

      strategy_type strategy (cacheBlockSize_, sizeOfScalar_);
      first_pass_type firstPass (applyType, Q, factorOutput.firstPassTauArrays,
				 factorOutput.topBlocks, C, strategy,
				 factorOutput.numPartitions(),
				 explicitQ, contiguousCacheBlocks);
      std::vector<view_type> topBlocksOfC (factorOutput.numPartitions());
      {
	typedef const std::pair<LocalOrdinal, LocalOrdinal> range_type;
	typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;
	const range_type range = 
	  cacheBlockIndexRange (C.nrows(), C.ncols(), partitionIndex,
				numPartitions_, strategy);
	blocker_type C_blocker (C.nrows(), C.ncols(), strategy);
	for (int partitionIndex = 0; 
	     partitionIndex < numPartitions_; 
	     ++partitionIndex)
	  {
	    if (range.first == range.second)
	      topBlocksOfC[partitionIndex] = view_type (0, 0, NULL, 0);
	    else
	      {
		topBlocksOfC[partitionIndex] = 
		  C_blocker.get_cache_block (C, range.first, 
					     contiguousCacheBlocks);
	      }
	  }
      }

      if (applyType.transposed())
	{
	  // parallel_for wants an exclusive range.
	  pNode_->parallel_for (0, numPartitions_, firstPass);
	  applySecondPass (applyType, factorOutput, topBlocksOfC, strategy, explicitQ_)
	}
      else
	{
	  applySecondPass (applyType, factorOutput, topBlocksOfC, strategy, explicitQ_);
	  // parallel_for wants an exclusive range.
	  pNode_->parallel_for (0, numPartitions_, firstPass);
	}
    }

    void
    factorPair (const MatView<LocalOrdinal, Scalar>& R_top,
		const MatView<LocalOrdinal, Scalar>& R_bot,
		std::vector<Scalar>& tau) const
    {
      // Our convention for such helper methods is for the immediate
      // parent to allocate workspace and other resources.
      if (R_top.ncols() > 0)
	combine_.factor_pair (R_top.ncols(), R_top.get(), R_top.lda(),
			      R_bot.get(), R_bot.lda(), &tau[0], &work_[0]);
    }

    void
    factorSecondPass (std::vector<MatView<LocalOrdinal, Scalar> >& topBlocks,
		      std::vector<std::vector<Scalar> > tauArrays,
		      const LocalOrdinal numPartitions) const
    {
      if (numPartitions <= 1)
	return; // Done!

      const LocalOrdinal numCols = topBlocks[0].ncols();
      typedef std::vector<MatView<LocalOrdinal, Scalar> >::size_type size_type;
      for (size_type k = 1; k < numPartitions; ++k)
	{
	  tau.resize (numCols);
	  factorPair (topBlocks[0], topBlocks[k], tauArrays[k]);
	}
    }

    void
    applyPair (const ApplyType& applyType,
	       const MatView<LocalOrdinal, Scalar>& R_bot,
	       const std::vector<Scalar>& tau,
	       const MatView<LocalOrdinal, Scalar>& C_top,
	       const MatView<LocalOrdinal, Scalar>& C_bot) const
    {
      if (R_top.ncols() > 0 && C_top.ncols() > 0)
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
      if (factorOutput.numPartitions() <= 1)
	return; // Done!
      const LocalOrdinal numCols = topBlocks[0].ncols();
      typedef std::vector<MatView<LocalOrdinal, Scalar> >::size_type size_type;
      typedef MatView<LocalOrdinal, Scalar> view_type;

      // Top blocks of C are the whole cache blocks.  We only want to
      // affect the top ncols x ncols part of each of those blocks in
      // this method.
      view_type C_top (topBlocksOfC[0].ncols(), topBlocksOfC[0].ncols(),
		       topBlocksOfC[0].get(), topBlocksOfC[0].lda());
      if (applyType.transposed())
	{
	  for (size_type k = 1; k < factorOutput.numPartitions(); ++k)
	    {
	      view_type C_cur (topBlocksOfC[k].ncols(), topBlocksOfC[k].ncols(),
			       topBlocksOfC[k].get(), topBlocksOfC[k].lda());
	      // If explicitQ: We've already done the first pass and
	      // filled the top blocks of C.
	      applyPair (applyType,
			 factorOutput.topBlocks[k],
			 factorOutput.secondPassTauArrays[k],
			 C_top, C_cur);
	    }
	}
      else
	{
	  // In non-transpose mode, when computing the first
	  // C.ncols() columns of the explicit Q factor, intranode
	  // TSQR would run after internode TSQR (i.e., DistTsqr)
	  // (even if only running on a single node in non-MPI
	  // mode).  Therefore, internode TSQR is responsible for
	  // filling the top block of this node's part of the C
	  // matrix.
	  for (size_type k = factorOutput.numPartitions-1; k >= 0; --k)
	    {
	      view_type C_cur (topBlocksOfC[k].ncols(), topBlocksOfC[k].ncols(),
			       topBlocksOfC[k].get(), topBlocksOfC[k].lda());
	      // The first pass will fill the rest of these top
	      // blocks; for now, we just fill the top n x n part of
	      // the top blocks with zeros.
	      if (explicitQ)
		C_cur.fill (Scalar(0));
	      applyPair (applyType,
			 factorOutput.topBlocks[k],
			 factorOutput.secondPassTauArrays[k],
			 C_top, C_cur);
	    }
	}
    }
  };
} // namespace TSQR

#endif // __TSQR_KokkosNodeTsqr_hpp
