//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

/// \file Tsqr_KokkosNodeTsqr.hpp
/// \brief Parallel intranode TSQR implemented using Kokkos::parallel_for.

#ifndef __TSQR_KokkosNodeTsqr_hpp
#define __TSQR_KokkosNodeTsqr_hpp

#include "Tsqr_CacheBlocker.hpp"
#include "Tsqr_Combine.hpp"
#include "Tsqr_NodeTsqr.hpp"
#include "Tsqr_Impl_SystemBlas.hpp"

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Kokkos_Core.hpp"

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
      using LO = LocalOrdinal;
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
      if (partitionIndex >= numPartitions) {
        return {0, 0};
      }

      const LO numRowsCacheBlock =
        strategy.cache_block_num_rows (numCols);
      const LO numCacheBlocks =
        strategy.num_cache_blocks (numRows, numCols, numRowsCacheBlock);

      // Figure out how many cache blocks my partition contains.  If
      // the number of partitions doesn't evenly divide the number
      // of cache blocks, we spread out the remainder among the
      // first few threads.
      const LO quotient = numCacheBlocks / numPartitions;
      const LO remainder = numCacheBlocks - quotient * numPartitions;
      const LO myNumCacheBlocks = (partitionIndex < remainder) ?
        (quotient + 1) : quotient;

      // If there are no cache blocks, there is nothing to factor.
      // Return an empty cache block range to indicate this.
      if (myNumCacheBlocks == 0) {
        return {0, 0};
      }

      // Index of my first cache block (inclusive).
      const LO myFirstCacheBlockIndex = (partitionIndex < remainder) ?
        partitionIndex * (quotient+1) :
        remainder * (quotient+1) + (partitionIndex - remainder) * quotient;
      // Index of my last cache block (exclusive).
      const LO myLastCacheBlockIndex = (partitionIndex+1 < remainder) ?
        (partitionIndex+1) * (quotient+1) :
        remainder * (quotient+1) + (partitionIndex+1 - remainder) * quotient;
      TEUCHOS_TEST_FOR_EXCEPTION
        (myLastCacheBlockIndex <= myFirstCacheBlockIndex,
         std::logic_error, "Partition " << (partitionIndex+1) << " of "
         << numPartitions << ":  My range of cache block indices ["
         << myFirstCacheBlockIndex << ", " << myLastCacheBlockIndex
         << ") is empty.");
      return {myFirstCacheBlockIndex, myLastCacheBlockIndex};
    }


    /// \class FactorFirstPass
    /// \brief First pass of KokkosNodeTsqr's factorization.
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class FactorFirstPass {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view_type;

    private:
      mat_view_type A_;
      // While tauArrays_ is shared among tasks (i.e., partitions),
      // there are no race conditions among entries, since each
      // partition writes its own entry.  Ditto for topBlocks_.
      std::vector<std::vector<Scalar> >& tauArrays_;
      std::vector<mat_view_type>& topBlocks_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      std::vector<Scalar>
      factorFirstCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
                             const mat_view_type& A_top,
                             std::vector<Scalar>& work) const
      {
        std::vector<Scalar> tau (A_top.extent(1));

        // We should only call this if A_top.extent(1) > 0 and therefore
        // work.size() > 0, but we've already checked for that, so we
        // don't have to check again.
        combine.factor_first (A_top, tau.data(), work.data());
        return tau;
      }

      std::vector<Scalar>
      factorCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
                        const mat_view_type& A_top,
                        const mat_view_type& A_cur,
                        std::vector<Scalar>& work) const
      {
        std::vector<Scalar> tau (A_top.extent(1));

        // We should only call this if A_top.extent(1) > 0 and therefore
        // tau.size() > 0 and work.size() > 0, but we've already
        // checked for that, so we don't have to check again.
        combine.factor_inner (A_top, A_cur, tau.data(), work.data());
        return tau;
      }

      /// \brief Factor the given cache block range using sequential TSQR.
      ///
      /// \param cbIndices [in] Half-exclusive range of cache block indices.
      /// \param partitionIndex [in] Zero-based index of my partition.
      ///
      /// \return A view of the top block of the cache block range.
      mat_view_type
      factor (const std::pair<LocalOrdinal, LocalOrdinal> cbIndices,
              const int partitionIndex) const
      {
        const char suffix[] = "  Please report this bug to the Tpetra developers.";
        using cb_range_type = CacheBlockRange<mat_view_type>;

        // Workspace is created here, because it must not be shared
        // among threads.
        std::vector<Scalar> work (A_.extent(1));

        // Range of cache blocks to factor.
        cb_range_type cbRange (A_, strategy_, cbIndices.first,
                               cbIndices.second, contiguousCacheBlocks_);
        // Iterator in the forward direction over the range of cache
        // blocks to factor.
        typedef typename CacheBlockRange<mat_view_type>::iterator range_iter_type;
        range_iter_type cbIter = cbRange.begin();

        // Remember the top (first) block.
        mat_view_type A_top = *cbIter;
        if (A_top.empty ()) {
          return A_top;
        }
        TEUCHOS_TEST_FOR_EXCEPTION
          (cbIndices.first >= cbIndices.second, std::logic_error,
           "FactorFirstPass::factor: A_top is not empty, but the "
           "cache block index range " << cbIndices.first << ","
           << cbIndices.second << " is empty." << suffix);

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
        while (cbIter != cbEnd) {
          mat_view_type A_cur = *cbIter;
          // Iteration over cache blocks of a partition should
          // always result in nonempty cache blocks.
          TEUCHOS_TEST_FOR_EXCEPTION
            (A_cur.empty (), std::logic_error, "FactorFirstPass::factor: "
             "The current cache block (the " << count << "-th to factor in the "
             "range [" << cbIndices.first << "," << cbIndices.second << ") of "
             "cache block indices) in partition " << (partitionIndex+1) << " "
             "(out of " << numPartitions_ << " partitions) is empty." << suffix);
          TEUCHOS_TEST_FOR_EXCEPTION
            (static_cast<size_t>(curTauIdx) >= tauArrays_.size(),
             std::logic_error, "FactorFirstPass::factor: curTauIdx (= "
             << curTauIdx << ") >= tauArrays_.size() (= "
             << tauArrays_.size() << ")." << suffix);
          tauArrays_[curTauIdx++] =
            factorCacheBlock (combine, A_top, A_cur, work);
          ++count;
          ++cbIter;
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
      FactorFirstPass (const mat_view_type& A,
                       std::vector<std::vector<Scalar> >& tauArrays,
                       std::vector<mat_view_type>& topBlocks,
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
        TEUCHOS_TEST_FOR_EXCEPTION(A_.empty(), std::logic_error,
                           "TSQR::FactorFirstPass constructor: A is empty.  "
                           "Please report this bug to the Kokkos developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(numPartitions < 1, std::logic_error,
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
      void operator() (const int partitionIndex) const
      {
        if (partitionIndex < 0 || partitionIndex >= numPartitions_ || A_.empty ()) {
          return;
        }
        else {
          const std::pair<LocalOrdinal, LocalOrdinal> cbIndices =
            cacheBlockIndexRange (A_.extent(0), A_.extent(1), partitionIndex,
                                  numPartitions_, strategy_);
          // It's legitimate, though suboptimal, for some partitions
          // not to get any work to do (in this case, not to get any
          // cache blocks to factor).
          if (cbIndices.second <= cbIndices.first) {
            return;
          } else {
            topBlocks_[partitionIndex] = factor (cbIndices, partitionIndex);
          }
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
    public:
      using const_mat_view_type = MatView<LocalOrdinal, const Scalar>;
      using mat_view_type = MatView<LocalOrdinal, Scalar>;

    private:
      ApplyType applyType_;
      const_mat_view_type Q_;
      const std::vector<std::vector<Scalar> >& tauArrays_;
      const std::vector<mat_view_type>& topBlocks_;
      mat_view_type C_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool explicitQ_, contiguousCacheBlocks_;

      void
      applyFirstCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
                            const ApplyType& applyType,
                            const const_mat_view_type& Q_top,
                            const std::vector<Scalar>& tau,
                            const mat_view_type& C_top,
                            std::vector<Scalar>& work) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(tau.size() < static_cast<size_t> (Q_top.extent(1)),
                           std::logic_error,
                           "ApplyFirstPass::applyFirstCacheBlock: tau.size() "
                           "(= " << tau.size() << ") < number of columns "
                           << Q_top.extent(1) << " in the Q factor.  Please "
                           "report this bug to the Kokkos developers.");

        // If we get this far, it's fair to assume that we have
        // checked whether tau and work have nonzero lengths.
        combine.apply_first (applyType, Q_top, tau.data(),
                             C_top, work.data());
      }

      void
      applyCacheBlock (Combine<LocalOrdinal, Scalar>& combine,
                       const ApplyType& applyType,
                       const const_mat_view_type& Q_cur,
                       const std::vector<Scalar>& tau,
                       const mat_view_type& C_top,
                       const mat_view_type& C_cur,
                       std::vector<Scalar>& work) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (tau.size() < static_cast<size_t> (Q_cur.extent(1)),
           std::logic_error, "ApplyFirstPass::applyCacheBlock: tau.size() "
           "(= " << tau.size() << ") < number of columns "
           << Q_cur.extent(1) << " in the Q factor."
           "  Please report this bug to the Tpetra developers.");

        // If we get this far, it's fair to assume that we have
        // checked whether tau and work have nonzero lengths.
        combine.apply_inner (applyType, C_cur.extent(0), C_cur.extent(1),
                             Q_cur.extent(1), Q_cur.data(), Q_cur.stride(1),
                             tau.data(),
                             C_top.data(), C_top.stride(1),
                             C_cur.data(), C_cur.stride(1),
                             work.data());
      }

      /// \fn apply
      /// \brief Apply the sequential part of the implicit Q factor to C.
      ///
      /// \param applyType [in] Whether we are applying Q, Q^T, or Q^H.
      /// \param cbIndices [in] Half-exclusive range of cache block
      ///   indices.
      /// \param partitionIndex [in] The argument to \c operator(); the
      ///   index of the partition which instance of ApplyFirstPass
      ///   is currently processing.
      void
      apply (const ApplyType& applyType,
             const std::pair<LocalOrdinal, LocalOrdinal> cbIndices,
             const int partitionIndex) const
      {
        using const_range_type = CacheBlockRange<const_mat_view_type>;
        using range_type = CacheBlockRange<mat_view_type>;
        const char suffix[] = "  Please report this bug to the Tpetra developers.";

        if (cbIndices.first >= cbIndices.second) {
          return; // My range of cache blocks is empty; nothing to do
        }

        // Q_range: Range of cache blocks in the Q factor.
        // C_range: Range of cache blocks in the matrix C.
        const_range_type Q_range (Q_, strategy_,
                                  cbIndices.first, cbIndices.second,
                                  contiguousCacheBlocks_);
        range_type C_range (C_, strategy_,
                            cbIndices.first, cbIndices.second,
                            contiguousCacheBlocks_);
        TEUCHOS_TEST_FOR_EXCEPTION
          (Q_range.empty(), std::logic_error,
           "Q_range is empty, but the range of cache block "
           "indices [" << cbIndices.first << ", "
           << cbIndices.second << ") is not empty." << suffix);
        TEUCHOS_TEST_FOR_EXCEPTION
          (C_range.empty(), std::logic_error,
           "C_range is empty, but the range of cache block "
           "indices [" << cbIndices.first << ", "
           << cbIndices.second << ") is not empty." << suffix);

        // Task-local workspace array of length C_.extent(1).  Workspace
        // must be per task, else there will be race conditions as
        // different tasks attempt to write to and read from the same
        // workspace simultaneously.
        std::vector<Scalar> work (C_.extent(1));

        Combine<LocalOrdinal, Scalar> combine;
        if (applyType.transposed ()) {
          auto Q_rangeIter = Q_range.begin();
          auto C_rangeIter = C_range.begin();
          TEUCHOS_TEST_FOR_EXCEPTION
            (Q_rangeIter == Q_range.end(), std::logic_error,
             "The Q cache block range claims to be nonempty, "
             "but the iterator range is empty." << suffix);
          TEUCHOS_TEST_FOR_EXCEPTION
            (C_rangeIter == C_range.end(), std::logic_error,
             "The C cache block range claims to be nonempty, "
             "but the iterator range is empty." << suffix);

          // Q_top: Topmost cache block in the cache block range of Q.
          // C_top: Topmost cache block in the cache block range of C.
          const_mat_view_type Q_top = *Q_rangeIter;
          mat_view_type C_top = *C_rangeIter;
          if (explicitQ_) {
            deep_copy (C_top, Scalar {});
            if (partitionIndex == 0) {
              for (LocalOrdinal j = 0; j < C_top.extent(1); ++j) {
                C_top(j,j) = Scalar (1.0);
              }
            }
          }
          LocalOrdinal curTauIndex = cbIndices.first;

          // Apply the first block.
          applyFirstCacheBlock (combine, applyType, Q_top,
                                tauArrays_[curTauIndex++], C_top, work);

          // Apply the rest of the blocks, if any.
          ++Q_rangeIter;
          ++C_rangeIter;
          while (Q_rangeIter != Q_range.end ()) {
            TEUCHOS_TEST_FOR_EXCEPTION
              (C_rangeIter == C_range.end(), std::logic_error,
               "When applying Q^T or Q^H to C: The Q cache "
               "block iterator is not yet at the end, but "
               "the C cache block iterator is." << suffix);
            const_mat_view_type Q_cur = *Q_rangeIter;
            mat_view_type C_cur = *C_rangeIter;
            ++Q_rangeIter;
            ++C_rangeIter;
            if (explicitQ_) {
              deep_copy (C_cur, Scalar {});
            }
            applyCacheBlock (combine, applyType, Q_cur,
                             tauArrays_[curTauIndex++],
                             C_top, C_cur, work);
          }
        }
        else {
          // Q_top: Topmost cache block in the cache block range of Q.
          // C_top: Topmost cache block in the cache block range of C.
          const_mat_view_type Q_top = *(Q_range.begin());
          mat_view_type C_top = *(C_range.begin());

          if (explicitQ_) {
            // We've already filled the top ncols x ncols block of
            // C_top with data (that's the result of applying the
            // internode part of the Q factor via DistTsqr).  However,
            // we still need to fill the rest of C_top (everything but
            // the top ncols rows of C_top) with zeros.
            mat_view_type C_top_rest (C_top.extent(0) - C_top.extent(1),
                                      C_top.extent(1),
                                      C_top.data() + C_top.extent(1),
                                      C_top.stride(1));
            deep_copy (C_top_rest, Scalar {});
          }
          LocalOrdinal curTauIndex = cbIndices.second-1;

          // When applying Q (rather than Q^T or Q^H), we apply the
          // cache blocks in reverse order.
          typename const_range_type::iterator Q_rangeIter = Q_range.rbegin();
          typename range_type::iterator C_rangeIter = C_range.rbegin();
          TEUCHOS_TEST_FOR_EXCEPTION
            (Q_rangeIter == Q_range.rend(), std::logic_error,
             "The Q cache block range claims to be nonempty, "
             "but the iterator range is empty." << suffix);
          TEUCHOS_TEST_FOR_EXCEPTION
            (C_rangeIter == C_range.rend(), std::logic_error,
             "The C cache block range claims to be nonempty, "
             "but the iterator range is empty." << suffix);

          // Equality of cache block range iterators only tests the
          // cache block index, not reverse-ness.  This means we can
          // compare a reverse-direction iterator (Q_rangeIter) with
          // a forward-direction iterator (Q_range.begin()).
          //
          // We do this because we need to handle the topmost block
          // of Q_range separately (applyFirstCacheBlock(), rather
          // than applyCacheBlock()).
          while (Q_rangeIter != Q_range.begin ()) {
            const_mat_view_type Q_cur = *Q_rangeIter;
            mat_view_type C_cur = *C_rangeIter;

            if (explicitQ_) {
              deep_copy (C_cur, Scalar {});
            }
            TEUCHOS_TEST_FOR_EXCEPTION
              (curTauIndex < cbIndices.first, std::logic_error,
               "curTauIndex=" << curTauIndex << " out of valid "
               "range [" << cbIndices.first << ","
               << cbIndices.second << ")." << suffix);
            applyCacheBlock (combine, applyType, Q_cur,
                             tauArrays_[curTauIndex--],
                             C_top, C_cur, work);
            ++Q_rangeIter;
            ++C_rangeIter;
          }
          TEUCHOS_TEST_FOR_EXCEPTION
            (curTauIndex < cbIndices.first, std::logic_error,
             "curTauIndex=" << curTauIndex << " out of valid range "
             "[" << cbIndices.first << "," << cbIndices.second << ")."
             << suffix);
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
                      const const_mat_view_type& Q,
                      const std::vector<std::vector<Scalar>>& tauArrays,
                      const std::vector<mat_view_type>& topBlocks,
                      const mat_view_type& C,
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
      void operator() (const int partitionIndex) const
      {
        const char prefix[] = "TSQR::ApplyFirstPass::operator(): ";
        const char suffix[] = "  Please report this bug to the Tpetra developers.";

        if (partitionIndex < 0 || partitionIndex >= numPartitions_ ||
            Q_.empty () || C_.empty ()) {
          return;
        }

        // We use the same cache block indices for Q and for C.
        std::pair<LocalOrdinal, LocalOrdinal> cbIndices =
          cacheBlockIndexRange (Q_.extent(0), Q_.extent(1), partitionIndex,
                                numPartitions_, strategy_);
        if (cbIndices.second <= cbIndices.first)
          return;
        {
          std::pair<size_t, size_t> cbInds (size_t (cbIndices.first),
                                            size_t (cbIndices.second));
          TEUCHOS_TEST_FOR_EXCEPTION
            (cbIndices.first < LocalOrdinal(0), std::logic_error,
             prefix << "cacheBlockIndexRange(" << Q_.extent (0) << ", "
             << Q_.extent(1) << ", " << partitionIndex << ", "
             << numPartitions_ << ", strategy) returned a cache block "
             "range " << cbIndices.first << "," << cbIndices.second <<
             " with negative starting index." << suffix);
          TEUCHOS_TEST_FOR_EXCEPTION
            (cbInds.second > tauArrays_.size (), std::logic_error,
             prefix << "cacheBlockIndexRange(" << Q_.extent (0) << ", "
             << Q_.extent(1) << ", " << partitionIndex << ", "
             << numPartitions_ << ", strategy) returned a cache block "
             "range" << cbIndices.first << "," << cbIndices.second <<
             " with starting index larger than the number of tau "
             "arrays " << tauArrays_.size () << "." << suffix);
        }
        apply (applyType_, cbIndices, partitionIndex);
      }
    };

    /// \class CacheBlockFunctor
    /// \brief Kokkos functor for KokkosNodeTsqr's (un_)cache_block() methods.
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class CacheBlockFunctor {
    private:
      using const_mat_view_type = MatView<LocalOrdinal, const Scalar>;
      using mat_view_type = MatView<LocalOrdinal, Scalar>;
      using const_range_type = CacheBlockRange<const_mat_view_type>;
      using range_type = CacheBlockRange<mat_view_type>;

      const_mat_view_type A_in_;
      mat_view_type A_out_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool unblock_;

      /// \brief Copy one range of cache blocks into another.
      ///
      /// \param cbInputRange [in] Range of input cache blocks.
      /// \param cbOutputRange [out] Range of output cache blocks.
      void copyRange (const_range_type& cbInputRange,
                      range_type& cbOutputRange) const
      {
        typedef typename const_range_type::iterator input_iter_type;
        typedef typename range_type::iterator output_iter_type;

        input_iter_type inputIter = cbInputRange.begin();
        output_iter_type outputIter = cbOutputRange.begin();

        input_iter_type inputEnd = cbInputRange.end();
        // TODO (mfh 29 Jun 2012) In a debug build, check in the loop
        // below whether outputIter == cbOutputRange.end().  If so,
        // throw std::logic_error.  Don't declare outputEnd unless
        // we're in a debug build, because otherwise the compiler may
        // report warnings (gcc 4.5 doesn't; gcc 4.6 does).
        // output_iter_type outputEnd = cbOutputRange.end();

        while (inputIter != inputEnd) {
          const_mat_view_type A_in_cur = *inputIter;
          mat_view_type A_out_cur = *outputIter;
          deep_copy (A_out_cur, A_in_cur);
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
      CacheBlockFunctor (const const_mat_view_type A_in,
                         const mat_view_type A_out,
                         const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
                         const int numPartitions,
                         const bool unblock) :
        A_in_ (A_in),
        A_out_ (A_out),
        strategy_ (strategy),
        numPartitions_ (numPartitions),
        unblock_ (unblock)
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (A_in_.extent(0) != A_out_.extent(0) ||
           A_in_.extent(1) != A_out_.extent(1),
           std::invalid_argument,
           "A_in and A_out do not have the same dimensions: "
           "A_in is " << A_in_.extent(0) << " by "
           << A_in_.extent(1) << ", but A_out is "
           << A_out_.extent(0) << " by "
           << A_out_.extent(1) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION
          (numPartitions_ < 1, std::invalid_argument,
           "The number of partitions " << numPartitions_
           << " is not a positive integer.");
      }

      /// \brief Method called by Kokkos::parallel_for.
      ///
      /// \param partitionIndex [in] Zero-based index of the partition
      ///   of the matrix.  We parallelize over partitions.
      ///   Partitions respect cache blocks.
      void operator() (const int partitionIndex) const
      {
        if (partitionIndex < 0 || partitionIndex >= numPartitions_ ||
            A_in_.empty()) {
          return;
        }
        else {
          using index_range_type = std::pair<LocalOrdinal, LocalOrdinal>;
          const index_range_type cbIndices =
            cacheBlockIndexRange (A_in_.extent (0), A_in_.extent (1),
                                  partitionIndex, numPartitions_, strategy_);
          // It's perfectly legal for a partitioning to assign zero
          // cache block indices to a particular partition.  In that
          // case, this task has nothing to do.
          if (cbIndices.first >= cbIndices.second) {
            return;
          }
          else {
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

    /// \class MultFunctor
    /// \brief Kokkos functor for \c KokkosNodeTsqr::Q_times_B().
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class MultFunctor {
    private:
      using const_mat_view_type = MatView<LocalOrdinal, const Scalar>;
      using mat_view_type = MatView<LocalOrdinal, Scalar>;
      using range_type = CacheBlockRange<mat_view_type>;

      mat_view_type Q_;
      const_mat_view_type B_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      // This uses SystemBlas for now.
      // In the future, we may want to use a TPL.
      // That means we could switch to RawBlas.
      void
      multBlock (Impl::SystemBlas<Scalar>& blas,
                 const mat_view_type& Q_cur,
                 Matrix<LocalOrdinal, Scalar>& Q_temp) const
      {
        using Teuchos::NO_TRANS;
        const LocalOrdinal numCols = Q_cur.extent (1);

        // GEMM doesn't like aliased arguments, so we use a copy.  We
        // only copy the current cache block, rather than all of Q;
        // this saves memory.
        Q_temp.reshape (Q_cur.extent (0), numCols);
        deep_copy (Q_temp, Q_cur);

        // Q_cur := Q_temp * B.
        blas.GEMM (NO_TRANS, NO_TRANS, Q_cur.extent(0), numCols, numCols,
                   Scalar (1.0),
                   Q_temp.data(), Q_temp.stride(1), B_.data(), B_.stride(1),
                   Scalar(0), Q_cur.data(), Q_cur.stride(1));
      }

      /// \brief Multiply (in place) each cache block in the range by B_.
      ///
      /// \param cbRange [in/out] Range of cache blocks.
      void multRange (range_type& cbRange) const
      {
        typedef typename range_type::iterator iter_type;
        iter_type iter = cbRange.begin();
        iter_type end = cbRange.end();

        // Temporary storage for the BLAS' matrix-matrix multiply
        // routine (which forbids aliasing of any input argument and
        // the output argument).
        Matrix<LocalOrdinal, Scalar> Q_temp;
        Impl::SystemBlas<Scalar> blas;
        while (iter != end) {
          mat_view_type Q_cur = *iter;
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
      MultFunctor (const mat_view_type Q,
                   const const_mat_view_type B,
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
      void operator() (const int partitionIndex) const
      {
        if (partitionIndex < 0 || partitionIndex >= numPartitions_ ||
            Q_.empty ()) {
          return;
        }
        else {
          typedef std::pair<LocalOrdinal, LocalOrdinal> index_range_type;
          const index_range_type cbIndices =
            cacheBlockIndexRange (Q_.extent (0), Q_.extent (1), partitionIndex,
                                  numPartitions_, strategy_);
          if (cbIndices.first >= cbIndices.second) {
            return;
          }
          else {
            range_type range (Q_, strategy_, cbIndices.first,
                              cbIndices.second, contiguousCacheBlocks_);
            multRange (range);
          }
        }
      }
    };

    /// \class FillFunctor
    /// \brief Kokkos functor for \c KokkosNodeTsqr::fill_with_zeros().
    /// \author Mark Hoemmen
    template<class LocalOrdinal, class Scalar>
    class FillFunctor {
    private:
      using mat_view_type = MatView<LocalOrdinal, Scalar>;
      using range_type = CacheBlockRange<mat_view_type>;

      mat_view_type A_;
      CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
      const Scalar value_;
      int numPartitions_;
      bool contiguousCacheBlocks_;

      //! Fill (in place) each cache block in the range with value.
      void fillRange (range_type& cbRange, const Scalar value) const
      {
        typedef typename range_type::iterator iter_type;
        iter_type iter = cbRange.begin();
        iter_type end = cbRange.end();
        while (iter != end) {
          mat_view_type A_cur = *iter;
          deep_copy (A_cur, value);
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
      FillFunctor (const mat_view_type A,
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
      void operator() (const int partitionIndex) const
      {
        if (partitionIndex < 0 || partitionIndex >= numPartitions_ ||
            A_.empty ()) {
          return;
        }
        else {
          typedef std::pair<LocalOrdinal, LocalOrdinal> index_range_type;
          const index_range_type cbIndices =
            cacheBlockIndexRange (A_.extent(0), A_.extent(1), partitionIndex,
                                  numPartitions_, strategy_);
          if (cbIndices.first >= cbIndices.second) {
            return;
          }
          else {
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
    typedef MatView<LocalOrdinal, Scalar> mat_view_type;

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
      TEUCHOS_TEST_FOR_EXCEPTION(theNumPartitions < 1, std::invalid_argument,
                         "TSQR::KokkosNodeTsqrFactorOutput: Invalid number of "
                         "partitions " << theNumPartitions << "; number of "
                         "partitions must be a positive integer.");
      // If there's only one partition, we don't even need a second
      // pass (it's just sequential TSQR), and we don't need a TAU
      // array for the top partition.
      secondPassTauArrays.resize (size_t (theNumPartitions-1));
      topBlocks.resize (size_t (theNumPartitions));
    }

    //! Total number of cache blocks in the matrix (over all partitions).
    int numCacheBlocks() const { return firstPassTauArrays.size(); }

    //! Number of partitions of the matrix; max available parallelism.
    int numPartitions() const { return topBlocks.size(); }

    //! TAU arrays from the first pass; one per cache block.
    std::vector<std::vector<Scalar>> firstPassTauArrays;

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
    std::vector<std::vector<Scalar>> secondPassTauArrays;

    /// \brief Views of the topmost cache blocks in each partition.
    ///
    /// One entry for each partition.
    std::vector<mat_view_type> topBlocks;
  };

  /// \class KokkosNodeTsqr
  /// \brief Intranode (within an MPI process) TSQR parallelized using
  ///   Kokkos::DefaultHostExecutionSpace.
  /// \author Mark Hoemmen
  ///
  /// \tparam LocalOrdinal The type of indices in the (node-local)
  ///   matrix.
  ///
  /// \tparam Scalar The type of entries in the (node-local) matrix.
  ///
  /// This implementation of the intranode part of TSQR factors the
  /// matrix in two passes.  The first pass parallelizes over
  /// partitions, doing Sequential TSQR over each partition.  The
  /// second pass combines the R factors from the partitions, and is
  /// not currently parallel.  Thus, the overall algorithm is similar
  /// to that of TbbTsqr, except that:
  /// <ul>
  /// <li> TbbTsqr partitions differently; KokkosNodeTsqr's partitions
  ///      use the same layout of cache blocks as SequentialTsqr,
  ///      whereas TbbTsqr uses a different layout. </li>
  /// <li> TbbTsqr reduces the R factors in parallel; it only needs
  ///      one "pass." </li>
  /// </ul>
  template<class LocalOrdinal, class Scalar>
  class KokkosNodeTsqr :
    public NodeTsqr<LocalOrdinal, Scalar, KokkosNodeTsqrFactorOutput<LocalOrdinal, Scalar>>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef Scalar scalar_type;

    using const_mat_view_type = MatView<LocalOrdinal, const Scalar>;
    using mat_view_type = MatView<LocalOrdinal, Scalar>;

    /// \typedef FactorOutput
    /// \brief Part of the implicit Q representation returned by factor().
    typedef typename NodeTsqr<LocalOrdinal, Scalar, KokkosNodeTsqrFactorOutput<LocalOrdinal, Scalar> >::factor_output_type FactorOutput;

    /// \brief Constructor (with user-specified parameters).
    ///
    /// \param params [in/out] List of parameters.  Missing parameters
    ///   will be filled in with default values.
    KokkosNodeTsqr (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
    {
      setParameterList (params);
    }

    /// \brief Whether this object is ready to perform computations.
    bool ready() const {
      return true;
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
         << ">: \"Cache Size Hint\"=" << strategy_.cache_size_hint()
         << ", \"Size of Scalar\"=" << strategy_.size_of_scalar()
         << ", \"Num Tasks\"=" << numPartitions_;
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
      if (paramList.is_null()) {
        plist = rcp (new ParameterList (*getValidParameters ()));
      }
      else {
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
        numPartitions = plist->get<int> ("Num Tasks");
      }
      catch (Teuchos::Exceptions::InvalidParameter& e) {
        std::ostringstream os;
        os << "Failed to read default parameters after setting defaults.  Pleas"
          "e report this bug to the Kokkos developers.  Original exception mess"
          "age: " << e.what();
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, os.str());
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
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters() const
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      if (defaultParams_.is_null()) {
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
        params->set ("Num Tasks",
                     defaultNumPartitions (),
                     std::string ("Number of partitions; the maximum available pa"
                                  "rallelelism in intranode TSQR.  Slight oversub"
                                  "scription is OK; undersubscription may have a "
                                  "performance cost."));
        defaultParams_ = params;
      }
      return defaultParams_;
    }

    FactorOutput
    factor (const LocalOrdinal numRows,
            const LocalOrdinal numCols,
            Scalar A[],
            const LocalOrdinal lda,
            Scalar R[],
            const LocalOrdinal ldr,
            const bool contiguousCacheBlocks) const
    {
      mat_view_type A_view (numRows, numCols, A, lda);
      mat_view_type R_view (numCols, numCols, R, ldr);
      return factorImpl (A_view, R_view, contiguousCacheBlocks);
    }

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
      const_mat_view_type Q_view (nrows, ncols_Q, Q, ldq);
      mat_view_type C_view (nrows, ncols_C, C, ldc);
      applyImpl (applyType, Q_view, factorOutput, C_view,
                 false, contiguousCacheBlocks);
    }

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
      const_mat_view_type Q_view (nrows, ncols_Q, Q, ldq);
      mat_view_type C_view (nrows, ncols_C, C, ldc);
      applyImpl (ApplyType::NoTranspose, Q_view, factorOutput,
                 C_view, true, contiguousCacheBlocks);
    }

    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return combine_.QR_produces_R_factor_with_nonnegative_diagonal ();
    }

    size_t cache_size_hint() const {
      return strategy_.cache_size_hint();
    }

    void
    fill_with_zeros (const LocalOrdinal nrows,
                     const LocalOrdinal ncols,
                     Scalar A[],
                     const LocalOrdinal lda,
                     const bool contiguousCacheBlocks) const
    {
      mat_view_type A_view (nrows, ncols, A, lda);

      using functor_type = details::FillFunctor<LocalOrdinal, Scalar>;
      const Scalar ZERO {};
      functor_type functor (A_view, strategy_, ZERO, numPartitions_,
                            contiguousCacheBlocks);
      using execution_space = Kokkos::DefaultHostExecutionSpace;
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<int>>
        range (0, numPartitions_);
      Kokkos::parallel_for ("KokkosNodeTsqr::fill_with_zeros", range, functor);
    }

    void
    cache_block (const LocalOrdinal nrows,
                 const LocalOrdinal ncols,
                 Scalar A_out[],
                 const Scalar A_in[],
                 const LocalOrdinal lda_in) const
    {
      const_mat_view_type A_in_view (nrows, ncols, A_in, lda_in);

      // The leading dimension of A_out doesn't matter here, since its
      // cache blocks are to be stored contiguously.  We set it
      // arbitrarily to a sensible value.
      mat_view_type A_out_view (nrows, ncols, A_out, nrows);

      using functor_type = details::CacheBlockFunctor<LocalOrdinal, Scalar>;
      functor_type functor (A_in_view, A_out_view, strategy_,
                            numPartitions_, false);
      using execution_space = Kokkos::DefaultHostExecutionSpace;
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<int>>
        range (0, numPartitions_);
      Kokkos::parallel_for ("KokkosNodeTsqr::cache_block", range, functor);
    }

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
      const_mat_view_type A_in_view (nrows, ncols, A_in, nrows);
      mat_view_type A_out_view (nrows, ncols, A_out, lda_out);

      using functor_type = details::CacheBlockFunctor<LocalOrdinal, Scalar>;
      functor_type functor (A_in_view, A_out_view, strategy_,
                            numPartitions_, true);
      using execution_space = Kokkos::DefaultHostExecutionSpace;
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<int>>
        range (0, numPartitions_);
      Kokkos::parallel_for ("KokkosNodeTsqr::un_cache_block", range, functor);
    }

    void
    Q_times_B (const LocalOrdinal nrows,
               const LocalOrdinal ncols,
               Scalar Q[],
               const LocalOrdinal ldq,
               const Scalar B[],
               const LocalOrdinal ldb,
               const bool contiguousCacheBlocks) const
    {
      mat_view_type Q_view (nrows, ncols, Q, ldq);
      const_mat_view_type B_view (ncols, ncols, B, ldb);

      using functor_type = details::MultFunctor<LocalOrdinal, Scalar>;
      functor_type functor (Q_view, B_view, strategy_, numPartitions_,
                            contiguousCacheBlocks);
      using execution_space = Kokkos::DefaultHostExecutionSpace;
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<int>>
        range (0, numPartitions_);
      Kokkos::parallel_for ("KokkosNodeTsqr::Q_times_B", range, functor);
    }

  private:
    //! Implementation of fundamental TSQR kernels.
    Combine<LocalOrdinal, Scalar> combine_;

    //! Workspace for Combine operations.
    mutable std::vector<Scalar> work_;

    //! Cache blocking strategy.
    CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;

    /// \brief Number of partitions; max available parallelism.
    ///
    /// The number of partitions is an int rather than a LocalOrdinal,
    /// to ensure that it is always stored in the ParameterList with
    /// the same type, despite the type of LocalOrdinal.  Besides,
    /// Kokkos wants an int anyway.
    int numPartitions_;

    //! Default parameter list (set by \c getValidParameters()).
    mutable Teuchos::RCP<const Teuchos::ParameterList> defaultParams_;

    //! Default number of partitions.
    int
    defaultNumPartitions () const
    {
      return Kokkos::DefaultHostExecutionSpace::concurrency ();
    }

    FactorOutput
    factorImpl (mat_view_type A,
                mat_view_type R,
                const bool contiguousCacheBlocks) const
    {
      const char prefix[] = "KokkosNodeTsqr::factorImpl: ";
      const char suffix[] = "  Please report this bug to the Tpetra developers.";
      using LO = LocalOrdinal;
      using execution_space = Kokkos::DefaultHostExecutionSpace;
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<int>>
        range (0, numPartitions_);

      if (A.empty ()) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (! R.empty (), std::logic_error, prefix << "A is empty, "
           "but R is not." << suffix);
        return FactorOutput (0, 0);
      }
      const LO numRowsPerCacheBlock =
        strategy_.cache_block_num_rows (A.extent(1));
      const LO numCacheBlocks =
        strategy_.num_cache_blocks (A.extent(0), A.extent(1), numRowsPerCacheBlock);
      //
      // Compute the first factorization pass (over partitions).
      //
      FactorOutput result (numCacheBlocks, numPartitions_);
      using first_pass_type = details::FactorFirstPass<LO, Scalar>;
      first_pass_type firstPass (A, result.firstPassTauArrays,
                                 result.topBlocks, strategy_,
                                 numPartitions_, contiguousCacheBlocks);
      Kokkos::parallel_for ("KokkosNodeTsqr::factorImpl::firstPass",
                            range, firstPass);

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
      const mat_view_type& R_top = result.topBlocks[0];
      TEUCHOS_TEST_FOR_EXCEPTION
        (R_top.empty (), std::logic_error, prefix << "After "
         "factorSecondPass: result.topBlocks[0] is an empty view."
         << suffix);
      mat_view_type R_top_square (R_top.extent(1), R_top.extent(1),
                                  R_top.data(), R_top.stride(1));
      deep_copy (R, Scalar {});
      // Only copy the upper triangle of R_top into R.
      copy_upper_triangle (R.extent(1), R.extent(1), R.data(), R.stride(1),
                           R_top.data(), R_top.stride(1));
      return result;
    }

    void
    applyImpl (const ApplyType& applyType,
               const const_mat_view_type& Q,
               const FactorOutput& factorOutput,
               const mat_view_type& C,
               const bool explicitQ,
               const bool contiguousCacheBlocks) const
    {
      const char prefix[] = "KokkosNodeTsqr::applyImpl: ";
      const char suffix[] = "  Please report this bug to the Tpetra developers.";
      using LO = LocalOrdinal;
      using details::cacheBlockIndexRange;
      using first_pass_type = details::ApplyFirstPass<LO, Scalar>;
      using execution_space = Kokkos::DefaultHostExecutionSpace;

      TEUCHOS_TEST_FOR_EXCEPTION
        (numPartitions_ != factorOutput.numPartitions(),
         std::invalid_argument, prefix << "KokkosNodeTsqr's number "
         "of partitions " << numPartitions_ << " does not match the "
         "given factorOutput's number of partitions "
         << factorOutput.numPartitions() << ".  This likely means "
         "that the given factorOutput object comes from a different "
         "instance of KokkosNodeTsqr." << suffix);
      const int numParts = numPartitions_;
      first_pass_type firstPass (applyType, Q,
                                 factorOutput.firstPassTauArrays,
                                 factorOutput.topBlocks, C, strategy_,
                                 numParts, explicitQ,
                                 contiguousCacheBlocks);
      // Get a view of each partition's top block of the C matrix.
      std::vector<mat_view_type> topBlocksOfC (numParts);
      {
        using index_range_type = std::pair<LO, LO>;
        using blocker_type = CacheBlocker<LO, Scalar>;
        blocker_type C_blocker (C.extent(0), C.extent(1), strategy_);

        // For each partition, collect its top block of C.
        for (int partIdx = 0; partIdx < numParts; ++partIdx) {
          const index_range_type cbIndices =
            cacheBlockIndexRange (C.extent(0), C.extent(1), partIdx,
                                  numParts, strategy_);
          if (cbIndices.first >= cbIndices.second) {
            topBlocksOfC[partIdx] = mat_view_type (0, 0, nullptr, 0);
          } else {
            topBlocksOfC[partIdx] =
              C_blocker.get_cache_block (C, cbIndices.first,
                                         contiguousCacheBlocks);
          }
        }
      }

      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<int>>
        range(0, numPartitions_);
      if (applyType.transposed ()) {
        Kokkos::parallel_for ("KokkosNodeTsqr::applyImpl::firstPass",
                              range, firstPass);
        applySecondPass (applyType, factorOutput, topBlocksOfC,
                         strategy_, explicitQ);
      }
      else {
        applySecondPass (applyType, factorOutput, topBlocksOfC,
                         strategy_, explicitQ);
        Kokkos::parallel_for ("KokkosNodeTsqr::applyImpl::firstPass",
                              range, firstPass);
      }
    }

    std::vector<Scalar>
    factorPair (const mat_view_type& R_top,
                const mat_view_type& R_bot) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (R_top.empty (), std::logic_error, "R_top is empty!");
      TEUCHOS_TEST_FOR_EXCEPTION
        (R_bot.empty(), std::logic_error, "R_bot is empty!");
      TEUCHOS_TEST_FOR_EXCEPTION
        (work_.size() == 0, std::logic_error,
         "Workspace array work_ has length zero.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (work_.size() < size_t (R_top.extent(1)), std::logic_error,
         "Workspace array work_ has length = " << work_.size()
         << " < R_top.extent(1) = " << R_top.extent(1) << ".");

      std::vector<Scalar> tau (R_top.extent (1));

      // Our convention for such helper methods is for the immediate
      // parent to allocate workspace (the work_ array in this case).
      //
      // The statement below only works if R_top and R_bot have a
      // nonzero (and the same) number of columns, but we have already
      // checked that above.
      combine_.factor_pair (R_top, R_bot, tau.data(), work_.data());
      return tau;
    }

    void
    factorSecondPass (std::vector<mat_view_type >& topBlocks,
                      std::vector<std::vector<Scalar> >& tauArrays,
                      const int numPartitions) const
    {
      const char prefix[] = "KokkosNodeTsqr::factorSecondPass: ";
      const char suffix[] = "  Please report this bug to the Tpetra developers.";

      if (numPartitions <= 1)
        return; // Done!
      TEUCHOS_TEST_FOR_EXCEPTION
        (topBlocks.size () < size_t (numPartitions), std::logic_error,
         prefix << "topBlocks.size() (= " << topBlocks.size() << ") "
         "< numPartitions (= " << numPartitions << ")." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION
        (tauArrays.size () < size_t (numPartitions-1),
         std::logic_error, prefix << "topBlocks.size() (= "
         << topBlocks.size() << ") < numPartitions-1 (= "
         << (numPartitions-1) << ")." << suffix);
      // The top partition (partition index zero) should always be
      // nonempty if we get this far, so its top block should also be
      // nonempty.
      TEUCHOS_TEST_FOR_EXCEPTION
        (topBlocks[0].empty(), std::logic_error,
         prefix << "topBlocks[0] is empty." << suffix);
      // However, other partitions besides the top one might be empty,
      // in which case their top blocks will be empty.  We skip over
      // the empty partitions in the loop below.
      work_.resize (size_t (topBlocks[0].extent(1)));
      for (int partIdx = 1; partIdx < numPartitions; ++partIdx) {
        if (! topBlocks[partIdx].empty ()) {
          tauArrays[partIdx-1] = factorPair (topBlocks[0], topBlocks[partIdx]);
        }
      }
    }

    void
    applyPair (const ApplyType& applyType,
               const mat_view_type& R_bot,
               const std::vector<Scalar>& tau,
               const mat_view_type& C_top,
               const mat_view_type& C_bot) const
    {
      // Our convention for such helper methods is for the immediate
      // parent to allocate workspace (the work_ array in this case).
      //
      // The statement below only works if C_top, R_bot, and C_bot
      // have a nonzero (and the same) number of columns, but we have
      // already checked that above.
      combine_.apply_pair (applyType, C_top.extent(1), R_bot.extent(1),
                           R_bot.data(), R_bot.stride(1), tau.data(),
                           C_top.data(), C_top.stride(1),
                           C_bot.data(), C_bot.stride(1), work_.data());
    }

    void
    applySecondPass (const ApplyType& applyType,
                     const FactorOutput& factorOutput,
                     std::vector<mat_view_type >& topBlocksOfC,
                     const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy,
                     const bool explicitQ) const
    {
      const char prefix[] = "KokkosNodeTsqr::applySecondPass: ";
      const char suffix[] = "  Please report this bug to the Tpetra developers.";

      const int numParts = factorOutput.numPartitions();
      if (numParts <= 1)
        return; // Done!
      TEUCHOS_TEST_FOR_EXCEPTION
        (topBlocksOfC.size () != size_t (numParts), std::logic_error,
         prefix << "topBlocksOfC.size() (= " << topBlocksOfC.size()
         << ") != number of partitions (= " << numParts << ")."
         << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION
        (factorOutput.secondPassTauArrays.size () != size_t (numParts-1),
         std::logic_error, prefix <<
         "factorOutput.secondPassTauArrays.size() (= "
         << factorOutput.secondPassTauArrays.size()
         << ") != number of partitions minus 1 (= "
         << (numParts-1) << ")." << suffix);
      const LocalOrdinal numCols = topBlocksOfC[0].extent(1);
      work_.resize (size_t (numCols));

      // Top blocks of C are the whole cache blocks.  We only want to
      // affect the top ncols x ncols part of each of those blocks in
      // this method.
      mat_view_type C_top_square (numCols, numCols, topBlocksOfC[0].data(),
                                  topBlocksOfC[0].stride(1));
      if (applyType.transposed ()) {
        // Don't include the topmost (index 0) partition in the
        // iteration; that corresponds to C_top_square.
        for (int partIdx = 1; partIdx < numParts; ++partIdx) {
          // It's legitimate for some partitions not to have any
          // cache blocks.  In that case, their top block will be
          // empty, and we can skip over them.
          const mat_view_type& C_cur = topBlocksOfC[partIdx];
          if (! C_cur.empty()) {
            mat_view_type C_cur_square (numCols, numCols, C_cur.data (),
                                        C_cur.stride (1));
            // If explicitQ: We've already done the first pass and
            // filled the top blocks of C.
            applyPair (applyType, factorOutput.topBlocks[partIdx],
                       factorOutput.secondPassTauArrays[partIdx-1],
                       C_top_square, C_cur_square);
          }
        }
      } else {
        // In non-transposed mode, when computing the first
        // C.extent(1) columns of the explicit Q factor, intranode
        // TSQR would run after internode TSQR (i.e., DistTsqr)
        // (even if only running on a single node in non-MPI mode).
        // Therefore, internode TSQR is responsible for filling the
        // top block of this node's part of the C matrix.
        //
        // Don't include the topmost partition in the iteration;
        // that corresponds to C_top_square.
        for (int partIdx = numParts - 1; partIdx > 0; --partIdx) {
          // It's legitimate for some partitions not to have any
          // cache blocks.  In that case, their top block will be
          // empty, and we can skip over them.
          const mat_view_type& C_cur = topBlocksOfC[partIdx];
          if (! C_cur.empty()) {
            mat_view_type C_cur_square (numCols, numCols,
                                        C_cur.data (),
                                        C_cur.stride (1));
            // The "first" pass (actually the last, only named
            // "first" by analogy with factorFirstPass()) will
            // fill the rest of these top blocks.  For now, we
            // just fill the top n x n part of the top blocks
            // with zeros.
            if (explicitQ) {
              deep_copy (C_cur_square, Scalar {});
            }
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
    const_mat_view_type
    const_top_block (const const_mat_view_type& C,
                     const bool contiguous_cache_blocks) const
    {
      typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;
      blocker_type blocker (C.extent(0), C.extent(1), strategy_);

      // C_top_block is a view of the topmost cache block of C.
      // C_top_block should have >= ncols rows, otherwise either cache
      // blocking is broken or the input matrix C itself had fewer
      // rows than columns.
      const_mat_view_type C_top = blocker.top_block (C, contiguous_cache_blocks);
      return C_top;
    }
  };
} // namespace TSQR

#endif // __TSQR_KokkosNodeTsqr_hpp
