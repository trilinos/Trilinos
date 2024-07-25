// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CRSPADDING_HPP
#define TPETRA_DETAILS_CRSPADDING_HPP

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Util.hpp"
#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace Tpetra {
  namespace Details {

    /// \brief Keep track of how much more space a CrsGraph or
    ///   CrsMatrix needs, when the graph or matrix is the target of a
    ///   doExport or doImport.
    template<class LocalOrdinal, class GlobalOrdinal>
    class CrsPadding {
    private:
      using LO = LocalOrdinal;
      using GO = GlobalOrdinal;

      enum class Phase {
        SAME,
        PERMUTE,
        IMPORT
      };

    public:
      CrsPadding(const int myRank,
                 const size_t /* numSameIDs */,
                 const size_t /* numPermutes */)
        : myRank_(myRank)
      {}

      CrsPadding(const int myRank,
                 const size_t /* numImports */)
        : myRank_(myRank)
      {}

      void
      update_same(
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        const LO whichSame = targetLocalIndex;
        update_impl(Phase::SAME, whichSame, targetLocalIndex,
                    tgtGblColInds, origNumTgtEnt, tgtIsUnique,
                    srcGblColInds, origNumSrcEnt, srcIsUnique);
      }

      void
      update_permute(
        const LO whichPermute, // index in permuteFrom/To
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        update_impl(Phase::PERMUTE, whichPermute, targetLocalIndex,
                    tgtGblColInds, origNumTgtEnt, tgtIsUnique,
                    srcGblColInds, origNumSrcEnt, srcIsUnique);
      }

      void
      update_import(
        const LO whichImport,
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        update_impl(Phase::IMPORT, whichImport, targetLocalIndex,
                    tgtGblColInds, origNumTgtEnt, tgtIsUnique,
                    srcGblColInds, origNumSrcEnt, srcIsUnique);
      }

      void print(std::ostream& out) const {
        const size_t maxNumToPrint =
          Details::Behavior::verbosePrintCountThreshold();
        const size_t size = entries_.size();
        out << "entries: [";
        size_t k = 0;
        for (const auto& keyval : entries_) {
          if (k > maxNumToPrint) {
            out << "...";
            break;
          }
          out << "(" << keyval.first << ", ";
          Details::verbosePrintArray(out, keyval.second,
            "Global column indices", maxNumToPrint);
          out << ")";
          if (k + size_t(1) < size) {
            out << ", ";
          }
          ++k;
        }
        out << "]";
      }

      struct Result {
        size_t numInSrcNotInTgt;
        bool found;
      };

      /// \brief For a given target matrix local row index, return the
      ///   number of unique source column indices to merge into that
      ///   row encountered thus far that are not already in the row,
      ///   and whether we've seen that row already.
      ///
      /// This method relies only on const methods of std::map, and
      /// thus should be thread safe (on host).
      Result
      get_result(const LO targetLocalIndex) const
      {
        auto it = entries_.find(targetLocalIndex);
        if (it == entries_.end()) {
          return {0, false};
        }
        else {
          return {it->second.size(), true};
        }
      }

    private:
      void
      update_impl(
        const Phase phase,
        const LO whichImport,
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        using std::endl;
        std::unique_ptr<std::string> prefix;
        if (verbose_) {
          prefix = createPrefix("update_impl");
          std::ostringstream os;
          os << *prefix << "Start: "
             << "targetLocalIndex=" << targetLocalIndex
             << ", origNumTgtEnt=" << origNumTgtEnt
             << ", origNumSrcEnt=" << origNumSrcEnt << endl;
          std::cerr << os.str();
        }

        // FIXME (08 Feb 2020) We only need to sort and unique
        // tgtGblColInds if we haven't already seen it before.
        size_t newNumTgtEnt = origNumTgtEnt;
        auto tgtEnd = tgtGblColInds + origNumTgtEnt;
        std::sort(tgtGblColInds, tgtEnd);
        if (! tgtIsUnique) {
          tgtEnd = std::unique(tgtGblColInds, tgtEnd);
          newNumTgtEnt = size_t(tgtEnd - tgtGblColInds);
          TEUCHOS_ASSERT( newNumTgtEnt <= origNumTgtEnt );
        }

        if (verbose_) {
          std::ostringstream os;
          os << *prefix << "finished src; process tgt" << endl;
          std::cerr << os.str();
        }

        size_t newNumSrcEnt = origNumSrcEnt;
        auto srcEnd = srcGblColInds + origNumSrcEnt;
        std::sort(srcGblColInds, srcEnd);
        if (! srcIsUnique) {
          srcEnd = std::unique(srcGblColInds, srcEnd);
          newNumSrcEnt = size_t(srcEnd - srcGblColInds);
          TEUCHOS_ASSERT( newNumSrcEnt <= origNumSrcEnt );
        }

        merge_with_current_state(phase, whichImport, targetLocalIndex,
                                 tgtGblColInds, newNumTgtEnt,
                                 srcGblColInds, newNumSrcEnt);
        if (verbose_) {
          std::ostringstream os;
          os << *prefix << "Done" << endl;
          std::cerr << os.str();
        }
      }

      std::vector<GO>&
      get_difference_col_inds(const Phase /* phase */,
                              const LO /* whichIndex */,
                              const LO tgtLclRowInd)
      {
        return entries_[tgtLclRowInd];
      }

      void
      merge_with_current_state(
        const Phase phase,
        const LO whichIndex,
        const LO tgtLclRowInd,
        const GO tgtColInds[], // sorted & merged
        const size_t numTgtEnt,
        const GO srcColInds[], // sorted & merged
        const size_t numSrcEnt)
      {
        using std::endl;
        std::unique_ptr<std::string> prefix;
        if (verbose_) {
          prefix = createPrefix("merge_with_current_state");
          std::ostringstream os;
          os << *prefix << "Start: "
             << "tgtLclRowInd=" << tgtLclRowInd
             << ", numTgtEnt=" << numTgtEnt
             << ", numSrcEnt=" << numSrcEnt << endl;
          std::cerr << os.str();
        }
        // We only need to accumulate those source indices that are
        // not already target indices.  This is because we always have
        // the target indices on input to this function, so there's no
        // need to store them here again.  That still could be a lot
        // to store, but it's better than duplicating target matrix
        // storage.
        //
        // This means that consumers of this data structure need to
        // treat entries_[tgtLclRowInd].size() as an increment, not as
        // the required new allocation size itself.
        //
        // We store
        //
        // difference(union(incoming source indices,
        //                  already stored source indices),
        //            target indices)

        auto tgtEnd = tgtColInds + numTgtEnt;
        auto srcEnd = srcColInds + numSrcEnt;

        // At least one input source index isn't in the target.
        std::vector<GO>& diffColInds =
          get_difference_col_inds(phase, whichIndex, tgtLclRowInd);
        const size_t oldDiffNumEnt = diffColInds.size();

        if (oldDiffNumEnt == 0) {
          if (verbose_) {
            std::ostringstream os;
            os << *prefix << "oldDiffNumEnt=0; call "
              "set_difference(src,tgt,diff)" << endl;
            std::cerr << os.str();
          }
          diffColInds.resize(numSrcEnt);
          auto diffEnd = std::set_difference(srcColInds, srcEnd,
                                             tgtColInds, tgtEnd,
                                             diffColInds.begin());
          const size_t newLen(diffEnd - diffColInds.begin());
          TEUCHOS_ASSERT( newLen <= numSrcEnt );
          diffColInds.resize(newLen);
        }
        else {
          // scratch = union(srcColInds, diffColInds);
          // diffColInds = difference(scratch, tgtColInds);

          const size_t maxUnionSize = numSrcEnt + oldDiffNumEnt;
          if (verbose_) {
            std::ostringstream os;
            os << *prefix << "oldDiffNumEnt=" << oldDiffNumEnt
               << ", maxUnionSize=" << maxUnionSize
               << "; call set_union(src,diff,union)" << endl;
            std::cerr << os.str();
          }
          if (scratchColInds_.size() < maxUnionSize) {
            scratchColInds_.resize(maxUnionSize);
          }
          auto unionBeg = scratchColInds_.begin();
          auto unionEnd = std::set_union(srcColInds, srcEnd,
             diffColInds.begin(), diffColInds.end(),
             unionBeg);
          const size_t unionSize(unionEnd - unionBeg);
          TEUCHOS_ASSERT( unionSize <= maxUnionSize );

          if (verbose_) {
            std::ostringstream os;
            os << *prefix << "oldDiffNumEnt=" << oldDiffNumEnt
               << ", unionSize=" << unionSize << "; call "
              "set_difference(union,tgt,diff)" << endl;
            std::cerr << os.str();
          }
          diffColInds.resize(unionSize);
          auto diffEnd = std::set_difference(unionBeg, unionEnd,
                                             tgtColInds, tgtEnd,
                                             diffColInds.begin());
          const size_t diffLen(diffEnd - diffColInds.begin());
          TEUCHOS_ASSERT( diffLen <= unionSize );
          diffColInds.resize(diffLen);
        }

        if (verbose_) {
          std::ostringstream os;
          os << *prefix << "Done" << endl;
          std::cerr << os.str();
        }
      }

      std::unique_ptr<std::string>
      createPrefix(const char funcName[])
      {
        std::ostringstream os;
        os << "Proc " << myRank_ << ": CrsPadding::" << funcName
           << ": ";
        return std::unique_ptr<std::string>(new std::string(os.str()));
      }

      // imports may overlap with sames and/or permutes, so it makes
      // sense to store them all in one map.
      std::map<LO, std::vector<GO> > entries_;
      std::vector<GO> scratchColInds_;
      int myRank_ = -1;
      bool verbose_ = Behavior::verbose("CrsPadding");
    };
  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CRSPADDING_HPP
