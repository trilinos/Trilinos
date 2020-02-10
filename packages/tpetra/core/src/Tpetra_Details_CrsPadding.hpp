#ifndef TPETRA_DETAILS_CRSPADDING_HPP
#define TPETRA_DETAILS_CRSPADDING_HPP

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Util.hpp"
#include <algorithm>
#include <type_traits>
#include <vector>

namespace Tpetra {
  namespace Details {

    template<class LocalOrdinal, class GlobalOrdinal, class DeviceType>
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
      struct create_from_sames_and_permutes_tag {};
      static constexpr create_from_sames_and_permutes_tag
        create_from_sames_and_permutes {};
      CrsPadding(create_from_sames_and_permutes_tag /* tag */,
                 const size_t /* numSameIDs */,
                 const size_t /* numPermutes */)
      {}

      struct create_from_imports_tag {};
      static constexpr create_from_imports_tag create_from_imports {};
      CrsPadding(create_from_imports_tag /* tag */,
                 const size_t /* numImports */)
      {}

      void
      update_same(
        size_t& tgtNumDups, // accumulator
        size_t& srcNumDups, // accumulator
        size_t& unionNumDups, // accumulator
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        const LO whichSame = targetLocalIndex;
        update_impl(Phase::SAME,
                    tgtNumDups, srcNumDups, unionNumDups,
                    whichSame, targetLocalIndex,
                    tgtGblColInds, origNumTgtEnt, tgtIsUnique,
                    srcGblColInds, origNumSrcEnt, srcIsUnique);
      }

      void
      update_permute(
        size_t& tgtNumDups, // accumulator
        size_t& srcNumDups, // accumulator
        size_t& unionNumDups, // accumulator
        const LO whichPermute, // index in permuteFrom/To
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        update_impl(Phase::PERMUTE,
                    tgtNumDups, srcNumDups, unionNumDups,
                    whichPermute, targetLocalIndex,
                    tgtGblColInds, origNumTgtEnt, tgtIsUnique,
                    srcGblColInds, origNumSrcEnt, srcIsUnique);
      }

      void
      update_import(
        size_t& tgtNumDups, // accumulator
        size_t& srcNumDups, // accumulator
        size_t& unionNumDups, // accumulator
        const LO whichImport,
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        update_impl(Phase::IMPORT,
                    tgtNumDups, srcNumDups, unionNumDups,
                    whichImport, targetLocalIndex,
                    tgtGblColInds, origNumTgtEnt, tgtIsUnique,
                    srcGblColInds, origNumSrcEnt, srcIsUnique);
      }

      void print(std::ostream& out) const {
        out << "increase: " << increase() << ", ";
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

      /// \brief Increase (increment in the number of entries) in
      ///   required allocation size.
      size_t increase() const {
        return increase_;
      }

      struct Result {
        size_t allocSize;
        bool found;
      };

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
        size_t& tgtNumDups,
        size_t& srcNumDups,
        size_t& unionNumDups,
        const LO whichImport,
        const LO targetLocalIndex,
        GO tgtGblColInds[],
        const size_t origNumTgtEnt,
        const bool tgtIsUnique,
        GO srcGblColInds[],
        const size_t origNumSrcEnt,
        const bool srcIsUnique)
      {
        // FIXME (08 Feb 2020) We only need to sort and unique
        // tgtGblColInds if we haven't already seen it before.
        size_t newNumTgtEnt = origNumTgtEnt;
        auto tgtEnd = tgtGblColInds + origNumTgtEnt;
        std::sort(tgtGblColInds, tgtEnd);
        if (! tgtIsUnique) {
          tgtEnd = std::unique(tgtGblColInds, tgtEnd);
          newNumTgtEnt = size_t(tgtEnd - tgtGblColInds);
        }
        tgtNumDups += (origNumTgtEnt - newNumTgtEnt);

        size_t newNumSrcEnt = origNumSrcEnt;
        auto srcEnd = srcGblColInds + origNumSrcEnt;
        std::sort(srcGblColInds, srcEnd);
        if (! srcIsUnique) {
          srcEnd = std::unique(srcGblColInds, srcEnd);
          newNumSrcEnt = size_t(srcEnd - srcGblColInds);
        }
        srcNumDups += (origNumSrcEnt - newNumSrcEnt);

        size_t unionNumEnt = 0;
        merge_with_current_state(phase, unionNumEnt,
                                 whichImport, targetLocalIndex,
                                 tgtGblColInds, newNumTgtEnt,
                                 srcGblColInds, newNumSrcEnt);
        unionNumDups += (newNumTgtEnt + newNumSrcEnt - unionNumEnt);
        if (unionNumEnt > origNumTgtEnt) {
          increase_ += (unionNumEnt - origNumTgtEnt);
        }
      }

      std::vector<GO>&
      get_union_col_inds(const Phase /* phase */,
                         const LO /* whichIndex */,
                         const LO tgtLclRowInd)
      {
        return entries_[tgtLclRowInd];
      }

      void
      merge_with_current_state(
        const Phase phase,
        size_t& unionNumEnt,
        const LO whichIndex,
        const LO tgtLclRowInd,
        const GO tgtColInds[], // sorted & merged
        const size_t numTgtEnt,
        const GO srcColInds[], // sorted & merged
        const size_t numSrcEnt)
      {
        std::vector<GO>& unionColInds =
          get_union_col_inds(phase, whichIndex, tgtLclRowInd);

        if (unionColInds.size() == 0) {
          auto tgtEnd = tgtColInds + numTgtEnt;
          auto srcEnd = srcColInds + numSrcEnt;
          const size_t numInCommon = Details::countNumInCommon(
            srcColInds, srcEnd, tgtColInds, tgtEnd);
          unionNumEnt = numTgtEnt + numSrcEnt - numInCommon;
          unionColInds.resize(unionNumEnt);
          (void) std::set_union(tgtColInds, tgtEnd,
                                srcColInds, srcEnd,
                                unionColInds.begin());
        }
        else {
          // We've already seen the target graph/matrix row before, so
          // we need not even look at tgtColInds.
          const size_t oldUnionSize = unionColInds.size();

          const size_t maxUnionSize = numSrcEnt + unionColInds.size();
          if (scratchColInds_.size() < maxUnionSize) {
            scratchColInds_.resize(maxUnionSize);
          }
          auto scratchEnd = std::set_union(
            srcColInds, srcColInds + numSrcEnt,
            unionColInds.begin(), unionColInds.end(),
            scratchColInds_.begin());
          unionNumEnt = size_t(scratchEnd - scratchColInds_.begin());
          unionColInds.resize(unionNumEnt);
          std::copy(scratchColInds_.begin(), scratchColInds_.end(),
                    unionColInds.begin());
        }
      }

      // imports may overlap with sames and/or permutes, so it makes
      // sense to store them all in one map.
      std::map<LO, std::vector<GO> > entries_;
      std::vector<GO> scratchColInds_;
      size_t increase_ = 0;
    };
  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CRSPADDING_HPP
