#ifndef NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETINFOTOINFOGRAPH_HPP_
#define NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETINFOTOINFOGRAPH_HPP_

#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <functional>

namespace independent_set
{

template <class INFO_TYPE>
class IndependentSetInfoToInfoGraph
{
public:
    IndependentSetInfoToInfoGraph(const std::vector<INFO_TYPE> &infos, const typename INFO_TYPE::ConflictFinder &conflictFinder)
    : mInfos(infos),
      mConflictFinder(conflictFinder)
    {
        build_overlapping_infos(infos, mConflictingIdsToContainingInfos);
    }

    std::set<size_t> get_conflicting_infos_for(size_t iInfo) const
    {
      const auto & info = mInfos[iInfo];
      std::set<size_t> conflictingInfos;

      for(const size_t iInfoNbr : mConflictFinder.get_other_conflicting_infos(info))
          conflictingInfos.insert(iInfoNbr);

      for(const typename INFO_TYPE::ExclusionIdentifierType &conflictingId : info.get_conflicting_ids())
          for(size_t iInfoNbr : mConflictingIdsToContainingInfos.at(conflictingId))
              if(iInfo != iInfoNbr)
                  conflictingInfos.insert(iInfoNbr);

      return conflictingInfos;
    }

    bool is_info_higher_priority_than_conflicting_neighbors(size_t iInfo, const std::function<bool(size_t, size_t)> & is_neighbor_higher_priority_than_info) const
    {
      const auto & info = mInfos[iInfo];

      for(const size_t iInfoNbr : mConflictFinder.get_other_conflicting_infos(info))
         if (is_neighbor_higher_priority_than_info(iInfoNbr, iInfo))
           return false;

      for(const typename INFO_TYPE::ExclusionIdentifierType &conflictingId : info.get_conflicting_ids())
          for(size_t iInfoNbr : mConflictingIdsToContainingInfos.at(conflictingId))
              if(iInfo != iInfoNbr && is_neighbor_higher_priority_than_info(iInfoNbr, iInfo))
                return false;

      return true;
    }

private:
    static void build_overlapping_infos(const std::vector<INFO_TYPE>& infos, std::unordered_map<typename INFO_TYPE::ExclusionIdentifierType, std::vector<size_t> > & conflictingIdsToContainingInfos)
    {
        for(size_t iInfo = 0; iInfo < infos.size(); ++iInfo)
        {
            const INFO_TYPE& oneInfo = infos[iInfo];
            for(const typename INFO_TYPE::ExclusionIdentifierType &conflictingId : oneInfo.get_conflicting_ids())
            {
                static constexpr size_t estimateMaxNumConflicts{20};
                std::vector<size_t> &infosConflictingWithId = conflictingIdsToContainingInfos[conflictingId];
                if(infosConflictingWithId.empty())
                    infosConflictingWithId.reserve(estimateMaxNumConflicts);
                infosConflictingWithId.push_back(iInfo);
            }
        }
    }

    const std::vector<INFO_TYPE> & mInfos;
    const typename INFO_TYPE::ConflictFinder & mConflictFinder;
    std::unordered_map<typename INFO_TYPE::ExclusionIdentifierType, std::vector<size_t> >  mConflictingIdsToContainingInfos;
};

}

#endif /* NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETINFOTOINFOGRAPH_HPP_ */
