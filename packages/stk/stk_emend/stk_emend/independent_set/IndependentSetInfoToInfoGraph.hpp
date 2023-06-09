#ifndef NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETINFOTOINFOGRAPH_HPP_
#define NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETINFOTOINFOGRAPH_HPP_

#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <functional>
#include <algorithm>

namespace independent_set
{

template <class INFO_TYPE>
class IndependentSetInfoToInfoGraph
{
public:
    IndependentSetInfoToInfoGraph(const std::vector<INFO_TYPE> &infos)
    : mInfos(infos)
    {
        build_overlapping_infos(infos, mConflictingIdsToContainingInfos);
    }

    const std::vector<size_t> &get_infos_with_conflicting_id(typename INFO_TYPE::ExclusionIdentifierType id) const
    {
        return mConflictingIdsToContainingInfos.at(id);
    }

    template <typename FUNC_TO_CALL>
    void for_each_conflicting_info(
                            const size_t iInfo,
                            const FUNC_TO_CALL &func) const
    {
        for(const typename INFO_TYPE::ExclusionIdentifierType &conflictingId : mInfos[iInfo].get_conflicting_ids())
            for(size_t iInfoNbr : mConflictingIdsToContainingInfos.at(conflictingId))
                if(iInfo != iInfoNbr)
                    func(iInfoNbr);
    }

    template <typename FUNC_TO_CALL>
    bool is_true_for_any_conflicting_info(
                            const size_t iInfo,
                            const FUNC_TO_CALL &func) const
    {
        for(const typename INFO_TYPE::ExclusionIdentifierType &conflictingId : mInfos[iInfo].get_conflicting_ids())
            for(size_t iInfoNbr : mConflictingIdsToContainingInfos.at(conflictingId))
                if(iInfo != iInfoNbr)
                    if(func(iInfoNbr))
                        return true;
        return false;
    }

    bool is_info_higher_priority_than_all_conflicting_neighbors(
                            size_t iInfo,
                            const typename INFO_TYPE::Comparator &infoComparator) const
    {
        if(is_true_for_any_conflicting_info(
                                iInfo,
                                [&](const size_t iAdjInfo)
                                {
                                    return !infoComparator.is_first_higher_priority_than_second(mInfos[iInfo], mInfos[iAdjInfo]);
                                }))
            return false;
        return true;
    }

    bool is_lower_priority_than_any_neighbor(
                            const size_t iInfo,
                            const typename INFO_TYPE::Comparator &infoComparator) const
    {
        return is_true_for_any_conflicting_info(
                                iInfo,
                                [&](const size_t iAdjInfo)
                                {
                                    return infoComparator.is_first_higher_priority_than_second(mInfos[iAdjInfo], mInfos[iInfo]);
                                });
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
    std::unordered_map<typename INFO_TYPE::ExclusionIdentifierType, std::vector<size_t> >  mConflictingIdsToContainingInfos;
};

}

#endif /* NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETINFOTOINFOGRAPH_HPP_ */
