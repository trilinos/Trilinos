#ifndef NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETFINDER_HPP_
#define NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETFINDER_HPP_
#include <stk_emend/independent_set/IndependentSetInfoToInfoGraph.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <vector>
#include <set>

namespace independent_set
{

template <class Info>
class IndependentSetFinder
{
public:
    static std::vector<Info> find_independent_set(const std::vector<Info> &infos,
                                                  const typename Info::Comparator &comparator,
                                                  const typename Info::ConflictFinder &conflictFinder,
                                                  const bool includeSetsWithHigherPriorityOutNeighbors,
                                                  stk::ParallelMachine comm);

private:
    enum InOutUnknownStatus
    {
        UNKNOWN_STATUS = 1,
        IN = 2,
        OUT = 3
    };

    const typename Info::Comparator &mInfoComparator;
    const bool mIncludeSetsWithHigherPriorityOutNeighbors;
    const std::vector<Info> &mInfos;
    const IndependentSetInfoToInfoGraph<Info> mConflictingInfoGraph;
    stk::ParallelMachine mComm;
    std::vector<size_t> mSortedIndicesForGlobalIdsOfOffProcInfos;

    IndependentSetFinder(const std::vector<Info> &info, const typename Info::Comparator &comparator, const typename Info::ConflictFinder &conflictFinder, const bool includeSetsWithHigherPriorityOutNeighbors, stk::ParallelMachine comm)
    : mInfoComparator{comparator},
      mIncludeSetsWithHigherPriorityOutNeighbors(includeSetsWithHigherPriorityOutNeighbors),
      mInfos{info},
      mConflictingInfoGraph{mInfos,conflictFinder},
      mComm{comm}
      {
          populate_sorted_info_indices();
      }

    std::vector<Info> find_independent_set() const;
    void assign_info_to_be_in_and_conflicting_neighbors_as_out(size_t iInfo, std::vector<InOutUnknownStatus> &inOutStatus) const;
    std::vector<InOutUnknownStatus> determine_independent_set_in_out_status_for_infos() const;
    std::vector<Info> cull_down_to_independent_infos(const std::vector<InOutUnknownStatus>& inOutStatus) const;
    bool determine_independent_set_in_out_status_for_infos_locally(std::vector<InOutUnknownStatus>& inOutStatus) const;
    void parallel_communicate_status_of_IN_infos(std::vector<InOutUnknownStatus> &inOutStatus) const;
    void pack_status_of_shared_infos(const std::vector<InOutUnknownStatus> &inOutStatus, stk::CommSparse &commSparse)const;
    void receive_status_of_shared_infos(std::vector<InOutUnknownStatus> &inOutStatus, stk::CommSparse &commSparse) const;
    void populate_sorted_info_indices();
    std::function<bool(size_t,size_t)> build_function_to_check_is_neighbor_higher_priority_than_info(const std::vector<InOutUnknownStatus>& inOutStatus) const;
};

template <class Info>
std::vector<Info> IndependentSetFinder<Info>::find_independent_set(const std::vector<Info> &infos,
                                                                                   const typename Info::Comparator &comparator,
                                                                                   const typename Info::ConflictFinder &conflictFinder,
                                                                                   const bool includeSetsWithHigherPriorityOutNeighbors,
                                                                                   stk::ParallelMachine comm)
{
    IndependentSetFinder finder(infos, comparator, conflictFinder, includeSetsWithHigherPriorityOutNeighbors, comm);
    return finder.find_independent_set();
}

template <class Info>
std::vector<Info> IndependentSetFinder<Info>::find_independent_set() const
{
    std::vector<InOutUnknownStatus> infoInOutStatus = determine_independent_set_in_out_status_for_infos();
    return cull_down_to_independent_infos(infoInOutStatus);
}

template <class Info>
void IndependentSetFinder<Info>::assign_info_to_be_in_and_conflicting_neighbors_as_out(size_t iInfo, std::vector<IndependentSetFinder<Info>::InOutUnknownStatus> &inOutStatus) const
{
    inOutStatus[iInfo] = IN;
    for(size_t iAdjInfo : mConflictingInfoGraph.get_conflicting_infos_for(iInfo))
    {
        inOutStatus[iAdjInfo] = OUT;
    }
}

template <class Info>
std::function<bool(size_t,size_t)> IndependentSetFinder<Info>::build_function_to_check_is_neighbor_higher_priority_than_info(const std::vector<InOutUnknownStatus>& inOutStatus) const
{
    auto is_neighbor_higher_priority_than_info = [&](size_t iInfoNbr, size_t iInfo)
    {
      ThrowRequireWithSierraHelpMsg(inOutStatus[iInfoNbr] != IN);
      const bool needToConsiderNeighbor = !mIncludeSetsWithHigherPriorityOutNeighbors || inOutStatus[iInfoNbr] != OUT;
      return needToConsiderNeighbor && mInfoComparator.is_first_higher_priority_than_second(mInfos[iInfoNbr], mInfos[iInfo]);
    };

    return is_neighbor_higher_priority_than_info;
}


template <class Info>
bool IndependentSetFinder<Info>::determine_independent_set_in_out_status_for_infos_locally(std::vector<InOutUnknownStatus>& infoInOutStatus) const
{
    const auto is_neighbor_higher_priority_than_info = build_function_to_check_is_neighbor_higher_priority_than_info(infoInOutStatus);
    int thisProc = stk::parallel_machine_rank(mComm);
    bool didAnythingChange{false};
    bool done=false;
    while(!done)
    {
        done = true;
        for(size_t iInfo = 0; iInfo < mInfos.size(); ++iInfo)
        {
            if( thisProc == mInfos[iInfo].get_owner() && infoInOutStatus[iInfo] == UNKNOWN_STATUS)
            {
                const bool isHighestPriority =
                    mConflictingInfoGraph.is_info_higher_priority_than_conflicting_neighbors(iInfo, is_neighbor_higher_priority_than_info);
                if(isHighestPriority)
                {
                    assign_info_to_be_in_and_conflicting_neighbors_as_out(iInfo, infoInOutStatus);
                    didAnythingChange=true;
                    done = !mIncludeSetsWithHigherPriorityOutNeighbors;
                }
                else if (!mIncludeSetsWithHigherPriorityOutNeighbors)
                {
                    infoInOutStatus[iInfo] = OUT;
                }
            }
        }
    }

    return didAnythingChange;
}

template <class Info>
void IndependentSetFinder<Info>::pack_status_of_shared_infos(const std::vector<InOutUnknownStatus> &infoInOutStatus, stk::CommSparse &commSparse) const
{
    stk::pack_and_communicate(commSparse,[&]()
    {
        for(size_t iInfo = 0; iInfo < mInfos.size(); ++iInfo)
        {
            if(infoInOutStatus[iInfo] != UNKNOWN_STATUS)
            {
                const Info& info = mInfos[iInfo];
                for ( const int procId : info.get_procs_that_need_to_know_about_this_info())
                {
                    if ( procId != commSparse.parallel_rank())
                    {
                        commSparse.send_buffer(procId).pack<typename Info::GlobalId>(info.get_unique_id());
                        commSparse.send_buffer(procId).pack<InOutUnknownStatus>(infoInOutStatus[iInfo]);
                    }
                }
            }
        }
    });
}

template <class Info>
size_t get_info_index(const std::vector<Info> &infos, const std::vector<size_t> & sortedIndices, const typename Info::GlobalId & infoId)
{
  auto it = std::lower_bound(sortedIndices.begin(), sortedIndices.end(), infoId,
            [&infos](size_t lhs, const typename Info::GlobalId & rhs) -> bool { return infos[lhs].get_unique_id() < rhs; });
  ThrowAssert(it != sortedIndices.end() && infos[*it].get_unique_id() == infoId);
  return *it;
}

template <class Info>
void IndependentSetFinder<Info>::receive_status_of_shared_infos(std::vector<InOutUnknownStatus> &inOutStatus, stk::CommSparse &commSparse) const
{
    const auto& sortedIndicesForGlobalIdsOfOffProcInfos = mSortedIndicesForGlobalIdsOfOffProcInfos;
    const auto& infos = mInfos;
    stk::unpack_communications(commSparse, [&](int procId)
     {
         typename Info::GlobalId infoId;
         commSparse.recv_buffer(procId).unpack<typename Info::GlobalId>(infoId);
         InOutUnknownStatus oneInOutStatus;
         commSparse.recv_buffer(procId).unpack<InOutUnknownStatus>(oneInOutStatus);

         const size_t iInfo = get_info_index(infos, sortedIndicesForGlobalIdsOfOffProcInfos, infoId);
         ThrowRequireWithSierraHelpMsg(inOutStatus[iInfo] == UNKNOWN_STATUS || inOutStatus[iInfo] == oneInOutStatus);
         inOutStatus[iInfo] = oneInOutStatus;
     });
}

template <class Info>
void IndependentSetFinder<Info>::parallel_communicate_status_of_IN_infos(std::vector<InOutUnknownStatus> &infosInOutStatus) const
{
    stk::CommSparse commSparse(mComm);
    pack_status_of_shared_infos(infosInOutStatus, commSparse);
    receive_status_of_shared_infos(infosInOutStatus, commSparse);
}

template <class Info>
std::vector<typename IndependentSetFinder<Info>::InOutUnknownStatus> IndependentSetFinder<Info>::determine_independent_set_in_out_status_for_infos() const
{
    std::vector<InOutUnknownStatus> infoInOutStatus(mInfos.size(), UNKNOWN_STATUS);
    bool done = false;
    while(!done)
    {
      done = true;
      const bool didAnythingChange=determine_independent_set_in_out_status_for_infos_locally(infoInOutStatus);
      if(stk::is_true_on_any_proc(mComm,didAnythingChange))
      {
          done = !mIncludeSetsWithHigherPriorityOutNeighbors;
          parallel_communicate_status_of_IN_infos(infoInOutStatus);
      }
    }

    return infoInOutStatus;
}

template <class Info>
std::vector<Info> IndependentSetFinder<Info>::cull_down_to_independent_infos(const std::vector<InOutUnknownStatus>& infosInOutStatus) const
{
    std::vector<Info> independentnodesToMoveToTarget;
    for(size_t iNodeToMoveToTarget = 0; iNodeToMoveToTarget < mInfos.size(); ++iNodeToMoveToTarget)
    {
        if(infosInOutStatus[iNodeToMoveToTarget] == IN)
            independentnodesToMoveToTarget.push_back(mInfos[iNodeToMoveToTarget]);
    }
    return independentnodesToMoveToTarget;
}

template <class Info>
void IndependentSetFinder<Info>::populate_sorted_info_indices()
{
    for(size_t iInfo = 0; iInfo < mInfos.size(); ++iInfo)
        if(mInfos[iInfo].get_procs_that_need_to_know_about_this_info().size()>1)
            mSortedIndicesForGlobalIdsOfOffProcInfos.push_back(iInfo);

    const auto & infos = mInfos;
    std::sort(mSortedIndicesForGlobalIdsOfOffProcInfos.begin(), mSortedIndicesForGlobalIdsOfOffProcInfos.end(),
        [&infos](size_t i1, size_t i2) {return infos[i1].get_unique_id() < infos[i2].get_unique_id();});
}

}

#endif /* NGS_INDEPENDENT_SET_LIB_INDEPENDENTSETFINDER_HPP_ */
