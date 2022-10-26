#ifndef stk_util_parallel_CommManagerData
#define stk_util_parallel_CommManagerData

#include <stddef.h>
#include <map>
#include <vector>
#include <memory>
#include <iostream>
#include "stk_util/parallel/MPITag.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, ThrowRequire
#include "stk_util/parallel/DeletionGroup.hpp"

namespace stk {
namespace impl {


class CommTagInUseList
{
  public: 
    CommTagInUseList(MPI_Comm comm, int minTag, int deletionGroupSize, int delayCount, int barrierTag) :
      m_comm(comm),
      m_entryCount(0),
      m_minFreeTag(minTag),
      m_deletionGroupSize(deletionGroupSize),
      m_delayCount(delayCount),
      m_barrierTag(barrierTag)
    {
      m_deletionGroups.emplace_back(comm, barrierTag);
    }

    ~CommTagInUseList();

    CommTagInUseList(const CommTagInUseList&) = delete;

    CommTagInUseList& operator=(const CommTagInUseList&) = delete;

    MPI_Comm get_comm() const { return m_comm; }

    void insert(std::shared_ptr<MPITagData> new_tag);

    void erase(MPITagData& tag);

    int get_min_free_tag();

    const std::map<int, std::weak_ptr<MPITagData>>& get_tags();

  private:
    using EntryCountInt = DeletionGroup::EntryCountInt;

    using MPITagInt = int;

    void check_tag_deletion_completion();

    EntryCountInt increment_entry_count();

    bool is_delay_satisfied(EntryCountInt startCount);

    void erase_internal(const std::vector<int>& tags);

    bool are_comms_identical(MPI_Comm comm1, MPI_Comm comm2);

    MPI_Comm m_comm;
    EntryCountInt m_entryCount;
    std::map<MPITagInt, std::weak_ptr<MPITagData>> m_tags;
    std::deque<DeletionGroup> m_deletionGroups;
    int m_minFreeTag;
    const size_t m_deletionGroupSize;
    const EntryCountInt m_delayCount;
    const int m_barrierTag;
};


}
}

#endif
