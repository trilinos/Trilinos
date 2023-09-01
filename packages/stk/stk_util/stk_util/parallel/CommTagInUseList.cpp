#include "stk_util/parallel/CommTagInUseList.hpp"
#include <limits>
#include <stdexcept>
#include <string>

namespace stk {
namespace impl {

CommTagInUseList::~CommTagInUseList()
{
  for (auto& deletionGroup : m_deletionGroups) {
    if (deletionGroup.is_barrier_in_progress()) {
      deletionGroup.finish_barrier();
      erase_internal(deletionGroup.get_tags());
    }
  }

  for (auto& weak_tag_ptr : m_tags) {
    if (auto tag_shared_ptr = weak_tag_ptr.second.lock()) {
      tag_shared_ptr->set_free();
    }
  }
}


void CommTagInUseList::insert(std::shared_ptr<MPITagData> newTag)
{
  #ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    assert(are_comms_identical(m_comm, newTag->get_comm_internal()));
  #else
    assert(are_comms_identical(m_comm, newTag->get_comm()));
  #endif

  int tagVal = newTag->get_tag();
  STK_ThrowRequireMsg(m_tags.count(tagVal) == 0, "Cannot create new tag with same value as existing tag");

  auto iteratorBoolPair = m_tags.insert(std::make_pair(tagVal, newTag));
  assert(iteratorBoolPair.second);

  if (newTag->get_tag() == m_minFreeTag)
  {
    auto it = iteratorBoolPair.first;
    int prevVal = it->first;
    bool found = false;
    while (it != m_tags.end())
    {
      int currentVal = it->first;
      if (currentVal - prevVal > 1)
      {
        m_minFreeTag = prevVal + 1;
        found = true;
        break;
      }
      prevVal = it->first;
      it++;
    }

    if (!found)
      m_minFreeTag = (std::prev(m_tags.end()))->first + 1;
  }

  increment_entry_count();
}


void CommTagInUseList::erase(MPITagData& tag)
{
  #ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    assert(are_comms_identical(m_comm, tag.get_comm_internal()));
  #else
    assert(are_comms_identical(m_comm, tag.get_comm()));
  #endif

  auto& deletionGroup  = m_deletionGroups.back();
  deletionGroup.insert_tag(tag.get_tag());

  if (deletionGroup.get_tags().size() == m_deletionGroupSize) {
    deletionGroup.start_barrier(m_entryCount);

    MPI_Comm comm;
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    comm = tag.get_comm_internal();
#else
    comm = tag.get_comm();
#endif
    m_deletionGroups.emplace_back(comm, m_barrierTag);
  }

  tag.set_free();
}


int CommTagInUseList::get_min_free_tag()
{ 
  check_tag_deletion_completion();
  return m_minFreeTag;
}


const std::map<int, std::weak_ptr<MPITagData>>& CommTagInUseList::get_tags()
{ 
  check_tag_deletion_completion();
  return m_tags;
}


void CommTagInUseList::check_tag_deletion_completion()
{
  int groupsFreed = 0;
  for ( auto& deletionGroup : m_deletionGroups)
  {
    bool barrierInProgress = deletionGroup.is_barrier_in_progress();
    if (barrierInProgress) {
      bool isDelaySatisfied = is_delay_satisfied(deletionGroup.get_entry_count_barrier_start());

      if (isDelaySatisfied) {
        deletionGroup.finish_barrier();
        erase_internal(deletionGroup.get_tags());
        groupsFreed++; 
      } else {
        deletionGroup.test_barrier();
      }
    } 
    
    if (!barrierInProgress){
      break;
    }
  }

  if (groupsFreed > 0)
  {
    auto begin = m_deletionGroups.begin();
    auto end   = m_deletionGroups.begin();
    std::advance(end, groupsFreed);
    m_deletionGroups.erase(begin, end);
  }
}


CommTagInUseList::EntryCountInt CommTagInUseList::increment_entry_count()
{
  if (std::numeric_limits<EntryCountInt>::max() == m_entryCount)
    throw std::runtime_error("Cannot call get_tag() more than " +
            std::to_string(std::numeric_limits<EntryCountInt>::max()) + " times");

  return ++m_entryCount;
}


bool CommTagInUseList::is_delay_satisfied(EntryCountInt startCount)
{
  auto currentCount = m_entryCount;
  assert(currentCount >= startCount);
  return currentCount - startCount >= m_delayCount;
}


void CommTagInUseList::erase_internal(const std::vector<int>& tagVals)
{
  for (auto& tag_val : tagVals) {
    STK_ThrowRequireMsg(m_tags.count(tag_val) == 1, "Cannot free tag this has not been assigned (possible double free)");
    m_tags.erase(tag_val);
    m_minFreeTag = std::min(m_minFreeTag, tag_val);
  }
}


bool CommTagInUseList::are_comms_identical(MPI_Comm comm1, MPI_Comm comm2)
{
  int result;
  MPI_Comm_compare(comm1, comm2, &result);
  return result == MPI_IDENT;
}


}
}