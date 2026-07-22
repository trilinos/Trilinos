#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/parallel/CouplingVersions.hpp"
#include <cassert>
#include "Parallel.hpp"

namespace stk {

MPITagManager::MPITagManager(int deletionGroupSize, int delayCount) :
  m_keyManager(std::make_shared<impl::MPIKeyManager>()),
  m_commData(impl::CommCompare(m_keyManager)),

  m_deletionGroupSize(deletionGroupSize),
  m_delayCount(delayCount)
{
  int isInitialized;
  MPI_Initialized(&isInitialized);
  STK_ThrowRequireMsg(isInitialized, "MPI must be initialized prior to constructing MPITagManager");

  int flag;
  int* val;
  MPI_Comm_get_attr(parallel_machine_world(), MPI_TAG_UB, &val, &flag);
  STK_ThrowRequireMsg(flag, "This MPI implementation is erroneous");
  STK_ThrowRequireMsg(*val >= m_tagMax, "MPI_TAG_UB must be at least " + std::to_string(m_tagMax));
  m_tagMax = *val - 1;

  m_callbackUID = m_keyManager->get_UID();
}

MPITagManager::~MPITagManager()
{
  m_keyManager->execute_callbacks_immediately(m_callbackUID);
}


MPITag MPITagManager::get_tag(MPI_Comm userComm, int tagHint)
{
  MPI_Comm comm;
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
  comm = m_commReplacer.get_copied_comm(userComm);
#else
  comm = userComm;
#endif

  if (!m_keyManager->has_key(comm))
  {
    m_keyManager->register_callback(userComm, m_callbackUID, std::bind(&MPITagManager::erase_comm, this, std::placeholders::_1));
    m_commData.emplace(std::piecewise_construct, std::make_tuple(comm), std::make_tuple(comm, m_tagMin, m_deletionGroupSize, m_delayCount, m_tagMax+1));
  }

  auto& commData = m_commData.at(comm);
  auto newTag = get_new_tag(commData, tagHint);

#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
  auto tagData = std::make_shared<impl::MPITagData>(this, userComm, comm, newTag);
#else
  auto tagData = std::make_shared<impl::MPITagData>(this, comm, newTag);
#endif
  commData.insert(tagData);

  return MPITag(tagData);
}


int MPITagManager::get_new_tag(impl::CommTagInUseList& commData, int tagHint)
{
  int newTag = -1;
  if (tagHint == MPI_ANY_TAG) {
    newTag = get_any_tag(commData);
  } else {
    newTag = new_tag_search(commData, tagHint);
  }

  assert(newTag >= 0);
  STK_ThrowRequireMsg(newTag <= m_tagMax, "New tag must be <= " + std::to_string(m_tagMax));
  check_same_value_on_all_procs_debug_only(commData.get_comm(), newTag);

  return newTag;
}


int MPITagManager::get_any_tag(impl::CommTagInUseList& commData)
{
  int newTag = commData.get_min_free_tag();
  STK_ThrowRequireMsg(newTag <= m_tagMax, std::string("MPI tag supply exhausted: there can only be ") + 
                                      std::to_string(m_tagMax) + " tags in use at any time");

  return newTag;
}


int MPITagManager::new_tag_search(impl::CommTagInUseList& commData, int tagHint)
{
  const auto& tags = commData.get_tags();

  int newTag = -1;
  auto it = tags.find(tagHint);
  if (it == tags.end()) {
    newTag = tagHint;
  } else {
    // find next available tag
    int prevVal = it->first;
    while (it != tags.end())
    {
      if ( (it->first - prevVal) > 1)
      {
        newTag = prevVal + 1;
        break;
      }
      prevVal = it->first;
      it++;
    }

    if (newTag == -1) { // if no spaces between existing tags found
      newTag = (std::prev(tags.end()))->first + 1;
    }
  }

  return newTag;
}


void MPITagManager::free_tag_local(impl::MPITagData& tag)
{

  MPI_Comm comm;
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
  comm = tag.get_comm_internal();
#else
  comm = tag.get_comm();
#endif
  m_commData.at(comm).erase(tag);
}


void MPITagManager::erase_comm(MPI_Comm origComm)
{
  MPI_Comm comm;
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
  comm = m_commReplacer.get_copied_comm(origComm);
#else
  comm = origComm;
#endif

  STK_ThrowRequireMsg(m_commData.count(comm) == 1, "Cannot free MPI Comm that is not assigned (possible double free)");
  m_commData.erase(comm);

#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
  m_commReplacer.delete_comm_pair(origComm);
#endif
}


void MPITagManager::check_same_value_on_all_procs_debug_only([[maybe_unused]] MPI_Comm comm, [[maybe_unused]] int newVal)
{
#ifndef NDEBUG
  int commSize, myrank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &myrank);

  const int root = 0;
  int recvSize = myrank == root ? commSize : 0;
  std::vector<int> recvBuf(recvSize);
  MPI_Gather(&newVal, 1, MPI_INT, recvBuf.data(), 1, MPI_INT, root, comm);

  if (myrank == root) {
    bool isSame = true;
    for (auto v : recvBuf) {
      isSame = isSame && v == recvBuf[0];
    }

    STK_ThrowRequireMsg(isSame, "Calls to MPICommManager must be collective");
  }

  MPI_Barrier(comm);
#endif
}


MPITagManager& get_mpi_tag_manager()
{
  static int delayCount = -1;
  if (delayCount < 0)
  {
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    // some MPI_Barrier algorithms use a tree based algorithm,
    // travering it down and them up again, which likely takes
    // 2 * log2(number of ranks)
    delayCount = std::max(2*std::ceil(std::log2(commSize)), 4.0);
  }

  int deletionGroupSize = std::max(33, 2 * delayCount);

  static MPITagManager tagManager(deletionGroupSize, delayCount);
  return tagManager;
}

}
