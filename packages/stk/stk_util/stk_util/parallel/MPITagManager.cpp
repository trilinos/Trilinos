#include "stk_util/parallel/MPITagManager.hpp"
#include <cassert>

namespace stk {

namespace Impl {

MPITagData::~MPITagData()
{
  if (!m_isFree)
    m_manager->free_tag(*this);
}

}  // namespace 

MPITagManager::MPITagManager() :
  m_keyManager(std::make_shared<MPIKeyManager>()),
  m_commCompare(m_keyManager),
  m_tags(m_commCompare)
{
  int isInitialized;
  MPI_Initialized(&isInitialized);
  ThrowRequireMsg(isInitialized, "MPI must be initialized prior to constructing MPITagManager");

  int flag;
  int* val;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &val, &flag);
  ThrowRequireMsg(flag, "This MPI implementation is erroneous");
  ThrowRequireMsg(*val >= m_tagMax, "MPI_TAG_UB must be at least " + std::to_string(m_tagMax));
  m_tagMax = *val;
}


MPITag MPITagManager::get_tag(MPI_Comm comm)
{
  return get_tag(comm, m_tagMin);
}

MPITag MPITagManager::get_tag(MPI_Comm comm, int tagHint)
{
  if (!m_keyManager->has_key(comm))
  {
    m_keyManager->get_key(comm);
    m_keyManager->register_callback(comm, std::bind(&MPITagManager::free_comm_keys, this, std::placeholders::_1));
  }

  auto& commRecord = m_tags[comm];
  auto& tags = commRecord.tags;

  auto it = tags.find(tagHint);

  int newTag = -1;
  if (it == tags.end()) {
    newTag = tagHint;
  } else {
    // find next available tag
    int prevVal = *it;
    while (it != tags.end())
    {
      if ( (*it - prevVal) > 1)
      {
        newTag = prevVal + 1;
        break;
      }
      prevVal = *it;
      it++;
    }

    if (newTag == -1) { // if no spaces between existing tags found
      newTag = *(std::prev(tags.end())) + 1;
    }
  }

  assert(newTag != -1);
  ThrowRequireMsg(newTag <= m_tagMax, "New tag must be <= " + std::to_string(m_tagMax));
  debug_check(comm, newTag);

  auto tagData = std::make_shared<Impl::MPITagData>(this, comm, newTag);
  commRecord.insert(tagData);

  return MPITag(tagData);
}


void MPITagManager::free_tag(MPITag& tag)
{
  free_tag(*(tag.m_data));
}


void MPITagManager::free_tag(Impl::MPITagData& tag)
{
  m_tags.at(tag.get_comm()).erase(tag);
  tag.set_free();
}


//-----------------------------------------------------------------------------
// MPI_Comm key management


void MPITagManager::free_comm_keys(MPI_Comm comm)
{
  ThrowRequireMsg(m_tags.count(comm) == 1, "Cannot free MPI Comm that is not assigned (possible double free)");

  auto& commRecords = m_tags[comm];
  for (auto& weak_tag_ptr : commRecords.tagPtrs) {
    if (auto tag_ptr = weak_tag_ptr.second.lock()) {
      tag_ptr->set_free();
    }
  }

  m_tags.erase(comm);
}

void MPITagManager::debug_check(MPI_Comm comm, int newVal)
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

    ThrowRequireMsg(isSame, "Calls to MPICommManager must be collective");
  }

  MPI_Barrier(comm);
#endif

}


//-----------------------------------------------------------------------------
// CommRecord

void MPITagManager::CommRecord::insert(std::shared_ptr<Impl::MPITagData> new_tag)
{
  int tag_val = new_tag->get_tag();
  ThrowRequireMsg(tags.count(tag_val) == 0, "Cannot create new tag with same value as existing tag");

  tags.insert(tag_val);
  tagPtrs[tag_val] = new_tag;
}

void MPITagManager::CommRecord::erase(Impl::MPITagData& tag)
{
  int tag_val = tag.get_tag();
  ThrowRequireMsg(tags.count(tag_val) == 1, "Cannot free tag this has not been assigned (possible double free");

  tags.erase(tag_val);
  tagPtrs.erase(tag_val);
}



MPITagManager& get_mpi_tag_manager()
{
  static MPITagManager tagManager;
  return tagManager;
}

}
