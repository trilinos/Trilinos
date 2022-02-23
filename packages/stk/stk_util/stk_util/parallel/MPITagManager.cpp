#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, ThrowRequire
#include <cassert>

namespace stk {

namespace Impl {

MPITagData::~MPITagData()
{
  if (!m_isFree)
    m_manager->free_tag(*this);
}

int delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state)
{
  auto comm_key_ptr          = reinterpret_cast<MPITagManager::CommKey*>(attribute_val);
  MPITagManager* tag_manager = reinterpret_cast<MPITagManager*>(extra_state);
  // the predefined comms sometimes (but not allways) get freed after the
  // MPICommManager static variable, so we don't need to (and can't)
  // unregister them
  bool isSpecial = comm == MPI_COMM_WORLD ||
                   comm == MPI_COMM_SELF  ||
                   comm == MPI_COMM_NULL;
  if (!isSpecial) {
    tag_manager->free_comm_key(comm_key_ptr);
  }

  return MPI_SUCCESS;
}

}  // namespace 

MPITagManager::MPITagManager()
{
  int isInitialized;
  MPI_Initialized(&isInitialized);
  ThrowRequireMsg(isInitialized, "MPI must be initialized prior to constructing MPITagManager");

  MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, &Impl::delete_mpi_comm_key, &m_mpiAttrKey, this);

  int flag;
  int* val;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &val, &flag);
  m_tagMax = *val;
  ThrowRequireMsg(flag, "This MPI implementation is erroneous");
}


MPITag MPITagManager::get_tag(MPI_Comm comm)
{
  return get_tag(comm, m_tagMin);
}

//TODO: make CommRecord non copyable
MPITag MPITagManager::get_tag(MPI_Comm comm, int tagHint)
{
  CommKey key = get_comm_key(comm);
  auto& commRecord = m_tags[key];
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
  ThrowRequireMsg(tags.count(newTag) == 0, "Error in get_tag(): selected tag is not unique");
  debug_check(comm, newTag);

  auto tagData = std::make_shared<Impl::MPITagData>(this, key, newTag);
  tags.insert(newTag);
  commRecord.tagPtrs[newTag] = tagData;

  return MPITag(tagData);
}


void MPITagManager::free_tag(MPITag& tag)
{
  free_tag(*(tag.m_data));
}


void MPITagManager::free_tag(Impl::MPITagData& tag)
{
  if (m_tags.find(tag.get_comm_key()) == m_tags.end()) {
    return;
  }

  auto& commRecord = m_tags.at(tag.get_comm_key());
  auto it = commRecord.tags.find(tag.get_tag());

  ThrowRequireMsg(it != commRecord.tags.end(), "Cannot free tag that is not assigned (possible double free)");
  commRecord.tags.erase(it);
  commRecord.tagPtrs.erase(tag.get_tag());
  tag.set_free();
}


//-----------------------------------------------------------------------------
// MPI_Comm key management

MPITagManager::CommKey MPITagManager::get_comm_key(MPI_Comm comm)
{
  const CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);

  if (!foundFlag) {
    commKey = generate_comm_key();
    MPI_Comm_set_attr(comm, m_mpiAttrKey, const_cast<CommKey*>(commKey));
  }

  return *commKey;
}

/*
MPITagManager::CommKey MPITagManager::get_comm_key(MPI_Comm comm)
{
  CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);
  ThrowRequireMsg(foundFlag, "MPI Comm key should already have been set");

  return *commKey;
}
*/


const MPITagManager::CommKey* MPITagManager::generate_comm_key()
{
  // find first unused key in range [0, int_max]
  int valPrev = -1, valFound = -1;
  if (m_usedCommKeys.size() == 0) {
    valFound = 0;
  } else {
    for (auto v : m_usedCommKeys)
    {
      if (v - valPrev > 1) {
        valFound = valPrev + 1;
        break;
      }

      valPrev = v;
    }
  }

  if (valFound == -1) {
    valFound = *(std::prev(m_usedCommKeys.end())) + 1;
  }

  ThrowRequireMsg(valFound >= 0, "generate_comm_key() failed or overflowed");

  auto p = m_usedCommKeys.insert(valFound);
  ThrowRequireMsg(p.second, "Error in generateCommKey()");

  return &(*p.first);
}


void MPITagManager::free_comm_key(CommKey* key)
{
  ThrowRequireMsg(m_tags.count(*key) == 1, "Cannot free MPI Comm that is not assigned (possible double free)");

  auto& commRecords = m_tags.at(*key);
  for (auto& pw : commRecords.tagPtrs) {
    if (auto p = pw.second.lock()) {
      p->set_free();
    }
  }

  m_tags.erase(*key);
  m_usedCommKeys.erase(*key);
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



MPITagManager& get_mpi_tag_manager()
{
  static MPITagManager tagManager;
  return tagManager;
}

}