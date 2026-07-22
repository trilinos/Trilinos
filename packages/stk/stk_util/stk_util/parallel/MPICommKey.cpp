#include <iostream>
#include <cassert>
#include <limits>
#include "stk_util/parallel/MPICommKey.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Parallel.hpp"


namespace stk {

namespace impl {

bool CommCompare::operator()(MPI_Comm comm1, MPI_Comm comm2)
{
    auto key1 = m_manager->get_key(comm1);
    auto key2 = m_manager->get_key(comm2);
    return std::less<MPIKeyManager::CommKey>{}(key1, key2);
}

bool CommCompare::operator()(MPI_Comm comm1, MPI_Comm comm2) const
{
    auto key1 = m_manager->get_key(comm1);
    auto key2 = m_manager->get_key(comm2);
    return std::less<MPIKeyManager::CommKey>{}(key1, key2);
}

namespace impl {

int delete_mpi_comm_key(MPI_Comm comm, int /*comm_keyval*/, void* /*attribute_val*/, void* extra_state)
{
  MPIKeyManager* key_manager = reinterpret_cast<MPIKeyManager*>(extra_state);
  key_manager->free_comm(comm);

  return MPI_SUCCESS;
}


}

MPIKeyManager::MPIKeyManager() : m_destructor([this]() { destructor(); })
{
  int isInitialized;
  MPI_Initialized(&isInitialized);
  STK_ThrowRequireMsg(isInitialized, "MPI must be initialized prior to constructing MPIKeyManager");

  MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, &impl::delete_mpi_comm_key, &m_mpiAttrKey, this);

  int myRank;
  MPI_Comm_rank(parallel_machine_world(), &myRank);
  m_currentCommKey = myRank;
}

MPIKeyManager::~MPIKeyManager()
{
  m_destructor.destructor();
}


MPIKeyManager::CommKey MPIKeyManager::get_key(MPI_Comm comm)
{
  const CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);

  if (!foundFlag) {
    commKey = generate_comm_key();
    MPI_Comm_set_attr(comm, m_mpiAttrKey, const_cast<CommKey*>(commKey));
    m_comms[*commKey] = comm;
  }


  return *commKey;
}

MPIKeyManager::CommKey MPIKeyManager::get_key(MPI_Comm comm) const
{
  const CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);

  STK_ThrowRequireMsg(foundFlag, "Key must be assigned before use in const version of get_key()");

  return *commKey;
}

bool MPIKeyManager::has_key(MPI_Comm comm)
{
  const CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);

  return foundFlag;
}


MPIKeyManager::CallerUID MPIKeyManager::get_UID()
{
  return m_currUID++;
}

void MPIKeyManager::register_callback(MPI_Comm comm, CallerUID uid, Callback func)
{
  assert(uid < m_currUID);
  auto key = get_key(comm);
  m_callbacks[key][uid].push_back(func);
}


void MPIKeyManager::execute_callbacks_immediately(CallerUID uid)
{
  for (auto& commKeyMapPair : m_callbacks)
  {
    MPI_Comm comm            = m_comms[commKeyMapPair.first];
    auto& uidCallbackVecMap  = commKeyMapPair.second;
    if (uidCallbackVecMap.count(uid) == 0)
      continue;

    for (auto& callback : uidCallbackVecMap[uid])
      callback(comm);

    uidCallbackVecMap.erase(uid);
  }
}


void MPIKeyManager::unregister_callbacks(CallerUID uid)
{
  for (auto& comm_map_pair : m_callbacks)
  {
    auto& uidCallbackVecMap = comm_map_pair.second;
    if (uidCallbackVecMap.count(uid) > 0) {
      uidCallbackVecMap.erase(uid);
    }
  }
}


MPIKeyManager::CommKey MPIKeyManager::get_next_comm_key()
{
  auto max_val = std::numeric_limits<CommKey>::max();
  if (m_currentCommKey == max_val)
    throw std::runtime_error(std::string("cannot have more than ") +
            std::to_string(max_val) + " communicators in a program");

  return m_currentCommKey++;
}


const MPIKeyManager::CommKey* MPIKeyManager::generate_comm_key()
{
  auto valFound = get_next_comm_key();

  auto p = m_usedCommKeys.insert(valFound);
  STK_ThrowRequireMsg(p.second, "Error in generateCommKey()");

  return &(*p.first);
}


void MPIKeyManager::free_comm(MPI_Comm comm)
{
  const auto& this_ref = *this;
  auto key = this_ref.get_key(comm);
  STK_ThrowRequireMsg(m_usedCommKeys.count(key) == 1, "Cannot free MPI Comm that is not assigned (possible double free)");
  for (auto& uidCallbackVecPair : m_callbacks[key]) {
    for (auto& callback : uidCallbackVecPair.second) {
      callback(comm);
    }
  }

  m_callbacks.erase(key);
  m_usedCommKeys.erase(key);
  m_comms.erase(key);
}


void MPIKeyManager::destructor()
{
  while (!m_comms.empty()) {
    auto& p = *(m_comms.begin());
    // delete the attribute to trigger callback to free_comm(), and
    // also ensure the callback does *not* get called again when
    // the comm is freed
    MPI_Comm_delete_attr(p.second, m_mpiAttrKey);
  }
}

}

}
