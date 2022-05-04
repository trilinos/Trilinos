
#include "stk_util/parallel/MPICommKey.hpp"
#include "ParallelComm.hpp"

namespace stk {

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

namespace Impl {

int delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state)
{
  //auto comm_key_ptr          = reinterpret_cast<MPIKeyManager::CommKey*>(attribute_val);
  MPIKeyManager* key_manager = reinterpret_cast<MPIKeyManager*>(extra_state);
  // the predefined comms sometimes (but not allways) get freed after the
  // MPICommManager static variable, so we don't need to (and can't)
  // unregister them
  bool isSpecial = comm == MPI_COMM_WORLD ||
                   comm == MPI_COMM_SELF  ||
                   comm == MPI_COMM_NULL;
  if (!isSpecial) {
    key_manager->free_comm(comm);
  }

  return MPI_SUCCESS;
}

}


MPIKeyManager::MPIKeyManager()
{
  int isInitialized;
  MPI_Initialized(&isInitialized);
  ThrowRequireMsg(isInitialized, "MPI must be initialized prior to constructing MPIKeyManager");

  MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, &Impl::delete_mpi_comm_key, &m_mpiAttrKey, this);
}


MPIKeyManager::CommKey MPIKeyManager::get_key(MPI_Comm comm)
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

MPIKeyManager::CommKey MPIKeyManager::get_key(MPI_Comm comm) const
{
  const CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);

  ThrowRequireMsg(foundFlag, "Key must be assigned before use in const version of get_key()");

  return *commKey;
}

bool MPIKeyManager::has_key(MPI_Comm comm)
{
  const CommKey* commKey;
  int foundFlag;
  MPI_Comm_get_attr(comm, m_mpiAttrKey, &commKey, &foundFlag);

  return foundFlag;
}


void MPIKeyManager::register_callback(MPI_Comm comm, Callback func)
{
  get_key(comm);
  m_callbacks[comm].push_back(func);
}


const MPIKeyManager::CommKey* MPIKeyManager::generate_comm_key()
{
  // find first unused key in range [0, int_max]
  CommKey valPrev = -1, valFound = -1;
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


void MPIKeyManager::free_comm(MPI_Comm comm)
{
  ThrowRequireMsg(m_usedCommKeys.count(get_key(comm)) == 1, "Cannot free MPI Comm that is not assigned (possible double free)");
  for ( auto& callback : m_callbacks[comm])
    callback(comm);

  m_callbacks.erase(comm);
  m_usedCommKeys.erase(get_key(comm));
  
}

}