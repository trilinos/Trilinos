#include "stk_util/parallel/MPITag.hpp"
#include "stk_util/parallel/MPITagManager.hpp"

namespace stk {
namespace impl {

MPITagData::~MPITagData()
{
  if (!m_isFree)
    m_manager->free_tag_local(*this);
}

MPIKeyManager::CommKey MPITagData::get_comm_key()
{ 
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    return m_manager->m_keyManager->get_key(m_commInternal); 
#else
    return m_manager->m_keyManager->get_key(m_comm); 
#endif
}

}  // namespace
}  // namespace