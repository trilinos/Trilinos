#include "stk_util/parallel/CommReplacer.hpp"
#include "stk_util/util/ReportHandler.hpp"   // for ThrowAssertMsg, ThrowRequire

namespace stk {
namespace impl {
      
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN

CommReplacer::CommReplacer()
{
  MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN, &m_mpiAttrKey, nullptr);
}


MPI_Comm CommReplacer::get_copied_comm(MPI_Comm origComm)
{
  int* unused = nullptr;
  int foundFlag;
  MPI_Comm_get_attr(origComm, m_mpiAttrKey, &unused, &foundFlag);

  if (!foundFlag) {
    delete_comm_pair(origComm);
    MPI_Comm_set_attr(origComm, m_mpiAttrKey, nullptr);


    MPI_Comm copyComm;
    MPI_Comm_dup(origComm, &copyComm);
    m_origComms.push_back(origComm);
    m_copyComms.push_back(copyComm);

    return copyComm;
  } else {
    for (unsigned int i=0; i < m_origComms.size(); ++i) {
      if (m_origComms[i] == origComm) {
        return m_copyComms[i];
      }
    }

    STK_ThrowRequireMsg(false, "Unable to find origComm");
  }

  return MPI_COMM_NULL;
}


MPI_Comm CommReplacer::get_orig_comm(MPI_Comm copyComm)
{
  for (unsigned int i=0; i < m_copyComms.size(); ++i) {
    if (m_copyComms[i] == copyComm) {
      return m_origComms[i];
    }
  }

  STK_ThrowRequireMsg(false, "Unable to find copyComm");
  return MPI_COMM_NULL;
}


void CommReplacer::delete_comm_pair(MPI_Comm origComm)
{
  for (unsigned int i=0; i < m_origComms.size(); ++i) {
    if (origComm == m_origComms[i])
    {
      MPI_Comm_free(&(m_copyComms[i]));
      auto itOrig = m_origComms.begin();
      auto itCopy = m_copyComms.begin();
      std::advance(itOrig, i);
      std::advance(itCopy, i);

      m_origComms.erase(itOrig);
      m_copyComms.erase(itCopy);
    }
  }
}

#endif

} // namespace
} // namespace