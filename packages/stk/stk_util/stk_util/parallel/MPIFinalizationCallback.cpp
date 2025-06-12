#include "MPIFinalizationCallback.hpp"
#include "stk_util/parallel/MPIFinalizationCallback.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {

int mpi_comm_destructor(MPI_Comm /*comm*/, int /*comm_keyval*/, void* attribute_val, void* /*extra_state*/)
{
  MPIFinalizationCallback* callback = reinterpret_cast<MPIFinalizationCallback*>(attribute_val);
  callback->destructor();

  return MPI_SUCCESS;
}

MPIFinalizationCallback::MPIFinalizationCallback() :
  MPIFinalizationCallback(std::function<void()>())
{}


MPIFinalizationCallback::MPIFinalizationCallback(std::function<void()> callback) :
  m_callback(callback)
{
  if (m_callback) {
    int isInitialized;
    MPI_Initialized(&isInitialized);
    STK_ThrowRequireMsg(isInitialized, "MPI must be initialized prior to constructing MPIFinalizationCallback");

    // The deleter function will be called when MPI_Finalize is invoked.  This
    // is the recommended way to execute callbacks, see Section 8.7.1 of the
    // MPI 3.1 Standard.
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, &mpi_comm_destructor, &m_destructorAttrKey, nullptr);
    MPI_Comm_set_attr(MPI_COMM_SELF, m_destructorAttrKey, this);
  }
}


MPIFinalizationCallback::~MPIFinalizationCallback()
{
  if (m_callback) {
    int isMpiFinalized = false;
    MPI_Finalized(&isMpiFinalized);
    if (isMpiFinalized) {
      assert(m_destructorCalled);
    }

    if (!isMpiFinalized)
    {
      // call destructor() via the MPI callback, and also ensure destructor()
      // doesn't get called again when MPI_Finalize is called
      MPI_Comm_delete_attr(MPI_COMM_SELF, m_destructorAttrKey);
    }
  }
}

void MPIFinalizationCallback::destructor() 
{ 
  if (!m_destructorCalled && m_callback)
  {
    m_callback(); 
    m_destructorCalled = true;
  }
}

}