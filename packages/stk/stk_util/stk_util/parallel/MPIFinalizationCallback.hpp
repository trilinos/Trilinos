#ifndef stk_util_parallel_MPIFinalizationCallback
#define stk_util_parallel_MPIFinalizationCallback

#include <functional>
#include <cassert>

#include "Parallel.hpp"  // for MPI


namespace stk {


// This class executes a callback when the first of two events occurs:
//   1. The destructor of this class gets called
//   2. MPI_Finalize is called
// The callback will be called exactly once
//
// The way to use this is:
// class NeedToDoCleanupWhenMPIFinalizeIsCalled
// {
//    public:
//      NeedToDoCleanupWhenMPIFinalizeIsCalled() :
//        m_destructor( [=](){cleanup_func();} )
//      {}
//
//      ~NeedToDoCleanupWhenMPIFinalizeIsCalled() { m_destructor.destructor()}
//    private:
//      void cleanup_func();  // does cleanup work
//      MPIFinalizationCallback m_destructor;
// };
// The reason the destructor needs to call m_destructor->destructor() is that
// the callback must run *before* the data members of NeedToDoCleanupWhenMPIFinalizeIsCalled
// are destroyed, which happens immediately after the destructor body executes
class MPIFinalizationCallback
{
  public:
    MPIFinalizationCallback();

    explicit MPIFinalizationCallback(std::function<void()> callback);

    MPIFinalizationCallback(const MPIFinalizationCallback&) = delete;
    MPIFinalizationCallback(MPIFinalizationCallback&& other) = delete;

    MPIFinalizationCallback& operator=(const MPIFinalizationCallback&) = delete;
    MPIFinalizationCallback& operator=(MPIFinalizationCallback&&) = delete;

    ~MPIFinalizationCallback();

    void destructor();

  private:
    std::function<void()> m_callback;
    int m_destructorAttrKey;
    bool m_destructorCalled = false;
};

}

#endif