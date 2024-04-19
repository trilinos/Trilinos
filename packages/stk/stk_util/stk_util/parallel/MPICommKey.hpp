#ifndef stk_util_parallel_MPICommKey
#define stk_util_parallel_MPICommKey

#include <functional>
#include <set>
#include <map>
#include <vector>
#include <deque>
#include <memory>
#include "Parallel.hpp"  // for MPI
#include "MPIFinalizationCallback.hpp"

// Detect if MPI has bug
#ifdef MPICH
  #define MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
#endif

namespace stk {

//TODO: get rid of this once MPICH/IntelMPI is fixed
namespace impl {

class MPIKeyManager;

class CommCompare
{
  public:
    explicit CommCompare(std::shared_ptr<MPIKeyManager> manager) :
      m_manager(manager)
    {}

    bool operator()(MPI_Comm comm1, MPI_Comm comm2);

    bool operator()(MPI_Comm comm1, MPI_Comm comm2) const;

  private:
    std::shared_ptr<MPIKeyManager> m_manager;
};


namespace impl {
  int delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);

  int destruct_mpi_key_manager(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);
}

// Assigns keys to MPI_Comm objects.  In conjunction with CommCompare, allows
// storing MPI_Comms in stl containers that require a comparison function.
// MPI_Comms sort in the same order on all processes.
// The key value is *not* guaranteed to be the same on all processes.
// Note: once the communicator is destroyed (ex. by MPI_Comm_free), it is no longer
//       possible to use a container that contains the destroyed MPI_Comm.  Use
//       the register_callback() function to get notified just before the MPI_Comm
//       gets destroyed so you can remove it from your data structures.
//       Callbacks will also be called when MPI_Finalize is called or when
//       the MPIKeyManager destructor runs.  The latter case can be a
//       problem if the the callback requires another object to still
//       be alive.  Use execute_callbacks_immediately() or
//       unregister_callbacks() to handle this.  One particular case
//       of this problem is:
//
//       class Foo
//       {
//         public:
//           Foo() :
//             m_keyManager(std::make_shared<MPIKeyManager>()),
//             m_map(CommCompare(m_keyManager))
//           {
//             m_callerUID = m_keyManager.get_UID();
//             m_keyManager.register_callback(..., m_callerUID, ...);
//           }
//
//           ~Foo()
//           {
//             m_keyManager.execute_callbacks_immediately(m_callerUID);
//           }
//
//
//
//         private:
//           std::shared_ptr<MPIKeyManager> m_keyManager;
//           std::map<MPI_Comm, int, CommCompare> m_map;
//           MPIKeyManager::CallerUID m_callerUID;
//       };
//
//       Notice the call to execute_callbacks_immediately() in the destructor.
//       Without this, the callbacks would execute when the m_keyManager
//       member is destructed, which happens *after* m_map is destructed.
//       So if the callback uses m_map, the program will have undefined
//       behavior.  For this reason, it is recommended to put a call
//       to either execute_callbacks_immediately() or unregister_callbacks()
//       in the destructor body.
class MPIKeyManager
{
  public:

  using CommKey = unsigned long long;

  using CallerUID = unsigned long long;

  using Callback = std::function<void(MPI_Comm)>;

  MPIKeyManager();

  ~MPIKeyManager();

  MPIKeyManager(const MPIKeyManager&) = delete;

  MPIKeyManager& operator=(const MPIKeyManager&) = delete;

  MPIKeyManager(MPIKeyManager&& ) = delete;

  MPIKeyManager& operator=(MPIKeyManager&&) = delete;

  CommKey get_key(MPI_Comm comm);

  CommKey get_key(MPI_Comm comm) const;

  bool has_key(MPI_Comm comm);

  CallerUID get_UID();

  // registers a callback that will be called just before the MPI_Comm is destroyed, or
  // MPI_Finalize is called, or this class gets destructed.
  // Function must be callable as func(MPI_Comm)
  void register_callback(MPI_Comm comm, CallerUID uid, Callback func);

  void execute_callbacks_immediately(CallerUID uid);

  void unregister_callbacks(CallerUID uid);

  private:
    CommKey get_next_comm_key();

    const CommKey* generate_comm_key();

    void free_comm(MPI_Comm comm);

    void destructor();

    int m_mpiAttrKey;
    CommKey m_currentCommKey = 0;
    CallerUID m_currUID   = 0;
    std::set<CommKey> m_usedCommKeys;
    std::map<CommKey, std::map<CallerUID, std::vector<Callback>>> m_callbacks;
    std::map<CommKey, MPI_Comm> m_comms;
    MPIFinalizationCallback m_destructor;

  friend int impl::delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);
};

}  // namespace

}  // namespace

#endif
