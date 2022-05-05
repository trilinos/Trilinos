#ifndef stk_util_parallel_MPICommKey
#define stk_util_parallel_MPICommKey

#include <functional>
#include <set>
#include <map>
#include <vector>
#include <memory>
#include "Parallel.hpp"  // for MPI

namespace stk {

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


namespace Impl {
  int delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);
}

// Assigns keys to MPI_Comm objects.  In conjunction with CommCompare, allows
// storing MPI_Comms in stl containers that require a comparison function.
// Note: once the communicator is destroyed (ex. by MPI_Comm_free), it is no longer
//       possible to use a container that contains the destroyed MPI_Comm.  Use
//       the register_callback() function to get notified just before the MPI_Comm
//       gets destroyed so you can remove it from your data structures.
class MPIKeyManager
{
  public:

  using CommKey = int;

  using Callback = std::function<void(MPI_Comm)>;

  MPIKeyManager();

  ~MPIKeyManager() = default;

  CommKey get_key(MPI_Comm comm);

  CommKey get_key(MPI_Comm comm) const;


  bool has_key(MPI_Comm comm);

  // registers a callback that will be called just before the MPI_Comm is destroyed.
  // Function must be callable as func(MPI_Comm)
  void register_callback(MPI_Comm comm, Callback func);

  private:

    const CommKey* generate_comm_key();

    void free_comm(MPI_Comm comm);


    int m_mpiAttrKey;
    std::set<CommKey> m_usedCommKeys;
    std::map<MPI_Comm, std::vector<Callback>> m_callbacks;

  friend int Impl::delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);

};

}

#endif
