#ifndef stk_util_parallel_CommReplacer
#define stk_util_parallel_CommReplacer

#include <vector>
#include "stk_util/parallel/Parallel.hpp"    // for MPI
#include "stk_util/parallel/MPICommKey.hpp"  // for MPIKeyManager

namespace stk {
namespace impl {

#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
  class CommReplacer
  {
    public:
      CommReplacer();

      MPI_Comm get_copied_comm(MPI_Comm origComm);

      MPI_Comm get_orig_comm(MPI_Comm copyComm);

      void delete_comm_pair(MPI_Comm origComm);

    private:
      std::vector<MPI_Comm> m_origComms;
      std::vector<MPI_Comm> m_copyComms;
      int m_mpiAttrKey;
  };
#endif

}  // namespace
}  // namespace

#endif