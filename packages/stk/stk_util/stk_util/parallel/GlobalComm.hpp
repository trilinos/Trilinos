#ifndef STK_UTIL_GLOBAL_COMM_hpp
#define STK_UTIL_GLOBAL_COMM_hpp

#include <mutex>

#ifdef HAVE_MPI

#include <mpi.h>

#else  // HAVE_MPI
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif // HAVE_MPI

namespace stk {

static std::mutex mpi_mutex;
static MPI_Comm Global_STK_Comm = MPI_COMM_WORLD;

inline void initialize_global_comm(MPI_Comm comm)
{
  std::lock_guard<std::mutex> guard(mpi_mutex);
  Global_STK_Comm = comm;
}

inline MPI_Comm get_global_comm()
{
  std::lock_guard<std::mutex> guard(mpi_mutex);
  return Global_STK_Comm;
}

}  // namespace stk

#endif  // STK_UTIL_GLOBAL_COMM_hpp