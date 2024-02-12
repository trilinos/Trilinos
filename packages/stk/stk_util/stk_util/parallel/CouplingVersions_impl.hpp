#ifndef stk_util_parallel_CouplingVersionsImpl_hpp
#define stk_util_parallel_CouplingVersionsImpl_hpp

#include "stk_util/parallel/Parallel.hpp"

namespace stk {
namespace util {

// these functions are for testing only.  Applications should *not* call them
namespace impl {

  void set_coupling_version(MPI_Comm comm, int version);

  void set_error_on_reset(bool val);

  void reset_global_max_coupling_version();

}
}
}

#endif
