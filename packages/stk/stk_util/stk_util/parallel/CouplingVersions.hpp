#ifndef stk_util_parallel_CouplingVersions_hpp
#define stk_util_parallel_CouplingVersions_hpp

#include "stk_util/stk_config.h"
#ifdef STK_HAS_MPI
#include "mpi.h"  // dont' include Parallel.hpp to avoid circular dependency
#endif

#include <string>

#define STK_MAX_COUPLING_VERSION 14
#define STK_MIN_COUPLING_VERSION 0                                                                               
 
namespace stk {

namespace util {

namespace impl {
constexpr int SHORT_TERM_STK_MAX_COUPLING_VERSION=1;
}


int get_common_coupling_version();

int get_local_max_coupling_version();

int get_local_min_coupling_version();

int get_global_max_coupling_version();

std::string get_deprecation_date(int version);

#ifdef STK_HAS_MPI
void set_coupling_version(MPI_Comm comm);
#endif

bool is_local_stk_coupling_deprecated();

void print_unsupported_version_warning(int version, int line, const char* file);

}

}

#endif
