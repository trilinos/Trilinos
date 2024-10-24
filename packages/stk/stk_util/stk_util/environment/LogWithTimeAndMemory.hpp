#ifndef STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_
#define STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_

#include "stk_util/parallel/Parallel.hpp"  // for MPI_Comm
#include "stk_util/environment/Env.hpp"
#include <string>                          // for string

namespace stk {

void log_with_time_and_memory(MPI_Comm communicator, const std::string &message, std::ostream& ostrm = sierra::Env::outputP0());

}

#endif /* STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_ */
