#ifndef STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_
#define STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_

#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <string>

namespace stk {

void log_with_time_and_memory(MPI_Comm communicator, const std::string &message);

}

#endif /* STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_ */
