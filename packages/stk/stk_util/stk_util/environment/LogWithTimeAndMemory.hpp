#ifndef STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_
#define STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_
#include <mpi.h>
#include <string>

namespace stk {

void log_with_time_and_memory(MPI_Comm communicator, const std::string &message);

}

#endif /* STK_STK_UTIL_STK_UTIL_ENVIRONMENT_LOGWITHTIMEANDMEMORY_HPP_ */
