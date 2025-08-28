/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
 
#ifndef STK_TRANSFER_UTIL_LOG_MESSAGE_HPP
#define STK_TRANSFER_UTIL_LOG_MESSAGE_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/environment/LogWithTimeAndMemory.hpp>

namespace stk {
namespace transfer_util {

inline
void log_message(MPI_Comm comm, const std::string& msg)
{
  stk::log_with_time_and_memory(comm, msg, stk::outputP0());
}

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_LOG_MESSAGE_HPP
