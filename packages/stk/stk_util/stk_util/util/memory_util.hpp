#ifndef STK_MEMORY_UTIL_H
#define STK_MEMORY_UTIL_H

#include <string>

namespace stk {

/*
 * Given a number of bytes, this produces a string of the form
 *    123 B or 123 KB or 123 MB or 123 GB
 * depending on the input number of bytes
 */

std::string human_bytes(size_t bytes);

/*
 *  return current memory usage ('now', resident set size) and high-water-mark ('hwm') in megabytes.
 *
 *  Note: Only systems with PROCFS (i.e., Linux) report the high-water-mark of the
 *  resident set size for a process. For those systems, *hwm may be larger than any
 *  observed value of *now. On all other systems, *hwm is the maximum value of *now
 *  that was observed through invoking this function so the true high-water-mark may
 *  be larger than that reported here.
 */
void get_memory_usage(size_t & now, size_t & hwm);

}

#endif /* STK_MEMORY_UTIL_H */
