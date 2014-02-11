#ifndef STK_MEMORY_UTIL_H
#define STK_MEMORY_UTIL_H

#include <string>
#include <mpi.h>

#include <stk_util/parallel/MPI.hpp>

namespace stk {

/*
 * Given a number of bytes, this produces a string of the form
 *    123 B or 123 KB or 123 MB or 123 GB
 * depending on the input number of bytes
 */

std::string human_bytes(size_t bytes);

/*
 *  return current memory usage ('now', resident set size) and high-water-mark ('hwm') in bytes.
 *
 *  Note: Only systems with PROCFS (i.e., Linux) report the high-water-mark of the
 *  resident set size for a process. For those systems, *hwm may be larger than any
 *  observed value of *now. On all other systems, *hwm is the maximum value of *now
 *  that was observed through invoking this function so the true high-water-mark may
 *  be larger than that reported here, depending on your timing.
 */
void get_memory_usage(size_t & now, size_t & hwm);

void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg);

template <typename T>
inline
void get_max_min_avg(MPI_Comm comm, T this_proc, T& max, T& min, T& avg)
{
  int num_procs = 1;
  MPI_Comm_size(comm, &num_procs);

  T this_proc_average = this_proc / num_procs;

  MPI_Datatype type = sierra::MPI::Datatype<T>::type();

  MPI_Allreduce(&this_proc,         &max, 1, type, MPI_MAX, comm);
  MPI_Allreduce(&this_proc,         &min, 1, type, MPI_MIN, comm);
  MPI_Allreduce(&this_proc_average, &avg, 1, type, MPI_SUM, comm);
}

}

#endif /* STK_MEMORY_UTIL_H */
