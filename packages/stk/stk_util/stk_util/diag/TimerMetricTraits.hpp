// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_diag_TimerMetricTraits_hpp
#define stk_util_diag_TimerMetricTraits_hpp

#include <string>

namespace stk {
namespace diag {

typedef unsigned long MetricsMask;

/**
 * @brief Member function <code>setTimerTimeFormat</code> sets the display format of time in the
 * output tables
 *
 * @param time_format    a <code>TimeFormat</code> variable...
 */
void setTimerTimeFormat(int time_format);

/**
 * @brief Member function <code>getTimerTimeFormat</code>
 *
 */
int getTimerTimeFormat();

/**
 * @brief Enumeration <code>Metrics</code> assigns a unique but for each type of timer.  The
 * METRICS_FORCE type allows a timer to force itself to be active, even is not enabled.
 *
 */
enum Metrics {
  METRICS_LAP_COUNT      = 0x0001,              ///< Count of timer starts
  METRICS_CPU_TIME       = 0x0002,              ///< CPU runtime 
  METRICS_WALL_TIME      = 0x0004,              ///< Wall clock time
  METRICS_MPI_COUNT      = 0x0008,              ///< MPI class count
  METRICS_MPI_BYTE_COUNT = 0x0010,              ///< MPI byte count
  METRICS_HEAP_ALLOC     = 0x0020,              ///< Heap allocation
  METRICS_ALL            = 0x7FFF,

  METRICS_FORCE          = 0x8000               ///< Force metrics to be acquired
};


struct LapCount {};                             ///< Lap counter metric tag
struct CPUTime {};                              ///< CPU runtime metric tag
struct WallTime {};                             ///< Wall clock metric tag
struct MPICount {};                             ///< MPI call count metric tag
struct MPIByteCount {};                         ///< MPI byte count metric tag
struct HeapAlloc {};                            ///< Heap allocation metric tag

template <class T>
struct MetricTraits;

template<>
struct MetricTraits<LapCount>
{
  typedef unsigned Type;
  enum {METRIC = METRICS_LAP_COUNT};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type time);
};

template<>
struct MetricTraits<CPUTime>
{
  typedef double Type;
  enum {METRIC = METRICS_CPU_TIME};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type time);
};

template<>
struct MetricTraits<WallTime>
{
  typedef double Type;
  enum {METRIC = METRICS_WALL_TIME};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type time);
};

template<>
struct MetricTraits<MPICount>
{
  typedef double Type;
  enum {METRIC = METRICS_MPI_COUNT};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type count);
};

template<>
struct MetricTraits<MPIByteCount>
{
  typedef double Type;
  enum {METRIC = METRICS_MPI_BYTE_COUNT};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type count);
};


template<>
struct MetricTraits<HeapAlloc>
{
  typedef double Type;
  enum {METRIC = METRICS_HEAP_ALLOC};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type count);
};

template <class T>
typename MetricTraits<T>::Type now() {
  return MetricTraits<T>::value_now();
}

} // namespace diag
} // namespace stk

#endif // stk_util_diag_TimerMetricTraits_hpp
