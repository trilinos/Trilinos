/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_THREAD_LOCAL_DATA_H
#define STK_THREAD_LOCAL_DATA_H

#include <stk_util/stk_config.h>
#include <iostream>
#include <vector>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace stk {

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Dec 2025
template<class dataType> class STK_DEPRECATED ThreadLocalData {

 public:

   std::size_t concurrency() const {
#if defined( _OPENMP )
    return omp_get_max_threads();
#else
    return 1;
#endif
  }

  std::size_t hardware_thread_id() const {
#if defined( _OPENMP )
    return omp_get_thread_num();
#else
    return 0;
#endif
  }

  ThreadLocalData()
      : m_threadScratchDataVector(concurrency()) {
  }

  explicit ThreadLocalData(dataType b)
      : m_threadScratchDataVector(concurrency(), b) {
  }

  //
  //  Total size of structure (generally number of threads available)
  //
  unsigned NumThreads() const {
    return m_threadScratchDataVector.size();
  }

  //
  //  Get entry for current thread.  MUST be called in a parallel section
  //
  dataType& getMyThreadEntry() {
    return m_threadScratchDataVector[hardware_thread_id()];
  }

  const dataType& getMyThreadEntry() const {
    return m_threadScratchDataVector[hardware_thread_id()];
  }

 private:
  std::vector<dataType> m_threadScratchDataVector;
};
#endif

}

#endif
