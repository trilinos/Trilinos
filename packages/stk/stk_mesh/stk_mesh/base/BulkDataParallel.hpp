// Copyright (c) 2014, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// (INCLUDING NEGLIGENCE OR OTHER ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef stk_mesh_BulkDataParallel_hpp
#define stk_mesh_BulkDataParallel_hpp

#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine

namespace stk {
namespace mesh {

class BulkDataParallel {

public: //methods
  BulkDataParallel(ParallelMachine parallel_in) :
    m_parallel_machine( parallel_in ),
    m_parallel_size( parallel_machine_size( parallel_in ) ),
    m_parallel_rank( parallel_machine_rank( parallel_in ) )
  {
    //no-op constructor right now
  }

  /** \brief  The parallel machine */
  ParallelMachine parallel() const { return m_parallel_machine ; }

  /** \brief  Size of the parallel machine */
  int parallel_size()   const { return m_parallel_size ; }

  /** \brief  Rank of the parallel machine's local processor */
  int parallel_rank()   const { return m_parallel_rank ; }

  ~BulkDataParallel() {}  //do nothing destructor

private: //data
  ParallelMachine m_parallel_machine;
  int m_parallel_size;
  int m_parallel_rank;

}; //class 

} //namespace mesh
} //namespace stk

#endif // stk_mesh_BulkDataParallel_hpp
