// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER
#ifndef PHX_MEMORY_POOL_HPP
#define PHX_MEMORY_POOL_HPP

#include "Phalanx_config.hpp"
#include <forward_list>
#include <memory>

namespace PHX {

  class MemoryPool {

    struct Tracker {
      bool is_used_;
      Kokkos::Impl::SharedAllocationTracker tracker_;
    };

    std::forward_list<PHX::MemoryPool::Tracker> trackers_;

    /// Tracks cloned memory pools so that all clones can get access to newly allocated fields from any individual memory pool.
    std::shared_ptr<std::forward_list<PHX::MemoryPool*>> shared_memory_pools_;

    void findMemoryAllocation() {}

  public:
    MemoryPool()
    {
      shared_memory_pools_ = std::make_shared<std::forward_list<PHX::MemoryPool*>>();
      shared_memory_pools_->push_front(this);
    }

    ~MemoryPool()
    {
      // remove self from list of memory pools
      shared_memory_pools_->remove(this);
    }

    /// Allocate a new memory pool re-using trackers from shared memory pools.
    MemoryPool(const MemoryPool& mp)
    {
      shared_memory_pools_ = mp.shared_memory_pools_;
      trackers_ = mp.trackers_;
      shared_memory_pools_->push_front(this);
    }

    /** \brief Clones MemoryPool to reuse tracker allocations with a separate FieldManager. */
    std::shared_ptr<PHX::MemoryPool> clone() const
    {return std::make_shared<PHX::MemoryPool>(*this);}

    /// Assigns memory to a view, allocates new memory if needed.
    template<class View>
    void bindViewMemory(const PHX::FieldTag& tag, View& view) {
      //const size_t bytes = view.span();

      // Find an unused memory allocation.


      // loop over other pools and register as free memory

    }

    /// Inserts tracker
    void insertTracker(Kokkos::Impl::SharedAllocationTracker& t) {

    }

  };

}

#endif
