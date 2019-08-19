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
#ifndef PHX_ALLOCATION_TRACKER_HPP
#define PHX_ALLOCATION_TRACKER_HPP

#include "Phalanx_config.hpp"
#include <forward_list>
#include <memory>

namespace PHX {

  /** \brief This object is siimilar to a memory pool in that allows
   *  for reuse of view allocations across the DAG and in other
   *  DataContainers and FieldManagers.
   *
   *  A field may only be used in a small section of the DAG. After
   *  topological sorting, we can find the span of evaluators in the
   *  sorted list that the field must exist over. Outside of this
   *  range, the view memory can be reused by other views that don't
   *  overlap within the same DAG.
   *
   *  An additional feature is that since only one evaluation type of
   *  one FieldManager is run at a time, phalanx can also reuse view
   *  allocations for different evaluation types in the same
   *  FieldManager and over all evaluation types in other
   *  FieldManagers. A special clone method exists that creates a new
   *  MemoryManager, pointing to the same allocations, but resetting the
   *  tracker objects for a new FieldManager or DataContainer.
   */
   class  MemoryManager {

    struct Allocation {
      /// Size of the allocation.
      std::size_t size_;
      /// Evaluator range where this allocation is being used.
      std::vector<std::pair<int,int>> use_ranges_;
      /// A reference counted memory allocation for a view.
      Kokkos::Impl::SharedAllocationTracker tracker_;
    };

    std::forward_list<PHX::MemoryManager::Allocation> allocations_;

    /// Tracks cloned MemoryManagers so that all clones can get
    /// access to newly allocated fields from any individual memory
    /// pool.
    std::shared_ptr<std::forward_list<PHX::MemoryManager*>> trackers_;

    void findMemoryAllocation() {}

  public:
    MemoryManager()
    {
      trackers_ = std::make_shared<std::forward_list<PHX::MemoryManager*>>();
      trackers_->push_front(this);
    }

    ~MemoryManager()
    {
      // remove self from list of memory pools
      trackers_->remove(this);
    }

    /// Allocate a new memory pool re-using allocations from other linked MemoryManagers.
    MemoryManager(const MemoryManager& mp)
    {
      trackers_ = mp.trackers_;
      trackers_ = mp.trackers_;
      trackers_->push_front(this);
    }

    /** \brief Clones MemoryManager to reuse tracker allocations with a separate FieldManager. */
    std::shared_ptr<PHX::MemoryManager> clone() const
    {return std::make_shared<PHX::MemoryManager>(*this);}

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
