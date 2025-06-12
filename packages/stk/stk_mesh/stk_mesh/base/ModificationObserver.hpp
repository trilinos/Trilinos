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

#ifndef ModificationObserver_hpp
#define ModificationObserver_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk
{
namespace mesh
{

enum ModificationObserverPriority : int {
  STK_INTERNAL_HIGH_PRIORITY = 0,
  STK_INTERNAL_LOW_PRIORITY  = 1,
  STK_INTERNAL               = STK_INTERNAL_LOW_PRIORITY,
  STK_TRANSITION             = 5,
  APPLICATION                = 9
};


class ModificationObserver
{
public:
    ModificationObserver(ModificationObserverPriority priority)
      : m_priority(priority)
    {}

    virtual ~ModificationObserver()
    {
    }

    void set_priority(ModificationObserverPriority priority) { m_priority = priority; }
    ModificationObserverPriority get_priority() const { return m_priority; }

    virtual void entity_parts_added(stk::mesh::Entity, const stk::mesh::OrdinalVector&)
    {
    }

    virtual void entity_parts_removed(stk::mesh::Entity, const stk::mesh::OrdinalVector&)
    {
    }

    virtual void entity_added(stk::mesh::Entity)
    {
    }

    virtual void entity_deleted(stk::mesh::Entity)
    {
    }

    virtual void modification_begin_notification()
    {
    }
  
    virtual void started_modification_end_notification()
    {
    }

    virtual void finished_modification_end_notification()
    {
    }

    virtual void elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec &)
    {
    }

    virtual void elements_moved_procs_notification(const stk::mesh::EntityProcVec &)
    {
    }

    virtual void local_entities_created_or_deleted_notification(stk::mesh::EntityRank)
    {
    }

    virtual void local_entity_comm_info_changed_notification(stk::mesh::EntityRank)
    {
    }

    virtual void local_buckets_changed_notification(stk::mesh::EntityRank)
    {
    }

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
    {
        valuesToReduce.clear();
    }

    virtual void set_reduced_values(const std::vector<size_t> &)
    {
    }

    virtual void relation_destroyed(Entity, Entity, ConnectivityOrdinal)
    {
    }

    virtual void relation_declared(Entity, Entity, ConnectivityOrdinal)
    {
    }

private:
    ModificationObserverPriority m_priority;

    ModificationObserver() = delete;
};

}
}

#endif
