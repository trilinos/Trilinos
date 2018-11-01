// Copyright (c) 2013, Sandia Corporation.
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

class ModificationObserver
{
public:
    virtual ~ModificationObserver()
    {
    }

    virtual bool should_notifier_delete() const
    {
        return false;
    }

    virtual void entity_parts_added(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts)
    {
    }

    virtual void entity_parts_removed(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts)
    {
    }

    virtual void entity_added(stk::mesh::Entity entity)
    {
    }

    virtual void entity_deleted(stk::mesh::Entity entity)
    {
    }

    virtual void started_modification_end_notification()
    {
    }

    virtual void finished_modification_end_notification()
    {
    }

    virtual void elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
    }

    virtual void elements_moved_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
    }

    virtual void local_entities_created_or_deleted_notification(stk::mesh::EntityRank rank)
    {
    }

    virtual void local_entity_comm_info_changed_notification(stk::mesh::EntityRank rank)
    {
    }

    virtual void local_buckets_changed_notification(stk::mesh::EntityRank rank)
    {
    }

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
    {
        valuesToReduce.clear();
    }

    virtual void set_reduced_values(const std::vector<size_t> &reducedValues)
    {
    }
};

}
}

#endif
