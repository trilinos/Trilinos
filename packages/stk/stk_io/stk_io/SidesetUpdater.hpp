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

#ifndef STK_SIDESET_UPDATER_HPP
#define STK_SIDESET_UPDATER_HPP

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequire
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/Env.hpp>
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector
#include "stk_io/StkIoUtils.hpp"

namespace stk { namespace io {

struct part_compare_by_ordinal {
  bool operator() (const stk::mesh::Part * const i, const stk::mesh::Part * const j) {
      if(i == j) return false;
      if(nullptr == i) return true;
      if(nullptr == j) return false;
      return (i->mesh_meta_data_ordinal() < j->mesh_meta_data_ordinal());
  }
};

enum OPERATION {ADDED, DELETED};

class SidesetUpdater : public stk::mesh::ModificationObserver
{
public:

    SidesetUpdater(stk::mesh::BulkData &bulk, stk::mesh::Selector active_selector)
    : bulkData(bulk), activeSelector(active_selector), internalSidesetWarningHasBeenIssued(false)
    {
    }

    bool should_notifier_delete() const override
    {
        return true;
    }

    void fill_sidesets_element_belongs_to(stk::mesh::Entity elem);

    virtual void entity_added(stk::mesh::Entity entity);

    virtual void entity_deleted(stk::mesh::Entity entity);

    void reconstruct_noninternal_sidesets(const std::vector<size_t> &reducedValues);

    void reconstruct_sidesets();

    virtual void finished_modification_end_notification();

    virtual void started_modification_end_notification();

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce);

    virtual void set_reduced_values(const std::vector<size_t> &reducedValues);

    virtual void elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove);

    virtual void elements_moved_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove);

    void tag_sideset(stk::mesh::Entity entity, const stk::mesh::Part& part, OPERATION op);

    void insert_parts(stk::mesh::Entity entity, const stk::mesh::ConstPartVector& parts, OPERATION op);

    void insert_parts(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts, OPERATION op);

    virtual void entity_parts_added(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts);

    virtual void entity_parts_removed(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts);

    void set_active(bool active);

private:
    stk::mesh::BulkData &bulkData;
    stk::mesh::Selector activeSelector;
    std::set<const stk::mesh::Part*, part_compare_by_ordinal> stkSideSets;
    bool isActive = true;
    bool internalSidesetWarningHasBeenIssued;
};

void toggle_sideset_updaters(stk::mesh::BulkData& bulk, bool flag);

}} // end stk mesh namespaces

#endif
