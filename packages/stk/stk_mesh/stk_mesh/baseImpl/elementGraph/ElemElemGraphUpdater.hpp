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

#ifndef STK_ELEM_ELEM_GRAPH_UPDATER_HPP
#define STK_ELEM_ELEM_GRAPH_UPDATER_HPP

#include <stddef.h>                     // for size_t
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity

#include "ElemElemGraph.hpp"

namespace stk { namespace mesh {

inline
void remove_created_and_deleted_elems(const stk::mesh::BulkData& bulk,
                                      stk::mesh::EntityVector& elems)
{
  stk::util::sort_and_unique(elems);
  elems.erase(
      std::remove_if(elems.begin(), elems.end(),
                     [&](stk::mesh::Entity elem) {
                       return !bulk.is_valid(elem) || bulk.state(elem) == Created;
                     }),
      elems.end()
  );
}

class ElemElemGraphUpdater : public stk::mesh::ModificationObserver
{
public:
    ElemElemGraphUpdater(stk::mesh::BulkData &bulk, stk::mesh::ElemElemGraph &elemGraph_)
    : ModificationObserver(ModificationObserverPriority::STK_INTERNAL_HIGH_PRIORITY),
      bulkData(bulk), elemGraph(elemGraph_), changeEntityOwnerInProgress(false)
    {
    }

    virtual void entity_added(stk::mesh::Entity entity)
    {
        if(bulkData.entity_rank(entity) == stk::topology::ELEM_RANK && !changeEntityOwnerInProgress)
        {
            elementsAdded.push_back(entity);
        }
    }

    virtual void entity_deleted(stk::mesh::Entity entity)
    {
        if(bulkData.entity_rank(entity) == stk::topology::ELEM_RANK && bulkData.bucket(entity).owned())
        {
            elementsDeleted.push_back({entity, bulkData.identifier(entity), bulkData.bucket(entity).topology()});
        }
    }

    virtual void relation_destroyed(Entity from, Entity to, ConnectivityOrdinal /*ordinal*/)
    {
      if (bulkData.entity_rank(from) == stk::topology::ELEM_RANK &&
          bulkData.entity_rank(to) == stk::topology::NODE_RANK &&
          bulkData.bucket(from).owned()) {
        elementsChangedConn.push_back(from);
      }
    }

    virtual void relation_declared(Entity from, Entity to, ConnectivityOrdinal /*ordinal*/)
    {
      if (bulkData.entity_rank(from) == stk::topology::ELEM_RANK &&
          bulkData.entity_rank(to) == stk::topology::NODE_RANK &&
          bulkData.bucket(from).owned()) {
        elementsChangedConn.push_back(from);
      }
    }

    virtual void finished_modification_end_notification()
    {
        if (changeEntityOwnerInProgress) {
            elementsAdded.clear();
            elementsDeleted.clear();
            elementsChangedConn.clear();
            elemGraph.fill_from_mesh();
            changeEntityOwnerInProgress = false;
            return;
        }

        if (maxNumAdded > 0) {
            elemGraph.add_elements(elementsAdded);
            elementsAdded.clear();
        }

        if (maxNumChanged > 0) {
            elemGraph.update_graph_edges(elementsChangedConn);
            elementsChangedConn.clear();
        }
    }

    virtual void started_modification_end_notification()
    {
        stk::util::sort_and_unique(elementsDeleted);
        remove_created_and_deleted_elems(bulkData, elementsChangedConn);

        constexpr size_t numValues = 2;
        size_t localValues[numValues] = {0}, globalValues[numValues] = {0};
        localValues[0] = elementsDeleted.size();
        localValues[1] = elementsChangedConn.size();
        stk::all_reduce_sum(bulkData.parallel(), localValues, globalValues, numValues);

        if (globalValues[0] > 0) {
            elemGraph.delete_elements(elementsDeleted);
            elementsDeleted.clear();
        }

        if (globalValues[1] > 0) {
            elemGraph.update_graph_edges(elementsChangedConn);
            elementsChangedConn.clear();
        }
    }

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
    {
        valuesToReduce.clear();
        stk::util::sort_and_unique(elementsAdded);
        remove_created_and_deleted_elems(bulkData, elementsChangedConn);
        unsigned value = any_added_elements_are_owned(elementsAdded) ? 1 : 0;
        if (value == 0) {
          elementsAdded.clear();
        }
        valuesToReduce.push_back(value);
        valuesToReduce.push_back(elementsChangedConn.size());
    }

    virtual void set_reduced_values(const std::vector<size_t> &reducedValues)
    {
        maxNumAdded = reducedValues[0];
        maxNumChanged = reducedValues[1];
    }

    virtual void elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec & /*elemProcPairsToMove*/)
    {
        changeEntityOwnerInProgress = true;
    }

private:
    bool any_added_elements_are_owned(stk::mesh::EntityVector& elems)
    {
      for(Entity& elem : elems) {
        if (bulkData.is_valid(elem) && bulkData.bucket(elem).owned()) {
          return true;
        }
      }
      return false;
    }

    stk::mesh::BulkData &bulkData;
    stk::mesh::ElemElemGraph &elemGraph;
    stk::mesh::impl::ParallelGraphInfo newParallelGraphEntries;
    stk::mesh::EntityVector elementsAdded;
    stk::mesh::impl::DeletedElementInfoVector elementsDeleted;
    stk::mesh::EntityVector elementsChangedConn;
    size_t maxNumAdded = 0;
    size_t maxNumChanged = 0;
    bool changeEntityOwnerInProgress;
};

}} // end stk mesh namespaces

#endif
