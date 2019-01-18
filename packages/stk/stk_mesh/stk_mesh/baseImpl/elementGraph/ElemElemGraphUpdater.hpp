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

#ifndef STK_ELEM_ELEM_GRAPH_UPDATER_HPP
#define STK_ELEM_ELEM_GRAPH_UPDATER_HPP

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequire
#include <stk_util/parallel/ParallelReduce.hpp>
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector

#include "ElemElemGraph.hpp"// for ElemElemGraph
#include "ElemElemGraphImpl.hpp"

namespace stk { namespace mesh {

class ElemElemGraphUpdater : public stk::mesh::ModificationObserver
{
public:
    ElemElemGraphUpdater(stk::mesh::BulkData &bulk, stk::mesh::ElemElemGraph &elemGraph_)
    : bulkData(bulk), elemGraph(elemGraph_), changeEntityOwnerInProgress(false)
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
            elementsDeleted.push_back({entity, bulkData.identifier(entity), bulkData.bucket(entity).topology().is_shell()});
        }
    }

    virtual void finished_modification_end_notification()
    {
        if (maxNumAdded > 0) {
            elemGraph.add_elements(elementsAdded);
            elementsAdded.clear();
        }
    }

    virtual void started_modification_end_notification()
    {
        if (get_global_sum(bulkData.parallel(), elementsDeleted.size()) > 0) {
            elemGraph.delete_elements(elementsDeleted);
            elementsDeleted.clear();
        }
    }

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
    {
        valuesToReduce.clear();
        valuesToReduce.push_back(elementsAdded.size());
    }

    virtual void set_reduced_values(const std::vector<size_t> &reducedValues)
    {
        maxNumAdded = reducedValues[0];
    }

    virtual void elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        changeEntityOwnerInProgress = true;
    }

    virtual void elements_moved_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        elemGraph.fill_from_mesh();
        changeEntityOwnerInProgress = false;
    }
private:
    stk::mesh::BulkData &bulkData;
    stk::mesh::ElemElemGraph &elemGraph;
    stk::mesh::impl::ParallelGraphInfo newParallelGraphEntries;
    stk::mesh::EntityVector elementsAdded;
    stk::mesh::impl::DeletedElementInfoVector elementsDeleted;
    size_t maxNumAdded = 0;
    bool changeEntityOwnerInProgress;
};

}} // end stk mesh namespaces

#endif
