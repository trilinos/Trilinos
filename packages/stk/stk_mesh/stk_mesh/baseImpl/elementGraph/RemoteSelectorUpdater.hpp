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

#ifndef STK_REMOTE_SELECTOR_UPDATER_HPP
#define STK_REMOTE_SELECTOR_UPDATER_HPP

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequire
#include <stk_util/parallel/ParallelReduce.hpp>
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector

#include "ElemElemGraph.hpp"// for ElemElemGraph
#include "ElemElemGraphImpl.hpp"

namespace stk { namespace mesh {
class BulkData;

class RemoteSelectorUpdater : public stk::mesh::ModificationObserver
{
public:
    RemoteSelectorUpdater(stk::mesh::BulkData &bulk, stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector_,
                          stk::mesh::Selector sel)
    : stk::mesh::ModificationObserver(stk::mesh::ModificationObserverPriority::APPLICATION)
    , bulkData(bulk), remoteActiveSelector(remoteActiveSelector_), selector(sel)
    {
    }

    virtual void entity_added(stk::mesh::Entity entity)
    {
        if(bulkData.entity_rank(entity) == stk::topology::ELEM_RANK)
        {
            numModifiedElems++;
        }
    }

    virtual void entity_deleted(stk::mesh::Entity entity)
    {
        if(bulkData.entity_rank(entity) == stk::topology::ELEM_RANK && bulkData.bucket(entity).owned())
        {
            numModifiedElems++;
        }
    }

    virtual void finished_modification_end_notification()
    {
        if (numModifiedElems > 0) {
            stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData,
                                                                         bulkData.get_face_adjacent_element_graph(),
                                                                         selector,
                                                                         remoteActiveSelector);
            numModifiedElems = 0;
        }
    }

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
    {
        valuesToReduce.clear();
        valuesToReduce.push_back(numModifiedElems);
    }

    virtual void set_reduced_values(const std::vector<size_t> &reducedValues)
    {
        numModifiedElems = reducedValues[0];
    }

    virtual void elements_moved_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData,
                                                                     bulkData.get_face_adjacent_element_graph(),
                                                                     selector,
                                                                     remoteActiveSelector);
        numModifiedElems = 0;
    }
private:
    stk::mesh::BulkData &bulkData;
    stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector;
    stk::mesh::Selector selector;
    size_t numModifiedElems = 0;
};

}} // end stk mesh namespaces

#endif
