#ifndef STK_REMOTE_SELECTOR_UPDATER_HPP
#define STK_REMOTE_SELECTOR_UPDATER_HPP

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequire
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

class RemoteSelectorUpdater : public stk::mesh::ModificationObserver
{
public:
    RemoteSelectorUpdater(stk::mesh::BulkData &bulk, stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector_,
                          stk::mesh::Selector sel)
    : bulkData(bulk), remoteActiveSelector(remoteActiveSelector_), selector(sel)
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
