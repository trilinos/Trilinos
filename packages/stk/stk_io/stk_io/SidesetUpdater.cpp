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
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/SidesetUpdater.hpp"


#if defined(UPDATER_DEBUG)
#define dout sierra::Env::outputP0()
#else
#define dout                                                                                                           \
  if(1)                                                                                                                \
    ;                                                                                                                  \
  else                                                                                                                 \
  sierra::Env::outputP0()
#endif

namespace stk { namespace io {

void SidesetUpdater::fill_sidesets_element_belongs_to(stk::mesh::Entity elem)
{
    std::vector<stk::mesh::SideSet *> sidesets = bulkData.get_sidesets();
    for (stk::mesh::SideSet* sideset : sidesets)
    {
        stk::mesh::Part* part = bulkData.mesh_meta_data().get_part(sideset->get_name());
        for (const stk::mesh::SideSetEntry& entry : *sideset)
        {
            if (entry.element == elem)
            {
                stkSideSets.insert(part);
                break;
            }
        }
    }
}

void SidesetUpdater::entity_added(stk::mesh::Entity entity)
{

}

void SidesetUpdater::entity_deleted(stk::mesh::Entity entity)
{
//    if(isActive) {
//        if ((bulkData.entity_rank(entity) == stk::topology::ELEMENT_RANK) && bulkData.bucket(entity).owned())
//        {
////            sierra::Env::outputP0() << "P[" << bulkData.parallel_rank() << "] Deleting element " << bulkData.identifier(entity) << std::endl;
//            fill_sidesets_element_belongs_to(entity);
//        }
//        if(bulkData.entity_rank(entity) == bulkData.mesh_meta_data().side_rank())
//        {
//            std::vector<const stk::mesh::Part*> sidesetParts = get_sideset_io_parts(bulkData, entity);
//            insert_parts(entity, sidesetParts, DELETED);
//        }
//    }
}

bool is_part_a_sideset(const stk::mesh::BulkData& bulkData, const stk::mesh::Part& part)
{
    bool isSideset = false;

    if(part.primary_entity_rank() == bulkData.mesh_meta_data().side_rank()) {
        const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(part);
        isSideset = bulkData.does_sideset_exist(parentPart);
    }

    return isSideset;
}

stk::mesh::Entity get_side(const stk::mesh::BulkData& bulk, const stk::mesh::SideSetEntry& entry)
{
    stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
    const stk::mesh::Entity* sides = bulk.begin(entry.element, sideRank);
    const stk::mesh::ConnectivityOrdinal* ordinals = bulk.begin_ordinals(entry.element, sideRank);
    unsigned numSides = bulk.num_connectivity(entry.element, sideRank);
    for(unsigned i=0; i<numSides; ++i)
    {
        if (ordinals[i] == entry.side)
        {
            return sides[i];
        }
    }
    return stk::mesh::Entity();
}

bool is_side_internal(const stk::mesh::BulkData& bulk, stk::mesh::Entity side,
                        const std::vector<const stk::mesh::Part*> &blocks, const stk::mesh::Selector& activeSelector)
{
    unsigned numElems = bulk.num_elements(side);
    if (numElems > 1)
    {
        const stk::mesh::Entity* elems = bulk.begin_elements(side);
        for(const stk::mesh::Part* block : blocks)
        {
            unsigned numSolidElementsInBlock = 0;

            for(unsigned i=0; i<numElems; ++i)
            {
                stk::mesh::Entity elem = elems[i];
                const stk::mesh::Bucket& bucket = bulk.bucket(elem);
                bool isShell = bucket.topology().is_shell();
                bool isInBlock = bucket.member(*block);
                bool isActive = activeSelector(bucket);
                if (!isShell && isInBlock && isActive)
                {
                    numSolidElementsInBlock++;
                }
            }

            if (numSolidElementsInBlock > 1)
            {
//                for(unsigned i=0; i<numElems; ++i)
//                {
//                    sierra::Env::outputP0()<<"for block "<<block->name()<<", elem "<<i<<": "<<bulk.identifier(elems[i])<<", topo: "<<bulk.bucket(elems[i]).topology()
//                    <<", parts: "<<std::endl;
//                    const stk::mesh::PartVector& parts = bulk.bucket(elems[i]).supersets();
//                    for(const stk::mesh::Part* part : parts) {
//                      sierra::Env::outputP0()<<"      "<<part->name()<<std::endl;
//                    }
//                }
                return true;
            }
        }
    }

    return false;
}

bool is_sideset_internal(const stk::mesh::BulkData& bulk, const stk::mesh::Part& surfacePart, const stk::mesh::Selector& activeSelector)
{
    if (is_part_a_sideset(bulk, surfacePart))
    {
        std::vector<const stk::mesh::Part*> blocks = bulk.mesh_meta_data().get_blocks_touching_surface(&surfacePart);

        stk::mesh::EntityVector sides;
        stk::mesh::get_selected_entities(surfacePart, bulk.buckets(bulk.mesh_meta_data().side_rank()), sides);

        for(stk::mesh::Entity side : sides) {
            if(is_side_internal(bulk, side, blocks, activeSelector)) {
                return true;
            }
        }
    }
    return false;
}

void issue_internal_sideset_warning(const std::string& sidesetName)
{
  std::cerr<<"WARNING, Internal sideset ("<<sidesetName<<") detected. STK doesn't support internal sidesets\n"
           <<"(i.e., sidesets between elements where both elements are in the same block)\n"
           <<"Execution will continue but correct results are not guaranteed. Contact sierra-help@sandia.gov"<<std::endl;
}

void SidesetUpdater::reconstruct_noninternal_sidesets(const std::vector<size_t> &reducedValues)
{
    if(bulkData.was_mesh_modified_since_sideset_creation())
    {
        std::vector<const stk::mesh::Part *> surfacesInMap = bulkData.mesh_meta_data().get_surfaces_in_surface_to_block_map();
        std::set<const stk::mesh::Part*, part_compare_by_ordinal> parents;
        ThrowRequireMsg(reducedValues.size() == surfacesInMap.size(), "ReducedValues wrong size!");

        for(unsigned i=0; i<surfacesInMap.size(); ++i) {
            const stk::mesh::Part* part = surfacesInMap[i];
            bool isInternal = reducedValues[i] == 1;
            if (!isInternal)
            {
                const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(*part);
                parents.insert(&parentPart);
            } else {
//                sierra::Env::outputP0() << "P[" << bulkData.parallel_rank() << "] Skipping internal sideset: " << part->name() << std::endl;
                if (!internalSidesetWarningHasBeenIssued && bulkData.parallel_rank()==0)
                {
                    issue_internal_sideset_warning(part->name());
                }
                internalSidesetWarningHasBeenIssued = true;
            }
        }

        for(auto part : parents) {
            bulkData.clear_sideset(*part);
        }

        for(size_t i = 0; i < surfacesInMap.size(); ++i)
        {
            const stk::mesh::Part* surfacePart = surfacesInMap[i];
            bool isInternal = reducedValues[i] == 1;
            if (!isInternal)
            {
//                sierra::Env::outputP0() << "P[" << bulkData.parallel_rank() << "] Constructing sideset: " << surfacePart->name() << std::endl;
                std::vector<const stk::mesh::Part *> touching_parts = bulkData.mesh_meta_data().get_blocks_touching_surface(surfacePart);

                stk::mesh::Selector elementSelector = stk::mesh::selectUnion(touching_parts);
                fill_sideset(*surfacePart, bulkData, elementSelector);
            }
        }
    }
}

void SidesetUpdater::reconstruct_sidesets()
{
    std::set<const stk::mesh::Part*, part_compare_by_ordinal> parents;
    std::set<const stk::mesh::Part*, part_compare_by_ordinal> subsets;

    for(auto part : stkSideSets) {
        const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(*part);
        parents.insert(&parentPart);

        if(part->subsets().empty()) {
            subsets.insert(part);
        }
    }

    for(auto part : parents) {
        bulkData.clear_sideset(*part);
    }

    for(auto part : subsets) {
        std::vector<const stk::mesh::Part *> touching_parts = bulkData.mesh_meta_data().get_blocks_touching_surface(part);

        stk::mesh::Selector elementSelector = stk::mesh::selectUnion(touching_parts);
        fill_sideset(*part, bulkData, elementSelector);
    }
}

void SidesetUpdater::finished_modification_end_notification()
{
    stkSideSets.clear();
}

void SidesetUpdater::started_modification_end_notification()
{
}

void SidesetUpdater::fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
{
    valuesToReduce.clear();
    if(bulkData.was_mesh_modified_since_sideset_creation())
    {
        std::vector<const stk::mesh::Part *> surfacesInMap = bulkData.mesh_meta_data().get_surfaces_in_surface_to_block_map();
        std::set<const stk::mesh::Part*, part_compare_by_ordinal> parents;

        valuesToReduce.assign(surfacesInMap.size(), 0);

        for(unsigned i=0; i<surfacesInMap.size(); ++i) {
            const stk::mesh::Part* part = surfacesInMap[i];
            bool isInternal = is_sideset_internal(bulkData, *part, activeSelector);
            if (isInternal)
            {
                valuesToReduce[i] = 1;
            }
        }
    }
}

void SidesetUpdater::set_reduced_values(const std::vector<size_t> &reducedValues)
{
    if(bulkData.was_mesh_modified_since_sideset_creation() && isActive)
    {
        reconstruct_noninternal_sidesets(reducedValues);
    }
}

void SidesetUpdater::elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
{
}

void SidesetUpdater::elements_moved_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
{
}


const std::string &op_string(OPERATION op)
{
    static std::string added("added");
    static std::string destroyed("destroyed");
    static std::string unknown("unknown");

    if(op == ADDED) return added;
    if(op == DELETED) return destroyed;

    return unknown;
}
void SidesetUpdater::tag_sideset(stk::mesh::Entity entity, const stk::mesh::Part& part, OPERATION op)
{
    if (bulkData.bucket_ptr(entity) != nullptr)
    {
        bool ownedOrShared = bulkData.bucket(entity).owned() || bulkData.bucket(entity).shared();
        if(is_part_a_sideset(bulkData, part) && ownedOrShared) {
            stkSideSets.insert(&part);
        }
    }
}

void SidesetUpdater::insert_parts(stk::mesh::Entity entity, const stk::mesh::ConstPartVector& parts, OPERATION op)
{
    for(const stk::mesh::Part* part : parts) {
        tag_sideset(entity, *part, op);
    }
}

void SidesetUpdater::insert_parts(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts, OPERATION op)
{
    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    for(stk::mesh::Ordinal partOrdinal : parts) {
        const stk::mesh::Part& part = *meta.get_parts()[partOrdinal];
        tag_sideset(entity, part, op);
    }
}

void SidesetUpdater::entity_parts_added(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts)
{
//    if(isActive) {
//        insert_parts(entity, parts, ADDED);
//    }
}

void SidesetUpdater::entity_parts_removed(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts)
{
//    if(isActive) {
//        insert_parts(entity, parts, DELETED);
//    }
}

void SidesetUpdater::set_active(bool active)
{
    isActive = active;
}


void toggle_sideset_updaters(stk::mesh::BulkData& bulk, bool flag)
{
    std::vector<SidesetUpdater*> updaters = bulk.get_observer_type<SidesetUpdater>();
    for(SidesetUpdater* updater : updaters) {
        updater->set_active(flag);
    }
}

} } //namespaces
