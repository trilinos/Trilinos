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
  if(isActive) {
    remove_element_entries_from_sidesets(bulkData, entity, &sidesetPartsWithDeletedEntries);
  }
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

std::pair<bool,bool> is_side_internal_and_modified(const stk::mesh::BulkData& bulk, stk::mesh::Entity side,
                        const std::vector<const stk::mesh::Part*> &blocks, const stk::mesh::Selector& activeSelector)
{
    bool isModified = !(bulk.state(side) == stk::mesh::Unchanged);
    bool isInternal = false;
    unsigned numElems = bulk.num_elements(side);
    if (numElems > 1) {
        const stk::mesh::Entity* elems = bulk.begin_elements(side);
        for (const stk::mesh::Part* block : blocks) {
            unsigned numSolidElementsInBlock = 0;

            for (unsigned i=0; i<numElems; ++i) {
                stk::mesh::Entity elem = elems[i];
                const stk::mesh::Bucket& bucket = bulk.bucket(elem);
                if (!bucket.topology().is_shell()) {
                    if (bucket.member(block->mesh_meta_data_ordinal())) {
                        if (activeSelector(bucket)) {
                            numSolidElementsInBlock++;
                        }
                    }
                }
                if (!isModified) {
                    isModified = !(bulk.state(elem) == stk::mesh::Unchanged);
                }
            }

            isInternal |= (numSolidElementsInBlock > 1);
        }
    }

    return std::make_pair(isInternal, isModified);
}

std::pair<bool,bool> is_sideset_internal_and_modified(const stk::mesh::BulkData& bulk, const stk::mesh::Part& surfacePart, const stk::mesh::Selector& activeSelector)
{
    std::pair<bool,bool> isInternalAndModified(false,false);
    if (is_part_a_sideset(bulk, surfacePart))
    {
        std::vector<const stk::mesh::Part*> blocks = bulk.mesh_meta_data().get_blocks_touching_surface(&surfacePart);

        const stk::mesh::BucketVector& buckets = bulk.get_buckets(bulk.mesh_meta_data().side_rank(), surfacePart);
        for(const stk::mesh::Bucket* bptr : buckets) {
            for(stk::mesh::Entity side : *bptr) {
                std::pair<bool,bool> internalAndModified = is_side_internal_and_modified(bulk, side, blocks, activeSelector);
                isInternalAndModified.first |= internalAndModified.first;
                isInternalAndModified.second |= internalAndModified.second;
                if (isInternalAndModified.first && isInternalAndModified.second) {
                  return isInternalAndModified;
                }
            }
        }
    }
    return isInternalAndModified;
}

void issue_internal_sideset_warning(const std::string& sidesetName, std::ostream& ostrm)
{
  ostrm <<"WARNING, Internal sideset ("<<sidesetName<<") detected. STK doesn't support internal sidesets\n"
        <<"(i.e., sidesets between elements where both elements are in the same block)\n"
        <<"Execution will continue but correct results are not guaranteed. Contact sierra-help@sandia.gov"<<std::endl;
}

void SidesetUpdater::update_sidesets_without_surface_block_mapping(stk::mesh::BulkData &bulk)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  std::vector<const stk::mesh::Part *> surfacesInMap = meta.get_surfaces_in_surface_to_block_map();
  std::set<const stk::mesh::Part*> alreadyUpdatedSidesetParts;
  std::set<const stk::mesh::Part*> difference;

  for(const stk::mesh::Part *part : surfacesInMap) {
    if(part->subsets().empty()) {
        const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(*part);

        bool sidesetExists = bulk.does_sideset_exist(parentPart);
        if(sidesetExists) {
          alreadyUpdatedSidesetParts.insert( &parentPart );
        }
    }
  }

  std::set_difference(sidesetPartsWithDeletedEntries.begin(), sidesetPartsWithDeletedEntries.end(),
                      alreadyUpdatedSidesetParts.begin(), alreadyUpdatedSidesetParts.end(), std::inserter(difference, difference.end()));

  for (const stk::mesh::Part* sidePart: difference)
  {
    bulk.clear_sideset(*sidePart);
    fill_sideset(*sidePart, bulk, meta.universal_part());
  }
}

bool value_is_from_internal_sideset(size_t value)
{
  return value > 0;
}

void SidesetUpdater::reconstruct_noninternal_sidesets(const std::vector<size_t> &reducedValues)
{
    bool reconstructed = false;
    if(bulkData.was_mesh_modified_since_sideset_creation())
    {
        std::vector<const stk::mesh::Part *> surfacesInMap = bulkData.mesh_meta_data().get_surfaces_in_surface_to_block_map();
        std::set<const stk::mesh::Part*, part_compare_by_ordinal> parents;
        ThrowRequireMsg(reducedValues.size() == surfacesInMap.size(), "ReducedValues wrong size!");

        for(unsigned i=0; i<surfacesInMap.size(); ++i) {
            const stk::mesh::Part* part = surfacesInMap[i];
            const bool isInternal = value_is_from_internal_sideset(reducedValues[i]);
            const bool isInternalAndModified = (reducedValues[i] == 2);
            if (!isInternal)
            {
                const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(*part);
                parents.insert(&parentPart);
            } else {
                if (!internalSidesetWarningHasBeenIssued &&
                    isInternalAndModified &&
                    bulkData.parallel_rank()==0 &&
                    outputStream != nullptr)
                {
                    issue_internal_sideset_warning(part->name(), *outputStream);
                }
                internalSidesetWarningHasBeenIssued = true;
            }
        }

        update_sidesets_without_surface_block_mapping(bulkData);

        for(auto part : parents) {
            bulkData.clear_sideset(*part);
        }

        for(size_t i = 0; i < surfacesInMap.size(); ++i)
        {
            const stk::mesh::Part* surfacePart = surfacesInMap[i];
            bool isInternal = value_is_from_internal_sideset(reducedValues[i]);
            if (!isInternal)
            {
                std::vector<const stk::mesh::Part *> touching_parts = bulkData.mesh_meta_data().get_blocks_touching_surface(surfacePart);

                stk::mesh::Selector elementSelector = stk::mesh::selectUnion(touching_parts);
                if(touching_parts.size() == 0) {
                  elementSelector = bulkData.mesh_meta_data().universal_part();
                }
                fill_sideset(*surfacePart, bulkData, elementSelector);
                reconstructed = true;
            }
        }

        if (reconstructed)
        {
            bulkData.synchronize_sideset_sync_count();
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
  sidesetPartsWithDeletedEntries.clear();
}

void SidesetUpdater::fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
{
    valuesToReduce.clear();
    if(bulkData.was_mesh_modified_since_sideset_creation())
    {
        std::vector<const stk::mesh::Part *> surfacesInMap = bulkData.mesh_meta_data().get_surfaces_in_surface_to_block_map();

        valuesToReduce.assign(surfacesInMap.size(), 0);

        for(unsigned i=0; i<surfacesInMap.size(); ++i) {
            const stk::mesh::Part* part = surfacesInMap[i];
            std::pair<bool,bool> isInternalAndModified = is_sideset_internal_and_modified(bulkData, *part, activeSelector);
            const bool isInternalSideset = isInternalAndModified.first;
            const bool isModified = isInternalAndModified.second;
            if (isInternalSideset)
            {
                if (isModified) {
                    valuesToReduce[i] = 2;
                }
                else {
                    valuesToReduce[i] = 1;
                }
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
    std::vector<std::shared_ptr<SidesetUpdater>> updaters = bulk.get_observer_type<SidesetUpdater>();
    bool changedFlag = false;
    for(std::shared_ptr<SidesetUpdater> updater : updaters) {
        if (updater->get_active_flag() != flag) {
            changedFlag = true;
            updater->set_active(flag);
        }
    }

    if (flag == true && changedFlag == true) {
        //We need to synchronize this so that the sideset updater won't try
        //to update anything for modifications that were done while the
        //updater was inactive.
        bulk.synchronize_sideset_sync_count();
    }
}

} } //namespaces
