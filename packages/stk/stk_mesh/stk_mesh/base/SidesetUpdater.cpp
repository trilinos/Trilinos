#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequire
#include <stk_util/environment/Env.hpp>
#include <algorithm>
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector
#include "stk_mesh/base/SideSetUtil.hpp"
#include "stk_mesh/base/SidesetUpdater.hpp"


#if defined(UPDATER_DEBUG)
#define dout sierra::Env::outputP0()
#else
#define dout                                                                                                           \
  if(1)                                                                                                                \
    ;                                                                                                                  \
  else                                                                                                                 \
  sierra::Env::outputP0()
#endif

namespace stk { namespace mesh {

void ReconstructionSidesetUpdater::fill_sidesets_element_belongs_to(Entity elem)
{
  std::vector<SideSet *> sidesets = m_bulkData.get_sidesets();
  for (SideSet* sideset : sidesets)
  {
    Part* part = m_metaData.get_part(sideset->get_name());
    for (const SideSetEntry& entry : *sideset)
    {
      if (entry.element == elem)
      {
        m_stkSideSets.insert(part);
        break;
      }
    }
  }
}

void ReconstructionSidesetUpdater::entity_added(Entity entity)
{

}

void ReconstructionSidesetUpdater::entity_deleted(Entity entity)
{
  if(m_isActive) {
    m_helper.remove_element_entries_from_sidesets(entity, &m_sidesetPartsWithDeletedEntries);
  }
}

bool is_part_a_sideset(const BulkData& bulkData, const Part& part)
{
  bool isSideset = false;

  if(part.primary_entity_rank() == bulkData.mesh_meta_data().side_rank()) {
    const Part &parentPart = get_sideset_parent(part);
    isSideset = bulkData.does_sideset_exist(parentPart);
  }

  return isSideset;
}

Entity get_side(const BulkData& bulk, const SideSetEntry& entry)
{
  EntityRank sideRank = bulk.mesh_meta_data().side_rank();
  const Entity* sides = bulk.begin(entry.element, sideRank);
  const ConnectivityOrdinal* ordinals = bulk.begin_ordinals(entry.element, sideRank);
  unsigned numSides = bulk.num_connectivity(entry.element, sideRank);
  for(unsigned i=0; i<numSides; ++i)
  {
    if (ordinals[i] == entry.side)
    {
      return sides[i];
    }
  }
  return Entity();
}

std::pair<bool,bool> is_side_internal_and_modified(const BulkData& bulk, Entity side,
    const std::vector<const Part*> &blocks,
    const Selector& activeSelector)
{
  bool isModified = !(bulk.state(side) == Unchanged);
  bool isInternal = false;
  unsigned numElems = bulk.num_elements(side);
  if (numElems > 1) {
    const Entity* elems = bulk.begin_elements(side);
    for (const Part* block : blocks) {
      unsigned numSolidElementsInBlock = 0;

      for (unsigned i=0; i<numElems; ++i) {
        Entity elem = elems[i];
        const Bucket& bucket = bulk.bucket(elem);
        if (!bucket.topology().is_shell()) {
          if (bucket.member(block->mesh_meta_data_ordinal())) {
            if (activeSelector(bucket)) {
              numSolidElementsInBlock++;
            }
          }
        }
        if (!isModified) {
          isModified = !(bulk.state(elem) == Unchanged);
        }
      }

      isInternal |= (numSolidElementsInBlock > 1);
    }
  }

  return std::make_pair(isInternal, isModified);
}

std::pair<bool,bool> is_sideset_internal_and_modified(const BulkData& bulk, const Part& surfacePart,
    const Selector& activeSelector)
{
  std::pair<bool,bool> isInternalAndModified(false,false);
  if (is_part_a_sideset(bulk, surfacePart))
  {
    std::vector<const Part*> blocks = bulk.mesh_meta_data().get_blocks_touching_surface(&surfacePart);

    const BucketVector& buckets = bulk.get_buckets(bulk.mesh_meta_data().side_rank(), surfacePart);
    for(const Bucket* bptr : buckets) {
      for(Entity side : *bptr) {
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

void ReconstructionSidesetUpdater::update_sidesets_without_surface_block_mapping()
{
  std::vector<const Part *> surfacesInMap = m_metaData.get_surfaces_in_surface_to_block_map();
  std::set<const Part*> alreadyUpdatedSidesetParts;
  std::set<const Part*> difference;

  for(const Part *part : surfacesInMap) {
    if(part->subsets().empty()) {
      const Part &parentPart = get_sideset_parent(*part);

      bool sidesetExists = m_bulkData.does_sideset_exist(parentPart);
      if(sidesetExists) {
        alreadyUpdatedSidesetParts.insert( &parentPart );
      }
    }
  }

  std::set_difference(m_sidesetPartsWithDeletedEntries.begin(), m_sidesetPartsWithDeletedEntries.end(),
      alreadyUpdatedSidesetParts.begin(), alreadyUpdatedSidesetParts.end(), std::inserter(difference, difference.end()));

  for (const Part* sidePart: difference)
  {
    m_bulkData.clear_sideset(*sidePart);
    fill_sideset(*sidePart, m_bulkData, m_metaData.universal_part());
  }
}

bool value_is_from_internal_sideset(size_t value)
{
  return value > 0;
}

void ReconstructionSidesetUpdater::reconstruct_noninternal_sidesets(const std::vector<size_t> &reducedValues)
{
  bool reconstructed = false;
  if(m_bulkData.was_mesh_modified_since_sideset_creation())
  {
    std::vector<const Part *> surfacesInMap = m_metaData.get_surfaces_in_surface_to_block_map();
    std::set<const Part*, part_compare_by_ordinal> parents;
    STK_ThrowRequireMsg(reducedValues.size() == surfacesInMap.size(), "ReducedValues wrong size!");

    for(unsigned i=0; i<surfacesInMap.size(); ++i) {
      const Part* part = surfacesInMap[i];
      const bool isInternal = value_is_from_internal_sideset(reducedValues[i]);
      const bool isInternalAndModified = (reducedValues[i] == 2);
      if (!isInternal)
      {
        const Part &parentPart = get_sideset_parent(*part);
        parents.insert(&parentPart);
      } else {
        if (!m_internalSidesetWarningHasBeenIssued &&
            isInternalAndModified &&
            m_bulkData.parallel_rank()==0 &&
            m_outputStream != nullptr &&
            m_warnAboutInternalSideset)
        {
          issue_internal_sideset_warning(part->name(), *m_outputStream);
        }
        m_internalSidesetWarningHasBeenIssued = true;
      }
    }

    update_sidesets_without_surface_block_mapping();

    for(auto part : parents) {
      m_bulkData.clear_sideset(*part);
    }

    for(size_t i = 0; i < surfacesInMap.size(); ++i)
    {
      const Part* surfacePart = surfacesInMap[i];
      bool isInternal = value_is_from_internal_sideset(reducedValues[i]);

      if (!isInternal)
      {
        std::vector<const Part *> touching_parts = m_metaData.get_blocks_touching_surface(surfacePart);
        Selector elementSelector = selectUnion(touching_parts);
        if(touching_parts.size() == 0) {
          elementSelector = m_metaData.universal_part();
        }
        fill_sideset(*surfacePart, m_bulkData, elementSelector);
        reconstructed = true;
      }
    }

    if (reconstructed)
    {
      m_bulkData.synchronize_sideset_sync_count();
    }
  }
}

void ReconstructionSidesetUpdater::reconstruct_sidesets()
{
  std::set<const Part*, part_compare_by_ordinal> parents;
  std::set<const Part*, part_compare_by_ordinal> subsets;

  for(auto part : m_stkSideSets) {
    const Part &parentPart = get_sideset_parent(*part);
    parents.insert(&parentPart);

    if(part->subsets().empty()) {
      subsets.insert(part);
    }
  }

  for(auto part : parents) {
    m_bulkData.clear_sideset(*part);
  }

  for(auto part : subsets) {
    std::vector<const Part *> touching_parts = m_metaData.get_blocks_touching_surface(part);

    Selector elementSelector = selectUnion(touching_parts);
    fill_sideset(*part, m_bulkData, elementSelector);
  }
}

void ReconstructionSidesetUpdater::modification_begin_notification()
{

}

void ReconstructionSidesetUpdater::finished_modification_end_notification()
{
  m_stkSideSets.clear();
  m_helper.reset_internal_sideset_detection(false);
}

void ReconstructionSidesetUpdater::started_modification_end_notification()
{
  m_sidesetPartsWithDeletedEntries.clear();
}

void ReconstructionSidesetUpdater::fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
{
  if (!m_isActive){return;}
  valuesToReduce.clear();
  if(m_bulkData.was_mesh_modified_since_sideset_creation())
  {
    std::vector<const Part *> surfacesInMap = m_metaData.get_surfaces_in_surface_to_block_map();

    valuesToReduce.assign(surfacesInMap.size(), 0);

    for(unsigned i=0; i<surfacesInMap.size(); ++i) {
      const Part* part = surfacesInMap[i];
      std::pair<bool,bool> isInternalAndModified = is_sideset_internal_and_modified(m_bulkData, *part, m_activeSelector);
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

void ReconstructionSidesetUpdater::set_reduced_values(const std::vector<size_t> &reducedValues)
{
  if(m_bulkData.was_mesh_modified_since_sideset_creation() && m_isActive)
  {
    reconstruct_noninternal_sidesets(reducedValues);
  }
}

void ReconstructionSidesetUpdater::elements_about_to_move_procs_notification(const EntityProcVec &elemProcPairsToMove)
{
}

void ReconstructionSidesetUpdater::elements_moved_procs_notification(const EntityProcVec &elemProcPairsToMove)
{
}

void ReconstructionSidesetUpdater::tag_sideset(Entity entity, const Part& part)
{
  if (m_bulkData.bucket_ptr(entity) != nullptr)
  {
    bool ownedOrShared = m_bulkData.bucket(entity).owned() || m_bulkData.bucket(entity).shared();
    if(is_part_a_sideset(m_bulkData, part) && ownedOrShared) {
      m_stkSideSets.insert(&part);
    }
  }
}

void ReconstructionSidesetUpdater::insert_parts(Entity entity, const ConstPartVector& parts)
{
  for(const Part* part : parts) {
    tag_sideset(entity, *part);
  }
}

void ReconstructionSidesetUpdater::insert_parts(Entity entity, const OrdinalVector& parts)
{
  const MetaData& meta = m_bulkData.mesh_meta_data();
  for(Ordinal partOrdinal : parts) {
    const Part& part = *meta.get_parts()[partOrdinal];
    tag_sideset(entity, part);
  }
}

void ReconstructionSidesetUpdater::entity_parts_added(Entity entity, const OrdinalVector& parts)
{
}

void ReconstructionSidesetUpdater::entity_parts_removed(Entity entity, const OrdinalVector& parts)
{
}

const Selector& IncrementalSidesetUpdater::get_cached_sideset_selector(const Part* part)
{
  std::map<const Part*, Selector>::iterator iter = m_cachedSidesetSelectors.find(part);

  if(iter == m_cachedSidesetSelectors.end()) {
    std::vector<const Part *> touchingBlocks = m_metaData.get_blocks_touching_surface(part);
    Selector selector = selectUnion(touchingBlocks);
    const auto entry = m_cachedSidesetSelectors.insert({part, selector});
    return entry.first->second;
  }

  return iter->second;
}

void IncrementalSidesetUpdater::fill_surfaces_touching_block(const Part* inputBlock, ConstPartVector& surfacesTouchingBlock)
{
  surfacesTouchingBlock.clear();

  std::vector<const Part *> allSurfaces = m_metaData.get_surfaces_in_surface_to_block_map();

  for (const Part * surface : allSurfaces) {
    std::vector<const Part*> blocks = m_metaData.get_blocks_touching_surface(surface);
    std::vector<const Part*>::iterator lowerBound = std::lower_bound(blocks.begin(), blocks.end(), inputBlock, PartLess());

    if (lowerBound != blocks.end() && *lowerBound == inputBlock) {
      surfacesTouchingBlock.push_back(surface);
    }
  }
}

const ConstPartVector& IncrementalSidesetUpdater::get_cached_surfaces_touching_block(const Part& part)
{
  BlockToSurfaceMapping blockToSurface(&part);

  std::vector<BlockToSurfaceMapping>::iterator lowerBound = std::lower_bound(m_cachedBlockToSurfaceMapping.begin(), m_cachedBlockToSurfaceMapping.end(), blockToSurface);
  if(lowerBound == m_cachedBlockToSurfaceMapping.end() || *lowerBound != blockToSurface) {
    fill_surfaces_touching_block(&part, blockToSurface.surfaces);

    auto iter = m_cachedBlockToSurfaceMapping.insert(lowerBound, blockToSurface);
    return iter->surfaces;
  }

  return lowerBound->surfaces;
}

void IncrementalSidesetUpdater::accumulate_element_block_part_changes(Entity entity, const OrdinalVector& partOrdinals)
{
  if(m_bulkData.entity_rank(entity) != stk::topology::ELEM_RANK || !m_bulkData.is_valid(entity)) return;

  PartChangeAccumulatorVector::iterator lowerBound = std::lower_bound(m_accumulatedElementPartChanges.begin(),
                                                                      m_accumulatedElementPartChanges.end(), entity);

  if(lowerBound == m_accumulatedElementPartChanges.end() || *lowerBound != entity) {
    PartChangeAccumulator accumulator(entity);
    lowerBound = m_accumulatedElementPartChanges.insert(lowerBound, accumulator);
  }

  for(auto partOrdinal : partOrdinals) {
    Part& part = m_metaData.get_part(partOrdinal);
    if(part.primary_entity_rank() == stk::topology::ELEM_RANK && part.id() != Part::INVALID_ID) {
      stk::util::insert_keep_sorted_and_unique(partOrdinal, lowerBound->partOrdinals);
    }
  }
}

void IncrementalSidesetUpdater::fill_info_for_element_affected_by_block_part_change(Entity entity,
                                                                                   SideSetSelectorVector& sidesetsAndSelectors,
                                                                                   std::vector<std::pair<Entity, unsigned>>& sideAndPartOrdinalVec)
{
  sideAndPartOrdinalVec.clear();

  if(m_bulkData.bucket_ptr(entity) == nullptr) { return; }

  const unsigned numSides = m_bulkData.num_sides(entity);

  if(sidesetsAndSelectors.size() > 0 && m_bulkData.entity_rank(entity) == stk::topology::ELEM_RANK && numSides > 0) {
    const stk::mesh::Entity* sides = m_bulkData.begin(entity, m_metaData.side_rank());

    for(unsigned i=0; i<numSides; i++) {
     for(SideSetSelector& sidesetSelector : sidesetsAndSelectors) {
       const Part& sidesetPart = sidesetSelector.part();
       unsigned sidesetPartOrdinal = sidesetPart.mesh_meta_data_ordinal();

       if(m_bulkData.bucket(sides[i]).member(sidesetPart)) {
         sideAndPartOrdinalVec.push_back(std::make_pair(sides[i], sidesetPartOrdinal));
       }
     }
    }
  }
}

void IncrementalSidesetUpdater::resolve_faces_affected_by_connected_element_part_change()
{
  SideSetSelectorVector sidesetsAndSelectors;
  std::vector<std::pair<Entity, unsigned>> sideAndPartOrdinals;

  for(const PartChangeAccumulator& accumulator : m_accumulatedElementPartChanges) {
    Entity element = accumulator.entity;
    const OrdinalVector& parts = accumulator.partOrdinals;

    fill_sidesets_and_selectors_from_blocks_touching_surfaces(parts, sidesetsAndSelectors);
    fill_info_for_element_affected_by_block_part_change(element, sidesetsAndSelectors, sideAndPartOrdinals);

    for(const auto& sideAndPartOrdinal : sideAndPartOrdinals) {
      Entity side = sideAndPartOrdinal.first;
      unsigned surfaceOrdinal = sideAndPartOrdinal.second;

      unsigned numElements = m_bulkData.num_elements(side);
      const stk::mesh::Entity* elements = m_bulkData.begin_elements(side);
      const ConnectivityOrdinal* ordinals = m_bulkData.begin_element_ordinals(side);

      const stk::mesh::Part& part = m_metaData.get_part(surfaceOrdinal);

      STK_ThrowAssert(part.primary_entity_rank() == m_metaData.side_rank());

      const stk::mesh::Part &parentPart = get_sideset_parent(part);
      bool sidesetExists = m_bulkData.does_sideset_exist(parentPart);

      if(sidesetExists) {
        SideSet& sideset = m_bulkData.get_sideset(parentPart);

        const Selector& elementSelector  = get_cached_sideset_selector(&part);

        for(unsigned i=0;i<numElements;++i) {
          bool isOwned = m_bulkData.bucket(elements[i]).owned();
          bool isSelected = elementSelector(m_bulkData.bucket(elements[i]));

          if(!isOwned || !isSelected) {
            SideSetEntry entry(elements[i], ordinals[i]);
            m_helper.remove_side_entry_from_sideset(&sideset, entry);
          }
        }
      }
    }
  }
}

void IncrementalSidesetUpdater::resolve_relation_updates()
{
  SideSetSelectorVector sidesets;

  if (m_relationUpdatesRemoved) {
    m_relationUpdates.erase(std::remove_if(m_relationUpdates.begin(), m_relationUpdates.end(),
                           [](const RelationUpdate& relUpdate){return relUpdate.removed;}),
                           m_relationUpdates.end());
    m_relationUpdatesRemoved = false;
  }

  if (!m_relationUpdatesAreSorted) {
    stk::util::sort_and_unique(m_relationUpdates);
    m_relationUpdatesAreSorted = true;
  }

  for(const RelationUpdate& update : m_relationUpdates) {
    const PartVector& parts = m_bulkData.bucket(update.face).supersets();
    fill_sidesets_and_selectors_from_parts(parts, sidesets);

    for (SideSetSelector& sideset : sidesets) {
      if(sideset.sideset()->get_accept_all_internal_non_coincident_entries()) {
        m_helper.add_element_side_entry_to_sideset(sideset, update.element, update.ordinal);
      } else {
        m_helper.add_coincident_element_side_entry_to_sideset(sideset, update.face, update.element, update.ordinal);
      }
    }
  }
}

void IncrementalSidesetUpdater::started_modification_end_notification()
{

}

void IncrementalSidesetUpdater::modification_begin_notification()
{

}

void IncrementalSidesetUpdater::finished_modification_end_notification()
{
  resolve_relation_updates();
  resolve_faces_affected_by_connected_element_part_change();

  m_helper.warn_internal_sideset_detection();

  m_sidesets.clear();
  m_cachedSidesetSelectors.clear();
  m_relationUpdates.clear();
  m_surfacesTouchingBlocks.clear();
  m_cachedBlockToSurfaceMapping.clear();
  m_accumulatedElementPartChanges.clear();

  std::vector<SideSet*> sidesets = m_bulkData.get_sidesets();

  for(SideSet* sideset : sidesets) {
    sideset->clear_modification_flag();
  }

  m_helper.reset_internal_sideset_detection(m_elemOrSideChangedRankedParts);
  m_elemOrSideChangedRankedParts = false;
}

void IncrementalSidesetUpdater::remove_relation(Entity element, Entity face, ConnectivityOrdinal ordinal)
{
  RelationUpdate relation(face, element, ordinal);

  if (!m_relationUpdatesAreSorted) {
    stk::util::sort_and_unique(m_relationUpdates);
    m_relationUpdatesAreSorted = true;
  }

  auto iter = std::lower_bound(m_relationUpdates.begin(), m_relationUpdates.end(), relation);
  if(iter != m_relationUpdates.end() && *iter == relation) {
    iter->removed = true;
    m_relationUpdatesRemoved = true;
  }
}

void IncrementalSidesetUpdater::add_relation(Entity element, Entity face, ConnectivityOrdinal ordinal)
{
  m_relationUpdates.emplace_back(face, element, ordinal);
  m_relationUpdatesAreSorted = false;
}

void IncrementalSidesetUpdater::entity_deleted(Entity entity)
{
  update_sideset_vector();
  m_helper.remove_entity_entries_from_sidesets(m_sidesets, entity, nullptr);
}

void IncrementalSidesetUpdater::relation_destroyed(Entity from, Entity to, ConnectivityOrdinal ordinal)
{
  if(m_bulkData.entity_rank(from) == stk::topology::ELEM_RANK && m_bulkData.entity_rank(to) == m_metaData.side_rank()) {
    update_sideset_vector();
    m_helper.remove_side_entry_from_sidesets(m_sidesets, from, ordinal, nullptr);
    remove_relation(from, to, ordinal);
  }
}

void IncrementalSidesetUpdater::relation_declared(Entity from, Entity to, ConnectivityOrdinal ordinal)
{
  if(m_bulkData.entity_rank(from) == stk::topology::ELEM_RANK && m_bulkData.entity_rank(to) == m_metaData.side_rank()) {
    add_relation(from, to, ordinal);
  }
}

void IncrementalSidesetUpdater::update_sideset_vector()
{
  if(m_bulkData.get_number_of_sidesets() != m_sidesets.size()) {
    m_sidesets = m_bulkData.get_sidesets();
  }
}

void IncrementalSidesetUpdater::update_surfaces_touching_blocks_vector()
{
  if(m_metaData.count_surfaces_in_surface_to_block_map() != m_surfacesTouchingBlocks.size()) {
    m_surfacesTouchingBlocks = m_metaData.get_surfaces_in_surface_to_block_map();
    std::sort(m_surfacesTouchingBlocks.begin(), m_surfacesTouchingBlocks.end(), PartLess());
  }
}

SideSet* IncrementalSidesetUpdater::get_or_create_connectivity_based_sideset(const Part& part, const bool checkSubsets)
{
  SideSet* sideset = nullptr;
  bool filterBySubsets = checkSubsets ? part.subsets().empty() : true;

  if(part.primary_entity_rank() == m_metaData.side_rank() && filterBySubsets) {
    update_surfaces_touching_blocks_vector();

    auto iter = std::lower_bound(m_surfacesTouchingBlocks.begin(), m_surfacesTouchingBlocks.end(), &part, PartLess());
    if(iter != m_surfacesTouchingBlocks.end() && (*iter)->mesh_meta_data_ordinal() == part.mesh_meta_data_ordinal()) {
      const Part &parentPart = get_sideset_parent(part);
      sideset = &m_bulkData.create_sideset(parentPart);
    }
  }

  return sideset;
}

void IncrementalSidesetUpdater::add_sidesets_and_selectors_from_part(const Part& part, const bool checkSubsets, SideSetSelectorVector& sidesets)
{
  SideSet* sideset = get_or_create_connectivity_based_sideset(part, checkSubsets);
  if(nullptr == sideset) {
    return;
  }

  if(sideset->is_automatically_updated() == false) {
    return;
  }
  const Selector& selector = get_cached_sideset_selector(&part);
  stk::util::insert_keep_sorted_and_unique(SideSetSelector(part, sideset, &selector), sidesets);

  if(m_metaData.count_blocks_touching_surface(&part) > 1) {
    sideset->set_accept_all_internal_non_coincident_entries(true);

    for(const Part* subset : part.subsets()) {
      if(m_bulkData.does_sideset_exist(*subset)) {
        SideSet& subsetSideSet = m_bulkData.get_sideset(*subset);
        subsetSideSet.set_accept_all_internal_non_coincident_entries(true);
      }
    }
  }
}

void IncrementalSidesetUpdater::fill_sidesets_and_selectors_from_part_ordinals(const OrdinalVector& parts, SideSetSelectorVector& sidesets)
{
  constexpr bool checkSubsets = true;
  sidesets.clear();
  sidesets.reserve(m_sidesets.size());

  for(auto partOrdinal : parts) {
    Part& part = m_metaData.get_part(partOrdinal);
    add_sidesets_and_selectors_from_part(part, checkSubsets, sidesets);
  }
}

void IncrementalSidesetUpdater::fill_sidesets_and_selectors_from_parts(const PartVector& parts, SideSetSelectorVector& sidesets)
{
  constexpr bool checkSubsets = false;
  sidesets.clear();
  sidesets.reserve(m_sidesets.size());

  for(auto part : parts) {
    add_sidesets_and_selectors_from_part(*part, checkSubsets, sidesets);
  }
}

void IncrementalSidesetUpdater::fill_sidesets_and_selectors(SideSetSelectorVector& sidesets)
{
  update_sideset_vector();
  sidesets.clear();
  sidesets.reserve(m_sidesets.size());

  for(SideSet* sideset : m_sidesets) {
    const Selector& elemSelector = get_cached_sideset_selector(sideset->get_part());
    stk::util::insert_keep_sorted_and_unique(SideSetSelector(*sideset->get_part(), sideset, &elemSelector), sidesets);
  }
}

void IncrementalSidesetUpdater::add_sidesets_from_part(const Part& part, std::vector<SideSet*>& sidesets)
{
  if(part.primary_entity_rank() == m_metaData.side_rank()) {
    const Part &parentPart = get_sideset_parent(part);
    if(m_bulkData.does_sideset_exist(parentPart)) {
      SideSet& sideset = m_bulkData.get_sideset(parentPart);
      stk::util::insert_keep_sorted_and_unique(&sideset, sidesets);
    }
  }
}

void IncrementalSidesetUpdater::fill_sidesets_from_part_ordinals(const OrdinalVector& parts, std::vector<SideSet*>& sidesets)
{
  update_sideset_vector();
  sidesets.clear();
  sidesets.reserve(m_sidesets.size());

  for(auto partOrdinal : parts) {
    Part& part = m_metaData.get_part(partOrdinal);
    add_sidesets_from_part(part, sidesets);
  }
}

void IncrementalSidesetUpdater::add_sidesets_from_parts(const PartVector& parts, std::vector<SideSet*>& sidesets)
{
  for(Part* part : parts) {
    add_sidesets_from_part(*part, sidesets);
  }
}

void IncrementalSidesetUpdater::add_sidesets_from_parts(const ConstPartVector& parts, std::vector<SideSet*>& sidesets)
{
  for(const Part* part : parts) {
    add_sidesets_from_part(*part, sidesets);
  }
}

void IncrementalSidesetUpdater::fill_sidesets_from_parts(const PartVector& parts, std::vector<SideSet*>& sidesets)
{
  update_sideset_vector();
  sidesets.clear();
  sidesets.reserve(m_sidesets.size());
  add_sidesets_from_parts(parts, sidesets);
}

namespace impl {
bool contains_elem_or_side_ranked_part(const MetaData& meta, const OrdinalVector& parts)
{
  const stk::mesh::EntityRank sideRank = meta.side_rank();

  for (unsigned ord : parts) {
    const stk::mesh::EntityRank partRank = meta.get_part(ord).primary_entity_rank();
    if ((partRank == stk::topology::ELEM_RANK) || (partRank == sideRank)) {
      return true;
    }
  }
  return false;
}
}

void IncrementalSidesetUpdater::entity_parts_removed(Entity entity, const OrdinalVector& parts)
{
  m_elemOrSideChangedRankedParts = m_elemOrSideChangedRankedParts || impl::contains_elem_or_side_ranked_part(m_metaData, parts);

  if(m_bulkData.entity_rank(entity) == m_metaData.side_rank()) {
    std::vector<SideSet*> sidesets;
    fill_sidesets_from_part_ordinals(parts, sidesets);
    m_helper.remove_side_entries_from_sidesets(sidesets, entity, nullptr);
  }

  if (m_bulkData.entity_rank(entity) == stk::topology::ELEM_RANK) {
    if (m_bulkData.bucket_ptr(entity) != nullptr && m_bulkData.num_sides(entity) > 0) {
      accumulate_element_block_part_changes(entity, parts);
    }
  }
}

void IncrementalSidesetUpdater::fill_sidesets_and_selectors_from_blocks_touching_surfaces(const OrdinalVector& parts,
                                                                                          SideSetSelectorVector& sidesetsAndSelectors)
{
  update_sideset_vector();
  sidesetsAndSelectors.clear();
  sidesetsAndSelectors.reserve(m_sidesets.size());
  constexpr bool checkSubsets = false;

  for(unsigned partOrdinal : parts) {
    Part& part = m_metaData.get_part(partOrdinal);

    if(part.primary_entity_rank() == stk::topology::ELEM_RANK) {
      const ConstPartVector& surfaces = get_cached_surfaces_touching_block(part);

      for(auto surface : surfaces) {
        add_sidesets_and_selectors_from_part(*surface, checkSubsets, sidesetsAndSelectors);
      }
    }
  }
}

void IncrementalSidesetUpdater::entity_parts_added(Entity entity, const OrdinalVector& parts)
{
  m_elemOrSideChangedRankedParts = m_elemOrSideChangedRankedParts || impl::contains_elem_or_side_ranked_part(m_metaData, parts);

  if(m_bulkData.entity_rank(entity) == m_metaData.side_rank()) {
    SideSetSelectorVector sidesets;
    fill_sidesets_and_selectors_from_part_ordinals(parts, sidesets);

    for (SideSetSelector& sideset : sidesets) {
      if(sideset.sideset()->get_accept_all_internal_non_coincident_entries()) {
        m_helper.add_side_entries_to_sideset(sideset, entity);
      } else {
        m_helper.add_coincident_side_entries_to_sideset(sideset, entity);
      }
    }
  }

  if(m_bulkData.entity_rank(entity) == stk::topology::ELEM_RANK) {
    if (m_bulkData.bucket_ptr(entity) != nullptr && m_bulkData.num_sides(entity) > 0) {
      SideSetSelectorVector sidesetsAndSelectors;
      fill_sidesets_and_selectors_from_blocks_touching_surfaces(parts, sidesetsAndSelectors);

      m_helper.add_sideset_entry_for_element_selected_by_sidesets(entity, sidesetsAndSelectors);
    }
  }
}

} } //namespaces
