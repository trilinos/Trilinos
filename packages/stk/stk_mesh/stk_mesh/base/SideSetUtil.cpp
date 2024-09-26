#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SidesetUpdater.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include "stk_mesh/baseImpl/SideSetUtilImpl.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/diag/StringUtil.hpp"           // for Type, etc
#include "stk_util/util/string_case_compare.hpp"
#include <stk_util/environment/RuntimeWarning.hpp>
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include <algorithm>
#include <utility>

namespace stk {
namespace mesh {

bool is_elem_side_pair_in_sideset(const stk::mesh::SideSet& sset, stk::mesh::Entity elem, stk::mesh::ConnectivityOrdinal ordinal)
{
    bool isPresent = false;
    for(const stk::mesh::SideSetEntry& entry : sset)
    {   
        if(entry.element == elem && entry.side == ordinal)
        {
            isPresent = true;
            break;
        }
    }   

    return isPresent;
}

std::pair<bool,bool> old_is_positive_sideset_polarity(const stk::mesh::BulkData &bulk,
                                                  const stk::mesh::Part& sideSetPart,
                                                  stk::mesh::Entity face,
                                                  const stk::mesh::Part* activePart)
{
    std::pair<bool,bool> returnValue(false,false);

    stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
    STK_ThrowRequire(bulk.entity_rank(face) == sideRank);
    STK_ThrowAssert(bulk.bucket(face).member(sideSetPart));

    const stk::mesh::Part &parentPart = stk::mesh::get_sideset_parent(sideSetPart);

    const stk::mesh::SideSet* ssetPtr = nullptr;
    try {
        ssetPtr = &bulk.get_sideset(parentPart);
    }
    catch(std::exception&) {
        return returnValue;
    }

    const stk::mesh::SideSet& sset = *ssetPtr;

    const stk::mesh::Entity* sideElements = bulk.begin_elements(face);
    const unsigned numSideElements = bulk.num_elements(face);

    int numFound = 0;

    bool foundElemWithPermutationZero = false;

    stk::mesh::Entity foundElem;

    stk::mesh::Selector activeSelector = (activePart == nullptr) ? bulk.mesh_meta_data().universal_part() : *activePart;

    for(unsigned i=0; i<numSideElements; ++i)
    {
        stk::mesh::Entity elem = sideElements[i];

        if (!activeSelector(bulk.bucket(elem))) {
            continue;
        }

        const stk::mesh::Entity * elem_sides = bulk.begin(elem, sideRank);
        stk::mesh::ConnectivityOrdinal const * side_ordinal = bulk.begin_ordinals(elem, sideRank);
        stk::mesh::Permutation const * side_permutations = bulk.begin_permutations(elem, sideRank);
        const unsigned num_elem_sides = bulk.num_connectivity(elem, sideRank);

        for(unsigned k = 0; k < num_elem_sides; ++k)
        {
            if(elem_sides[k] == face)
            {
                if(is_elem_side_pair_in_sideset(sset, elem, side_ordinal[k])) {
                    numFound++;
                    foundElem = elem;

                    if (side_permutations[k] == 0) {
                        foundElemWithPermutationZero = true;
                    }
                }
            }
        }
    }

    returnValue.first = true;

    if(numFound == 1) {
        returnValue.second = foundElemWithPermutationZero;
    } else if(numFound == 0 && numSideElements == 1) {

        stk::mesh::Entity elem = sideElements[0];
        const stk::mesh::Entity * elem_sides = bulk.begin(elem, sideRank);
        stk::mesh::Permutation const * side_permutations = bulk.begin_permutations(elem, sideRank);
        const size_t num_elem_sides = bulk.num_connectivity(elem, sideRank);

        for(size_t k = 0; k < num_elem_sides; ++k)
        {
            if(elem_sides[k] == face)
            {
                if (side_permutations[k] == 0) {
                    foundElemWithPermutationZero = true;
                }
            }
        }

        returnValue.second = !foundElemWithPermutationZero;
    }

    return returnValue;
}

stk::mesh::EntityVector get_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Part& sidesetPart)
{
    stk::mesh::MetaData &meta = bulkData.mesh_meta_data();
    stk::mesh::EntityVector sides;
    stk::mesh::Selector sideSelector = sidesetPart & ( meta.locally_owned_part() | meta.globally_shared_part());
    const bool sortByGlobalId = false;
    stk::mesh::EntityRank sideRank = sidesetPart.primary_entity_rank() == stk::topology::INVALID_RANK ? meta.side_rank() : sidesetPart.primary_entity_rank();
    stk::mesh::get_entities(bulkData, sideRank, sideSelector, sides, sortByGlobalId);
    return sides;
}

void fill_sideset(const stk::mesh::Part& sidesetPart, stk::mesh::BulkData& bulkData, const stk::mesh::Selector& elementSelector)
{
    stk::mesh::SideSet *sideSet = nullptr;
    if(sidesetPart.subsets().empty()) {
        const stk::mesh::Part &parentPart = get_sideset_parent(sidesetPart);

        bool sidesetExists = bulkData.does_sideset_exist(parentPart);
        if(sidesetExists)
            sideSet = &bulkData.get_sideset(parentPart);
        else
            sideSet = &bulkData.create_sideset(parentPart);

        stk::mesh::EntityVector sides = get_sides(bulkData, parentPart);
        std::vector<stk::mesh::SideSetEntry> newSides;
        newSides.reserve(sides.size());

        for(stk::mesh::Entity side : sides)
        {
            unsigned numElements = bulkData.num_elements(side);
            const stk::mesh::Entity* elements = bulkData.begin_elements(side);
            const stk::mesh::ConnectivityOrdinal *ordinals = bulkData.begin_element_ordinals(side);

            for(unsigned i=0;i<numElements;++i) {
                const stk::mesh::Bucket& bucket = bulkData.bucket(elements[i]);
                bool isOwned = bucket.owned();
                bool isSelected = elementSelector(bucket);

                if(isOwned && isSelected) {
                    stk::mesh::ConnectivityOrdinal sideOrdOffset = bulkData.entity_rank(side) == stk::topology::FACE_RANK ? 0 : bucket.topology().num_faces();
                    newSides.emplace_back(elements[i], ordinals[i] + sideOrdOffset);
                }
            }
        }

        sideSet->add(newSides);
    }
}

bool is_face_represented_in_sideset(const stk::mesh::BulkData& bulk, const stk::mesh::Entity face, const stk::mesh::SideSet& sset)
{
    stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
    STK_ThrowRequire(bulk.entity_rank(face) == sideRank);

    std::vector<stk::mesh::Entity> side_elements;

    stk::mesh::ConnectedEntities sideNodes = bulk.get_connected_entities(face, stk::topology::NODE_RANK);
    stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEM_RANK, sideNodes.size(), sideNodes.data(), side_elements);

    bool found = false;

    for(stk::mesh::Entity elem : side_elements)
    {
        const stk::mesh::Entity * elem_sides = bulk.begin(elem, sideRank);
        stk::mesh::ConnectivityOrdinal const * side_ordinal = bulk.begin_ordinals(elem, sideRank);
        const size_t num_elem_sides = bulk.num_connectivity(elem, sideRank);

        for(size_t k = 0; k < num_elem_sides; ++k)
        {
            if(elem_sides[k] == face)
            {
                if(sset.contains(elem, side_ordinal[k])) {
                    found = true;
                    break;
                }
            }
        }

        if(found) {
            break;
        }
    }

    return found;
}

bool should_reconstruct_sideset(const stk::mesh::BulkData& bulkData, const stk::mesh::Part& surfacePart)
{
    bool reconstructSideset = false;

    const stk::mesh::Part &parentPart = get_sideset_parent(surfacePart);
    if (bulkData.does_sideset_exist(parentPart))
    {
        const stk::mesh::SideSet& ss = bulkData.get_sideset(parentPart);
        if (ss.is_from_input())
        {
            reconstructSideset = false;
        }

        for(const stk::mesh::SideSetEntry &entry : ss)
        {
            stk::mesh::Entity elem = entry.element;

            if(!bulkData.is_valid(elem) || !bulkData.bucket(elem).owned())
            {
                reconstructSideset = true;
                break;
            }
        }

        stk::mesh::EntityVector faces;
        stk::topology::rank_t sideRank = bulkData.mesh_meta_data().side_rank();
        stk::mesh::get_selected_entities(surfacePart, bulkData.buckets(sideRank), faces);

        for(stk::mesh::Entity face : faces)
        {
            bool modifiedState    = (bulkData.state(face) != stk::mesh::Unchanged);
            bool emptySideset     = (ss.size() == 0);
            bool faceNotInSideset = (emptySideset || !is_face_represented_in_sideset(bulkData, face, ss));

            if(modifiedState && faceNotInSideset)
            {
                reconstructSideset = true;
                break;
            }
        }
    }
    else
    {
        reconstructSideset = true;
    }

    return reconstructSideset;
}

void reconstruct_sideset(stk::mesh::BulkData& bulkData, const stk::mesh::Part& surfacePart)
{
    bulkData.clear_sideset(surfacePart);
    std::vector<const stk::mesh::Part *> touching_parts = bulkData.mesh_meta_data().get_blocks_touching_surface(&surfacePart);

    stk::mesh::Selector elementSelector = stk::mesh::selectUnion(touching_parts);
    fill_sideset(surfacePart, bulkData, elementSelector);
}

void create_bulkdata_sidesets(stk::mesh::BulkData& bulkData)
{
    if(bulkData.was_mesh_modified_since_sideset_creation())
    {
        std::vector<const stk::mesh::Part *> surfacesInMap = bulkData.mesh_meta_data().get_surfaces_in_surface_to_block_map();

        bool reconstructSideset = false;
        for(size_t i=0;i<surfacesInMap.size();++i)
        {
            const stk::mesh::Part* surfacePart = surfacesInMap[i];
            reconstructSideset = should_reconstruct_sideset(bulkData, *surfacePart) || reconstructSideset;
        }
        if (reconstructSideset)
        {
            for(size_t i=0;i<surfacesInMap.size();++i)
            {
                const stk::mesh::Part* surfacePart = surfacesInMap[i];
                reconstruct_sideset(bulkData, *surfacePart);
            }
        }
    }
}

std::vector<const stk::mesh::Part*> get_sideset_io_parts(const stk::mesh::BulkData& bulkData, stk::mesh::Entity face)
{
    STK_ThrowRequire(bulkData.entity_rank(face) == bulkData.mesh_meta_data().side_rank());
    const stk::mesh::PartVector& parts = bulkData.bucket(face).supersets();

    std::vector<const stk::mesh::Part *> sidesetParts;

    for(const stk::mesh::Part *part : parts)
    {
        if(stk::mesh::is_side_set(*part))
        {
            sidesetParts.push_back(part);
        }
    }

    return sidesetParts;
}

bool are_elements_in_different_element_blocks(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements)
{
    STK_ThrowRequireWithSierraHelpMsg(elements.size()>1);

    stk::mesh::Selector elementBlockSelector = *get_element_block_part(bulkData, elements[0]);

    bool allElementsSameBlock = true;
    for(unsigned i=1;i<elements.size();++i)
    {
        if(!elementBlockSelector(bulkData.bucket(elements[i])))
        {
            allElementsSameBlock = false;
            break;
        }
    }

    return !allElementsSameBlock;
}


bool is_side_supported(const stk::mesh::BulkData &bulk, stk::mesh::Entity side, const stk::mesh::impl::ParallelPartInfo& parallelPartInfo)
{
    const stk::mesh::EntityVector elements(bulk.begin_elements(side), bulk.end_elements(side));
    const stk::mesh::OrdinalVector ordinals(bulk.begin_ordinals(side, stk::topology::ELEM_RANK), bulk.end_ordinals(side, stk::topology::ELEM_RANK));

    bool supported = true;

    int indexFirstOwnedElement = -1;
    for(size_t i=0;i<elements.size();++i)
    {
        if(bulk.bucket(elements[i]).owned())
        {
            indexFirstOwnedElement = i;
            break;
        }
    }

    if(indexFirstOwnedElement >= 0)
    {
        if(elements.size()>1)
            supported = are_elements_in_different_element_blocks(bulk, elements);

        const stk::mesh::ElemElemGraph &graph = bulk.get_face_adjacent_element_graph();
        const stk::mesh::Part* elementBlockPart = get_element_block_part(bulk, elements[indexFirstOwnedElement]);
        stk::mesh::PartOrdinal partOrdinal = elementBlockPart->mesh_meta_data_ordinal();
        stk::mesh::impl::LocalId elementId = graph.get_local_element_id(elements[indexFirstOwnedElement]);
        stk::mesh::Ordinal ord = ordinals[indexFirstOwnedElement];

        stk::mesh::GraphEdgesForElement graphEdges = graph.get_edges_for_element(elementId);
        for(size_t i = 0; i < graphEdges.size(); ++i)
        {
            const stk::mesh::GraphEdge& graphEdge =  graph.get_graph().get_edge_for_element(elementId, i);
            stk::mesh::Ordinal sideOrdinal = graphEdge.side1();
            if(!graph.is_connected_elem_locally_owned(elements[indexFirstOwnedElement], i) && ord == sideOrdinal)
            {
                const auto iter = parallelPartInfo.find(graphEdge.elem2());
                STK_ThrowRequireWithSierraHelpMsg(iter!=parallelPartInfo.end());
                const std::vector<stk::mesh::PartOrdinal>& other_element_part_ordinals = iter->second.elementPartOrdinals;
                bool found = std::binary_search(other_element_part_ordinals.begin(), other_element_part_ordinals.end(), partOrdinal);
                if(found)
                {
                    supported = false;
                    break;
                }
            }
        }
    }

    return supported;
}

bool does_not_contain_internal_sideset(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &sides, const stk::mesh::impl::ParallelPartInfo &parallelPartInfo)
{
    bool isSidesetSupported = true;
    for(stk::mesh::Entity side : sides)
       if(!is_side_supported(bulk, side, parallelPartInfo))
           isSidesetSupported = false;

   return stk::is_true_on_all_procs(bulk.parallel(), isSidesetSupported);
}


const stk::mesh::Part& get_sideset_parent(const stk::mesh::Part& sidesetPart)
{
    for(stk::mesh::Part * part : sidesetPart.supersets()) {
        if((part->id() != stk::mesh::Part::INVALID_ID) &&
           (part->id() == sidesetPart.id()) &&
           (part->primary_entity_rank() == sidesetPart.primary_entity_rank())) {
            return *part;
        }
    }

    return sidesetPart;
}



std::pair<bool,bool> is_positive_sideset_polarity(const stk::mesh::BulkData &bulk,
                                                  const stk::mesh::Part& sideSetPart,
                                                  stk::mesh::Entity face,
                                                  const stk::mesh::Part* activePart,
                                                  const stk::mesh::SideSet* inputSidesetPtr)
{
  STK_ThrowRequireMsg(bulk.has_face_adjacent_element_graph(), "ERROR: Bulk data does not have face adjacency graph enabled");

  const Selector activeSelector = (activePart == nullptr) ? bulk.mesh_meta_data().universal_part() : *activePart;

  std::vector<std::shared_ptr<SidesetUpdater>> updaters = bulk.get_observer_type<SidesetUpdater>();
  STK_ThrowRequireMsg(updaters.size() > 0, "ERROR: Bulk data does not have a SidesetUpdater enabled");

  impl::PolarityDataCache polarityDataCache(bulk, activeSelector, updaters[0]->get_parallel_part_info());

  return impl::internal_is_positive_sideset_polarity(sideSetPart, face, polarityDataCache, inputSidesetPtr);
}

std::pair<bool,bool> is_positive_sideset_face_polarity(const stk::mesh::BulkData &bulk, stk::mesh::Entity face,
                                                       const stk::mesh::Part* activePart)
{
  STK_ThrowRequireMsg(bulk.has_face_adjacent_element_graph(), "ERROR: Bulk data does not have face adjacency graph enabled");

  const Selector activeSelector = (activePart == nullptr) ? bulk.mesh_meta_data().universal_part() : *activePart;

  std::vector<std::shared_ptr<SidesetUpdater>> updaters = bulk.get_observer_type<SidesetUpdater>();
  STK_ThrowRequireMsg(updaters.size() > 0, "ERROR: Bulk data does not have a SidesetUpdater enabled");

  impl::PolarityDataCache polarityDataCache(bulk, activeSelector, updaters[0]->get_parallel_part_info());

  return impl::internal_is_positive_sideset_face_polarity(face, polarityDataCache);
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

}
}
