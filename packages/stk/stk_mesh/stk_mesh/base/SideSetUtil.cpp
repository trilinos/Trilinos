#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SidesetUpdater.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp"
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
    ThrowRequire(bulk.entity_rank(face) == sideRank);
    ThrowAssert(bulk.bucket(face).member(sideSetPart));

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
                bool isOwned = bulkData.bucket(elements[i]).owned();
                bool isSelected = elementSelector(bulkData.bucket(elements[i]));
                if(isOwned && isSelected) {
                    newSides.emplace_back(elements[i], ordinals[i]);
                }
            }
        }

        sideSet->add(newSides);
    }
}

bool is_face_represented_in_sideset(const stk::mesh::BulkData& bulk, const stk::mesh::Entity face, const stk::mesh::SideSet& sset)
{
    stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
    ThrowRequire(bulk.entity_rank(face) == sideRank);

    std::vector<stk::mesh::Entity> side_elements;
    std::vector<stk::mesh::Entity> side_nodes(bulk.begin_nodes(face), bulk.end_nodes(face));

    stk::mesh::get_entities_through_relations(bulk, side_nodes, stk::topology::ELEMENT_RANK, side_elements);

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
    ThrowRequire(bulkData.entity_rank(face) == bulkData.mesh_meta_data().side_rank());
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

const stk::mesh::Part* get_element_block_selector_for_element(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element)
{
    const stk::mesh::Part* elementBlockPart = nullptr;;
    const stk::mesh::PartVector& parts = bulkData.bucket(element).supersets();
    unsigned block_counter = 0;
    for(const stk::mesh::Part *part : parts)
    {
        if(stk::mesh::is_element_block(*part))
        {
            elementBlockPart = part;
            block_counter++;
        }
    }

    bool elementIsAssociatedWithOnlyOneElementBlock = block_counter == 1;
    ThrowRequireWithSierraHelpMsg(elementIsAssociatedWithOnlyOneElementBlock);
    ThrowRequireWithSierraHelpMsg(elementBlockPart != nullptr);
    return elementBlockPart;
}


bool are_elements_in_different_element_blocks(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements)
{
    ThrowRequireWithSierraHelpMsg(elements.size()>1);

    stk::mesh::Selector elementBlockSelector = *get_element_block_selector_for_element(bulkData, elements[0]);

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
        const stk::mesh::Part* elementBlockPart = get_element_block_selector_for_element(bulk, elements[indexFirstOwnedElement]);
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
                ThrowRequireWithSierraHelpMsg(iter!=parallelPartInfo.end());
                const std::vector<stk::mesh::PartOrdinal>& other_element_part_ordinals = iter->second;
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
    const stk::mesh::MeshIndex& meshIndex = bulk.mesh_index(face);
    const stk::mesh::Bucket& bucket = *meshIndex.bucket;

    stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
    ThrowRequire(bucket.entity_rank() == sideRank);
    ThrowAssert(bucket.member(sideSetPart));

    const stk::mesh::SideSet* ssetPtr = inputSidesetPtr;
    if (ssetPtr == nullptr) {
      const stk::mesh::Part &parentPart = stk::mesh::get_sideset_parent(sideSetPart);
      if(bulk.does_sideset_exist(parentPart))
      {
          ssetPtr = &bulk.get_sideset(parentPart);
      }
      else
      {
          return std::make_pair(false,false);
      }
    }

    const stk::mesh::SideSet& sset = *ssetPtr;

    const unsigned bucketOrdinal = meshIndex.bucket_ordinal;
    const unsigned numSideElements = bucket.num_elements(bucketOrdinal);
    const stk::mesh::Entity* sideElements = bucket.begin_elements(bucketOrdinal);
    const stk::mesh::ConnectivityOrdinal * sideOrdinals = bucket.begin_element_ordinals(bucketOrdinal);
    const stk::mesh::Permutation * sidePermutations = bucket.begin_element_permutations(bucketOrdinal);

    if (numSideElements == 1) {
      //hopefully the most common case, fast/early return:
      const bool isActive = (activePart == nullptr) || bulk.bucket(sideElements[0]).member(*activePart);
      const bool inSideset = isActive && sset.contains(sideElements[0], sideOrdinals[0]);
      const bool permutationZero = sidePermutations[0]==0;
      const bool positivePolarity = (inSideset) ? permutationZero : !permutationZero;
      return std::make_pair(true, positivePolarity);
    }

    int numFound = 0;
    bool foundElemWithPermutationZero = false;

    for(unsigned i=0; i<numSideElements; ++i)
    {
        stk::mesh::Entity elem = sideElements[i];

        if ((activePart != nullptr) && !bulk.bucket(elem).member(*activePart)) {
            continue;
        }

        if(sset.contains(elem, sideOrdinals[i])) {
            numFound++;

            if (sidePermutations[i] == 0) {
                foundElemWithPermutationZero = true;
            }
        }
    }

    return std::make_pair(true, (numFound==1 ? foundElemWithPermutationZero : false));
}

std::pair<bool,bool> is_positive_sideset_face_polarity(const stk::mesh::BulkData &bulk, stk::mesh::Entity face,
                                                      const stk::mesh::Part* activePart)
{
    std::pair<bool,bool> returnValue(false,false);

    std::vector<const stk::mesh::Part*> sidesetParts = get_sideset_io_parts(bulk, face);

    if(sidesetParts.size() == 0) {
        return returnValue;
    }

    returnValue = is_positive_sideset_polarity(bulk, *sidesetParts[0], face, activePart);

    if (sidesetParts.size() > 1) {
        stk::mesh::Selector activeSelector = activePart == nullptr ? bulk.mesh_meta_data().universal_part() : *activePart;
        for(unsigned i=1; i<sidesetParts.size(); ++i)
        {
            const stk::mesh::Part* sidesetPart = sidesetParts[i];
            std::pair<bool,bool> partPolarity = is_positive_sideset_polarity(bulk, *sidesetPart, face, activePart);

            ThrowRequireMsg(partPolarity == returnValue,
                            "Polarity for face: " << bulk.identifier(face) << " on sideset: "
                                                  << sidesetParts[0]->name()  << " does not match that on sideset: " << sidesetPart->name());
        }
    }

    return returnValue;
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
