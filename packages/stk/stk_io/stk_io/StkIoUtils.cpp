
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "StkIoUtils.hpp"
#include <algorithm>
#include <utility>
#include "Ioss_Field.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "StkMeshIoBroker.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPDecl.hpp"                                   // for RCP
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/SideSetEntry.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/environment/ReportHandler.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include "stk_util/util/SortAndUnique.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################




namespace stk {
namespace io {


stk::mesh::Selector internal_build_selector(const stk::mesh::Selector *subset_selector,
                                            const stk::mesh::Selector *output_selector,
                                            const stk::mesh::Selector *shared_selector,
                                            stk::mesh::Part &part,
                                            bool include_shared)
{
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);

    stk::mesh::Selector own_share = meta.locally_owned_part();
    if(include_shared) {
        own_share |= (shared_selector == nullptr) ? meta.globally_shared_part() : *shared_selector;
    }

    stk::mesh::Selector selector = part & own_share;
    if(subset_selector)
        selector &= *subset_selector;
    if(output_selector)
        selector &= *output_selector;

    return selector;
}

size_t get_entities_for_nodeblock(stk::io::OutputParams &params,
                    stk::mesh::Part &part,
                    stk::mesh::EntityRank type,
                    stk::mesh::EntityVector &entities,
                    bool include_shared)
{
    stk::mesh::Selector selector =  internal_build_selector(params.get_subset_selector(),
                                                            params.get_output_selector(),
                                                            params.get_shared_selector(),
                                                            part,
                                                            include_shared);

    get_selected_nodes(params.bulk_data(), params.io_region(), selector, entities);
    return entities.size();
}

size_t get_entities(stk::io::OutputParams &params,
                        stk::mesh::Part &part,
                        stk::mesh::EntityRank type,
                        stk::mesh::EntityVector &entities,
                        bool include_shared)
{
    stk::mesh::Selector selector =  internal_build_selector(params.get_subset_selector(),
                                                            params.get_output_selector(),
                                                            nullptr,
                                                            part,
                                                            include_shared);

    get_selected_entities(selector, params.bulk_data().buckets(type), entities);
    return entities.size();
}

stk::mesh::EntityRank part_primary_entity_rank(const stk::mesh::Part &part)
{
  if (stk::mesh::MetaData::get(part).universal_part() == part) {
    return stk::topology::NODE_RANK;
  }
  else {
    return part.primary_entity_rank();
  }
}

IossBlockMembership get_block_memberships(stk::io::StkMeshIoBroker& stkIo)
{
    IossBlockMembership blockMemberships;
    const Ioss::Region & io_region = *stkIo.get_input_io_region();
    for(auto sideset : io_region.get_sidesets())
    {
        std::vector<std::string> sideBlockNames;
        sideset->block_membership(sideBlockNames);
        blockMemberships[sideset->name()] = sideBlockNames;
        for (auto&& side_subset : sideset->get_side_blocks())
        {
          sideBlockNames.clear();
          side_subset->block_membership(sideBlockNames);
          blockMemberships[side_subset->name()] = sideBlockNames;
        }
    }

    return blockMemberships;
}

void fill_block_parts_given_names(const std::vector<std::string>& side_block_names,
                                              stk::mesh::MetaData& meta,
                                              std::vector<const stk::mesh::Part*>& blocks)
{
    for(const std::string& partName : side_block_names)
    {
        stk::mesh::Part* block = meta.get_part(partName);
        ThrowRequireWithSierraHelpMsg(block != nullptr);
        blocks.push_back(block);
    }
}

stk::mesh::EntityVector get_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Part& sidesetPart)
{
    stk::mesh::MetaData &meta = bulkData.mesh_meta_data();
    stk::mesh::EntityVector sides;
    stk::mesh::Selector sideSelector = sidesetPart & ( meta.locally_owned_part() | meta.globally_shared_part());
    stk::mesh::get_selected_entities(sideSelector, bulkData.buckets(meta.side_rank()), sides);
    return sides;
}


void fill_sideset(const stk::mesh::Part& sidesetPart, stk::mesh::BulkData& bulkData, stk::mesh::Selector elementSelector)
{
    stk::mesh::SideSet *sideSet = nullptr;
    if(sidesetPart.subsets().empty()) {
        const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(sidesetPart);

        bool sidesetExists = bulkData.does_sideset_exist(parentPart);
        if(sidesetExists)
            sideSet = &bulkData.get_sideset(parentPart);
        else
            sideSet = &bulkData.create_sideset(parentPart);

        stk::mesh::EntityVector sides = get_sides(bulkData, parentPart);
        for(stk::mesh::Entity side : sides)
        {
            unsigned numElements = bulkData.num_elements(side);
            const stk::mesh::Entity* elements = bulkData.begin_elements(side);
            const stk::mesh::ConnectivityOrdinal *ordinals = bulkData.begin_element_ordinals(side);

            for(unsigned i=0;i<numElements;++i)
                if(bulkData.bucket(elements[i]).owned() && elementSelector(bulkData.bucket(elements[i])))
                    (*sideSet).push_back(stk::mesh::SideSetEntry{elements[i], ordinals[i]});
        }

        if(sidesetExists)
            stk::util::sort_and_unique(*sideSet);
    }
}

void create_bulkdata_sidesets(stk::mesh::BulkData& bulkData)
{
    if(bulkData.was_mesh_modified_since_sideset_creation())
    {
        bulkData.clear_sidesets();

        std::vector<const stk::mesh::Part *> surfacesInMap = bulkData.mesh_meta_data().get_surfaces_in_surface_to_block_map();
        for(size_t i=0;i<surfacesInMap.size();++i)
        {
            std::vector<const stk::mesh::Part *> touching_parts = bulkData.mesh_meta_data().get_blocks_touching_surface(surfacesInMap[i]);

            stk::mesh::Selector elementSelector = stk::mesh::selectUnion(touching_parts);
            fill_sideset(*surfacesInMap[i], bulkData, elementSelector);
        }
    }
}

const stk::mesh::Part* getElementBlockSelectorForElement(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element)
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


bool areElementsInDifferentElementBlocks(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements)
{
    ThrowRequireWithSierraHelpMsg(elements.size()>1);

    stk::mesh::Selector elementBlockSelector = *getElementBlockSelectorForElement(bulkData, elements[0]);

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


bool isSideSupported(const stk::mesh::BulkData &bulk, stk::mesh::Entity side, const stk::mesh::impl::ParallelPartInfo& parallelPartInfo)
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
            supported = areElementsInDifferentElementBlocks(bulk, elements);

        const stk::mesh::ElemElemGraph &graph = bulk.get_face_adjacent_element_graph();
        const stk::mesh::Part* elementBlockPart = getElementBlockSelectorForElement(bulk, elements[indexFirstOwnedElement]);
        stk::mesh::PartOrdinal partOrdinal = elementBlockPart->mesh_meta_data_ordinal();
        stk::mesh::impl::LocalId elementId = graph.get_local_element_id(elements[indexFirstOwnedElement]);
        stk::mesh::Ordinal ord = ordinals[indexFirstOwnedElement];

        const stk::mesh::GraphEdgesForElement& graphEdges = graph.get_edges_for_element(elementId);
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

bool isSidesetSupported(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &sides, const stk::mesh::impl::ParallelPartInfo &parallelPartInfo)
{
    bool isSidesetSupported = true;
    for(stk::mesh::Entity side : sides)
       if(!isSideSupported(bulk, side, parallelPartInfo))
           isSidesetSupported = false;

   return stk::is_true_on_all_procs(bulk.parallel(), isSidesetSupported);
}


stk::mesh::FieldVector get_transient_fields(stk::mesh::MetaData &meta)
{
    stk::mesh::FieldVector fields;

    for(stk::mesh::FieldBase *field : meta.get_fields())
    {
        const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
        if(fieldRole != nullptr && *fieldRole == Ioss::Field::TRANSIENT)
            fields.push_back(field);
    }

    return fields;
}

stk::mesh::FieldVector get_transient_fields(stk::mesh::MetaData &meta, const stk::mesh::EntityRank rank)
{
    stk::mesh::FieldVector fields;

    for(stk::mesh::FieldBase *field : meta.get_fields())
    {
        if(field->entity_rank() == rank)
        {
            const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
            if(fieldRole != nullptr && *fieldRole == Ioss::Field::TRANSIENT)
                fields.push_back(field);
        }
    }

    return fields;
}


const stk::mesh::Part& get_sideset_parent(const stk::mesh::Part& sidesetPart)
{
    for(stk::mesh::Part * part : sidesetPart.supersets()) {
        bool hasSameId   = (part->id() == sidesetPart.id());
        bool hasSameRank = (part->primary_entity_rank() == sidesetPart.primary_entity_rank());
        bool isIoPart    = stk::io::is_part_io_part(*part);

        if(hasSameId && hasSameRank && isIoPart) {
            return *part;
        }
    }

    return sidesetPart;
}

template<typename DATA_TYPE>
void write_global_to_stk_io(stk::io::StkMeshIoBroker& stkIo, size_t dbIndex,
                            const std::string& externalName,
                            size_t component_count, const void* ptr)
{
    const DATA_TYPE* data = reinterpret_cast<const DATA_TYPE*>(ptr);
    std::vector<DATA_TYPE> vecData(data, data+component_count);
    stkIo.write_global(dbIndex, externalName, vecData);
}

template
void write_global_to_stk_io<int>(stk::io::StkMeshIoBroker& stkIo, size_t dbIndex,
                                 const std::string& externalName,
                                 size_t component_count, const void* ptr);

template
void write_global_to_stk_io<double>(stk::io::StkMeshIoBroker& stkIo, size_t dbIndex,
                                    const std::string& externalName,
                                    size_t component_count, const void* ptr);
}}

