
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
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/diag/StringUtil.hpp"           // for Type, etc
#include "stk_util/util/string_case_compare.hpp"
#include <stk_util/environment/RuntimeWarning.hpp>

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################




namespace stk {
namespace io {


stk::mesh::Selector internal_build_selector(const stk::mesh::Selector *subset_selector,
                                            const stk::mesh::Selector *output_selector,
                                            const stk::mesh::Selector *shared_selector,
                                            const stk::mesh::Part &part,
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
                    const stk::mesh::Part &part,
                    stk::mesh::EntityRank type,
                    stk::mesh::EntityVector &entities,
                    bool include_shared)
{
    stk::mesh::Selector selector =  internal_build_selector(params.get_subset_selector(),
                                                            params.get_output_selector(type),
                                                            params.get_shared_selector(),
                                                            part,
                                                            include_shared);

    get_selected_nodes(params, selector, entities);
    return entities.size();
}

size_t get_entities(stk::io::OutputParams &params,
                        const stk::mesh::Part &part,
                        stk::mesh::EntityRank type,
                        stk::mesh::EntityVector &entities,
                        bool include_shared)
{
    stk::mesh::Selector selector =  internal_build_selector(params.get_subset_selector(),
                                                            params.get_output_selector(type),
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
        if (block)
        {
            blocks.push_back(block);
        }
    }
}

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

            for(unsigned i=0;i<numElements;++i) {
                bool isOwned = bulkData.bucket(elements[i]).owned();
                bool isSelected = elementSelector(bulkData.bucket(elements[i]));
                if(isOwned && isSelected) {
                    (*sideSet).add(stk::mesh::SideSetEntry{elements[i], ordinals[i]});
                }
            }
        }

        if(sidesetExists)
            stk::util::sort_and_unique(*sideSet);
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
                if(is_elem_side_pair_in_sideset(sset, elem, side_ordinal[k])) {
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

    const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(surfacePart);
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

bool storage_type_is_general(const std::string &storage)
{
    bool value = false;

    if(stk::equal_case(storage,"scalar")          ||
       stk::equal_case(storage,"vector_2d")       ||
       stk::equal_case(storage,"vector_3d")       ||
       stk::equal_case(storage,"full_tensor_36")  ||
       stk::equal_case(storage,"full_tensor_32")  ||
       stk::equal_case(storage,"full_tensor_22")  ||
       stk::equal_case(storage,"full_tensor_12")  ||
       stk::equal_case(storage,"sym_tensor_33")   ||
       stk::equal_case(storage,"sym_tensor_31")   ||
       stk::equal_case(storage,"sym_tensor_21")   ||
       stk::equal_case(storage,"matrix_22")       ||
       stk::equal_case(storage,"matrix_33")) {
        value = true;
    }

    return value;
}

bool storage_type_is_real(const std::string &storage)
{
    bool value = false;

    if(storage_type_is_general(storage) ||
       stk::equal_case(storage,"real")) {
        value = true;
    } else {
        std::string str = storage.substr(0, 4);
        if(stk::equal_case(str,"real")) {
            value = true;
        }
    }

    return value;
}

bool storage_type_is_integer(const std::string &storage)
{
    bool value = false;

    if(storage_type_is_general(storage)   ||
       stk::equal_case(storage,"integer") ||
       stk::equal_case(storage,"integer64")) {
        value = true;
    } else {
        std::string str = storage.substr(0, 7);
        if(stk::equal_case(str,"integer")) {
            value = true;
        }
    }

    return value;
}

std::pair<size_t, stk::util::ParameterType::Type> parse_square_bracket_case(const std::string &storage,
                                                                            stk::util::ParameterType::Type scalar,
                                                                            stk::util::ParameterType::Type vector)
{
    std::pair<size_t, stk::util::ParameterType::Type> type = std::make_pair(0, stk::util::ParameterType::INVALID);

    std::string storageType = sierra::make_lower(storage);

    // Step 0:
    // See if the type contains '[' and ']'
    char const *typestr = storageType.c_str();
    char const *lbrace  = std::strchr(typestr, '[');
    char const *rbrace  = std::strrchr(typestr, ']');

    if (lbrace != nullptr && rbrace != nullptr) {
        // Step 1:
        // First, we split off the basename (REAL/INTEGER) from the component count
        // ([2])
        // and see if the basename is a valid variable type and the count is a
        // valid integer.

        char *base = std::strtok(const_cast<char *>(typestr), "[]");
        assert(base != nullptr);

        if(!strcmp(base, "real") || !strcmp(base, "integer")) {
            char *num_str = std::strtok(nullptr, "[]");
            assert(num_str != nullptr);

            std::istringstream os(num_str);
            int n = 0;
            os >> n;

            type = std::make_pair(size_t(n), vector);
        }
    }

    return type;
}

std::pair<size_t, stk::util::ParameterType::Type> get_parameter_type_from_storage(const std::string &storage,
                                                                                  stk::util::ParameterType::Type scalar,
                                                                                  stk::util::ParameterType::Type vector)
{
    std::pair<size_t, stk::util::ParameterType::Type> type = std::make_pair(0, stk::util::ParameterType::INVALID);

    if(stk::equal_case(storage,"scalar")) {
        type = std::make_pair(1, scalar);
    } else if(stk::equal_case(storage,"integer")) {
        type = std::make_pair(1, scalar);
    } else if(stk::equal_case(storage,"real")) {
        type = std::make_pair(1, scalar);
    }else if(stk::equal_case(storage,"vector_2d")) {
        type = std::make_pair(2, vector);
    } else if(stk::equal_case(storage,"vector_3d")) {
        type = std::make_pair(3, vector);
    } else if(stk::equal_case(storage,"full_tensor_36")) {
        type = std::make_pair(9, vector);
    } else if(stk::equal_case(storage,"full_tensor_32")) {
        type = std::make_pair(5, vector);
    } else if(stk::equal_case(storage,"full_tensor_22")) {
        type = std::make_pair(4, vector);
    } else if(stk::equal_case(storage,"full_tensor_12")) {
        type = std::make_pair(3, vector);
    } else if(stk::equal_case(storage,"sym_tensor_33")) {
        type = std::make_pair(6, vector);
    } else if(stk::equal_case(storage,"sym_tensor_31")) {
        type = std::make_pair(4, vector);
    } else if(stk::equal_case(storage,"sym_tensor_21")) {
        type = std::make_pair(3, vector);
    } else if(stk::equal_case(storage,"matrix_33")) {
        type = std::make_pair(9, vector);
    } else if(stk::equal_case(storage,"matrix_22")) {
        type = std::make_pair(4, vector);
    } else {
        type = parse_square_bracket_case(storage, scalar, vector);
    }

    return type;
}

std::pair<size_t, stk::util::ParameterType::Type> get_parameter_type_from_field_representation(const std::string &storage,
                                                                                               Ioss::Field::BasicType dataType,
                                                                                               int copies)
{
    std::pair<size_t, stk::util::ParameterType::Type> type = std::make_pair(0, stk::util::ParameterType::INVALID);

    stk::util::ParameterType::Type scalar = stk::util::ParameterType::INVALID;
    stk::util::ParameterType::Type vector = stk::util::ParameterType::INVALID;

    if((dataType == Ioss::Field::DOUBLE) && storage_type_is_real(storage)) {
        scalar = stk::util::ParameterType::DOUBLE;
        vector = stk::util::ParameterType::DOUBLEVECTOR;
        type = get_parameter_type_from_storage(storage, scalar, vector);
    } else if((dataType == Ioss::Field::INTEGER) && storage_type_is_integer(storage)) {
        scalar = stk::util::ParameterType::INTEGER;
        vector = stk::util::ParameterType::INTEGERVECTOR;
        type = get_parameter_type_from_storage(storage, scalar, vector);
    } else if((dataType == Ioss::Field::INT64) && storage_type_is_integer(storage)) {
        scalar = stk::util::ParameterType::INT64;
        vector = stk::util::ParameterType::INT64VECTOR;
        type = get_parameter_type_from_storage(storage, scalar, vector);
    }

    if(copies > 1) {
        type.second = vector;
    }

    type.first *= copies;

    return type;
}

std::pair<size_t, Ioss::Field::BasicType> get_io_parameter_size_and_type(const stk::util::ParameterType::Type type,
                                                                         const boost::any &value)
{
  try {
    switch(type)  {
    case stk::util::ParameterType::INTEGER: {
      return std::make_pair(1, Ioss::Field::INTEGER);
    }

    case stk::util::ParameterType::INT64: {
      return std::make_pair(1, Ioss::Field::INT64);
    }

    case stk::util::ParameterType::DOUBLE: {
      return std::make_pair(1, Ioss::Field::REAL);
    }

    case stk::util::ParameterType::DOUBLEVECTOR: {
      std::vector<double> vec = boost::any_cast<std::vector<double> >(value);
      return std::make_pair(vec.size(), Ioss::Field::REAL);
    }

    case stk::util::ParameterType::INTEGERVECTOR: {
      std::vector<int> vec = boost::any_cast<std::vector<int> >(value);
      return std::make_pair(vec.size(), Ioss::Field::INTEGER);
    }

    case stk::util::ParameterType::INT64VECTOR: {
      std::vector<int64_t> vec = boost::any_cast<std::vector<int64_t> >(value);
      return std::make_pair(vec.size(), Ioss::Field::INT64);
    }

    default: {
      return std::make_pair(0, Ioss::Field::INVALID);
    }
    }
  }
  catch(...) {
    std::cerr << "ERROR: Actual type of parameter does not match the declared type. Something went wrong."
              << " Maybe you need to call add_global_ref() instead of add_global() or vice-versa.";
    throw;
  }
}

std::pair<bool,bool> is_positive_sideset_polarity(const stk::mesh::BulkData &bulk,
                                                  const stk::mesh::Part& sideSetPart,
                                                  stk::mesh::Entity face)
{
    std::pair<bool,bool> returnValue(false,false);

    stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
    ThrowRequire(bulk.entity_rank(face) == sideRank);
    ThrowRequire(bulk.bucket(face).member(sideSetPart));

    const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(sideSetPart);

    if(!bulk.does_sideset_exist(parentPart))
    {
        stk::RuntimeWarning() << "Sideset parent part: " << parentPart.name() << " does not exist for surface part: " << sideSetPart.name();
        return returnValue;
    }

    const stk::mesh::SideSet& sset = bulk.get_sideset(parentPart);

    std::vector<stk::mesh::Entity> side_elements;
    std::vector<stk::mesh::Entity> side_nodes(bulk.begin_nodes(face), bulk.end_nodes(face));

    stk::mesh::get_entities_through_relations(bulk, side_nodes, stk::topology::ELEMENT_RANK, side_elements);

    int numFound = 0;

    bool foundElemWithPermutationZero = false;

    stk::mesh::Entity foundElem;

    for(stk::mesh::Entity elem : side_elements)
    {
        const stk::mesh::Entity * elem_sides = bulk.begin(elem, sideRank);
        stk::mesh::ConnectivityOrdinal const * side_ordinal = bulk.begin_ordinals(elem, sideRank);
        stk::mesh::Permutation const * side_permutations = bulk.begin_permutations(elem, sideRank);
        const size_t num_elem_sides = bulk.num_connectivity(elem, sideRank);

        for(size_t k = 0; k < num_elem_sides; ++k)
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
    } else if(numFound == 0 && side_elements.size() == 1) {

        stk::mesh::Entity elem = side_elements[0];
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

std::pair<bool,bool> is_positive_sideset_face_polarity(const stk::mesh::BulkData &bulk, stk::mesh::Entity face)
{
    std::pair<bool,bool> returnValue(false,false);
    std::vector<const stk::mesh::Part*> sidesetParts = stk::io::get_sideset_io_parts(bulk, face);

    if(sidesetParts.size() == 0) {
        return returnValue;
    }

    returnValue = stk::io::is_positive_sideset_polarity(bulk, *sidesetParts[0], face);

    for(unsigned i=1; i<sidesetParts.size(); ++i)
    {
        const stk::mesh::Part* sidesetPart = sidesetParts[i];
        std::pair<bool,bool> partPolarity = stk::io::is_positive_sideset_polarity(bulk, *sidesetPart, face);

        ThrowRequireMsg(partPolarity == returnValue,
                        "Polarity for face: " << bulk.identifier(face) << " on sideset: "
                                              << sidesetParts[0]->name()  << " does not match that on sideset: " << sidesetPart->name());
    }

    return returnValue;
}

void superset_mesh_parts(const stk::mesh::Part& part, stk::mesh::PartVector& supersetParts)
{
  bool is_io_part = stk::io::is_part_io_part(part);
  if(is_io_part) {
    for(auto& i : part.supersets()) {
      if(i == nullptr || stk::mesh::is_auto_declared_part(*i)) {
        continue;
      }

      supersetParts.push_back(i);

      stk::mesh::PartVector supersets2 = i->supersets();
      for(size_t j = 0; j < supersets2.size(); ++j) {
        if(supersets2[j] != nullptr && !stk::mesh::is_auto_declared_part(*supersets2[j])) {
          supersetParts.push_back(supersets2[j]);
        }
      }
    }
  }
}

stk::mesh::Selector construct_sideset_selector(stk::io::OutputParams &params)
{
    const mesh::BulkData &bulk_data = params.bulk_data();
    const stk::mesh::Selector *subset_selector = params.get_subset_selector();
    const stk::mesh::Selector *output_selector = params.get_output_selector(stk::topology::ELEM_RANK);

    stk::mesh::Selector selector = ( bulk_data.mesh_meta_data().locally_owned_part() | bulk_data.mesh_meta_data().globally_shared_part() );
    if(subset_selector)
        selector &= *subset_selector;
    if(output_selector)
        selector &= *output_selector;

    return selector;
}

}}

