
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

    const bool sortById = true;
    stk::mesh::get_entities(params.bulk_data(), type, selector, entities, sortById);
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

void throw_if_any_elem_block_has_invalid_topology(const stk::mesh::MetaData& meta,
                                                  const std::string& msgRegionName)
{
  const stk::mesh::PartVector& parts = meta.get_parts();
  for(const stk::mesh::Part* part : parts) {
    if (part->primary_entity_rank() == stk::topology::ELEM_RANK && stk::io::is_part_io_part(*part)) {
      if (part->topology() == stk::topology::INVALID_TOPOLOGY) {
        std::ostringstream msg;
        msg << " INTERNAL_ERROR when defining output for region '"<<msgRegionName
            <<"': Part "<<part->name()<<" returned INVALID from get_topology(). "
            <<"Please contact sierra-help@sandia.gov";
        throw std::runtime_error(msg.str());
      }
    }
  }
}

stk::mesh::FieldVector get_transient_fields(const stk::mesh::MetaData &meta)
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

stk::mesh::FieldVector get_transient_fields(const stk::mesh::MetaData &meta, const stk::mesh::EntityRank rank)
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

std::string construct_parallel_filename(const std::string &baseFilename, int numSubdomains, int subdomainIndex)
{
    int width = std::log10(static_cast<double>(numSubdomains))+1;
    std::ostringstream os;
    os << baseFilename << "." << numSubdomains << "." << std::setfill('0') << std::setw(width) << subdomainIndex;
    return os.str();
}

std::string construct_filename_for_serial_or_parallel(const std::string &baseFilename, int numSubdomains, int subdomainIndex)
{
    ThrowRequire(numSubdomains > 0);
    ThrowRequire(subdomainIndex >=0 && subdomainIndex<numSubdomains);
    if(numSubdomains == 1)
        return baseFilename;
    return stk::io::construct_parallel_filename(baseFilename, numSubdomains, subdomainIndex);
}

}}

