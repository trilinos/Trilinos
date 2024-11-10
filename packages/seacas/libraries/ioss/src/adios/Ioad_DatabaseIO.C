// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <tokenize.h>

#include "Ioss_CodeTypes.h"       // for HAVE_MPI
#include "Ioss_CommSet.h"         // for CommSet
#include "Ioss_DBUsage.h"         // for DatabaseUsage, etc
#include "Ioss_DatabaseIO.h"      // for DatabaseIO
#include "Ioss_EdgeBlock.h"       // for EdgeBlock
#include "Ioss_EdgeSet.h"         // for EdgeSet
#include "Ioss_ElementBlock.h"    // for ElementBlock
#include "Ioss_ElementSet.h"      // for ElementSet
#include "Ioss_ElementTopology.h" // for NameList
#include "Ioss_EntityType.h"      // for EntityType::ELEMENTBLOCK
#include "Ioss_FaceBlock.h"       // for FaceBlock
#include "Ioss_FaceSet.h"         // for FaceSet
#include "Ioss_FileInfo.h"        // for FileInfo
#include "Ioss_Map.h"             // for Map, MapContainer
#include "Ioss_NodeBlock.h"       // for NodeBlock
#include "Ioss_NodeSet.h"         // for NodeSet
#include "Ioss_ParallelUtils.h"   // for ParallelUtils, etc
#include "Ioss_Property.h"        // for Property
#include "Ioss_SerializeIO.h"     // for SerializeIO
#include "Ioss_SideBlock.h"       // for SideBlock
#include "Ioss_Utils.h"           // for Utils, IOSS_ERROR, etc
#include <exodusII.h>

#include "adios/Ioad_Constants.h"
#include "adios/Ioad_Helper.h"
#include "adios/Ioad_TemplateToValue.h"

#include <climits>

#include "adios/Ioad_DatabaseIO.h"

namespace Ioss {
  class PropertyManager;
}

namespace Ioad {

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, Ioss_MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props), rank(RankInit()),
        adios_wrapper(communicator, filename, is_input(), rank, props)
  {
    dbState = Ioss::STATE_UNKNOWN;
    // Always 64 bits
    dbIntSizeAPI = Ioss::USE_INT64_API;
    // Streaming initialization
    is_streaming            = adios_wrapper.IsStreaming();
    previous_time_streaming = -1.0;
  }

  // Used to force `rank` initialization before creating `adios_wrapper`.
  int DatabaseIO::RankInit()
  {
    Ioss::SerializeIO serializeIO_(this);
    number_proc = Ioss::SerializeIO::getSize();
    return Ioss::SerializeIO::getRank();
  }

  bool DatabaseIO::begin_nl(Ioss::State state)
  {
    // initialization
    Ioss::Region *this_region = get_region();

    // this_region->property_update("streaming_status", -1);

    dbState = state;
    if (state == Ioss::STATE_DEFINE_MODEL) {
      // Should `BeginStep()` only be if (!is_input() || is_streaming) ???

      adios2::StepStatus status = adios_wrapper.BeginStep();
      this_region->property_update("streaming_status", static_cast<int>(status));
    }
    return true;
  }

  void DatabaseIO::write_properties(const Ioss::GroupingEntity *const entity,
                                    const std::string                &encoded_name)
  {
    // Write field properties
    std::vector<std::string> property_list = properties_to_save(entity);
    for (auto &property_name : property_list) {
      Ioss::Property property = entity->get_property(property_name);

      std::string variable_name = get_property_variable_name(property_name);
      switch (property.get_type()) {
      case Ioss::Property::BasicType::REAL:
        adios_wrapper.PutMetaVariable<double>(variable_name, property.get_real(), encoded_name);
        break;
      case Ioss::Property::BasicType::INTEGER:
        adios_wrapper.PutMetaVariable<int>(variable_name, property.get_int(), encoded_name);
        break;
      case Ioss::Property::BasicType::STRING:
        adios_wrapper.PutMetaVariable<std::string>(variable_name, property.get_string(),
                                                   encoded_name);
        break;
      case Ioss::Property::BasicType::POINTER:
      case Ioss::Property::BasicType::INVALID:
      default:
        // Do not save properties with an invalid type.
        break;
      }
    }
  }

  // TODO: refactor and consolidate `write_meta_data_container` functions to avoid code duplication.
  template <typename T> int64_t DatabaseIO::write_meta_data_container(const T &entity_blocks)
  {
    int64_t count = 0;
    for (auto &entity_block : entity_blocks) {
      count += entity_block->entity_count();
      std::string entity_type  = entity_block->type_string();
      std::string entity_name  = entity_block->name();
      std::string encoded_name = encode_field_name({entity_type, entity_name});
      write_properties(entity_block, encoded_name);
      adios_wrapper.PutMetaVariable<std::string>(Topology_meta, entity_block->topology()->name(),
                                                 encoded_name);
    }
    return count;
  }

  template <>
  int64_t DatabaseIO::write_meta_data_container<Ioss::CommSetContainer>(
      const Ioss::CommSetContainer &entity_blocks)
  {
    for (auto &entity_block : entity_blocks) {
      std::string entity_type  = entity_block->type_string();
      std::string entity_name  = entity_block->name();
      std::string encoded_name = encode_field_name({entity_type, entity_name});
      write_properties(entity_block, encoded_name);
    }
    return 0;
  }

  std::pair<int64_t, int64_t>
  DatabaseIO::write_meta_data_sideblockcontainer(const Ioss::SideBlockContainer &entity_blocks)
  {
    int64_t count    = 0;
    int64_t df_count = 0;

    for (auto &entity_block : entity_blocks) {
      count += entity_block->entity_count();
      df_count += entity_block->get_property("distribution_factor_count").get_int();
      std::string entity_type  = entity_block->type_string();
      std::string entity_name  = entity_block->name();
      std::string encoded_name = encode_field_name({entity_type, entity_name});
      adios_wrapper.PutMetaVariable<std::string>(Topology_meta, entity_block->topology()->name(),
                                                 encoded_name);
      adios_wrapper.PutMetaVariable<std::string>(
          Parent_topology_meta, entity_block->parent_element_topology()->name(), encoded_name);
      write_properties(entity_block, encoded_name);
    }
    return {count, df_count};
  }

  template <>
  int64_t
  DatabaseIO::write_meta_data_container<Ioss::SideSetContainer>(const Ioss::SideSetContainer &ssets)
  {
    int64_t count = 0;

    for (auto &sset : ssets) {
      count += sset->entity_count();
      const std::string entity_type = sset->type_string();
      const std::string entity_name = sset->name();
      write_properties(sset, encode_field_name({entity_type, entity_name}));

      // side blocks
      const Ioss::SideBlockContainer &sblocks = sset->get_side_blocks();
      std::string encoded_name                = encode_sideblock_name(entity_type, entity_name);

      std::string stringified_sblock_names = stringify_side_block_names(sblocks);
      adios_wrapper.InquireAndPut<std::string>(encoded_name, &stringified_sblock_names);

      int sblock_count                 = 0;
      int df_count                     = 0;
      std::tie(sblock_count, df_count) = write_meta_data_sideblockcontainer(sblocks);
      Ioss::SideSet *new_entity        = const_cast<Ioss::SideSet *>(sset);
      new_entity->property_add(Ioss::Property("entity_count", sblock_count));
      new_entity->property_add(Ioss::Property("distribution_factor_count", df_count));
    }
    return count;
  }

  std::string DatabaseIO::encoded_coordinate_frame_name(Ioss::CoordinateFrame coordinate_frame)
  {
    int64_t id = coordinate_frame.id();

    char        tag = coordinate_frame.tag();
    std::string tag_str(1, tag);
    return encode_field_name({coordinate_frame_name, std::to_string(id), tag_str});
  }

  template <>
  int64_t DatabaseIO::write_meta_data_container<Ioss::CoordinateFrameContainer>(
      const Ioss::CoordinateFrameContainer &coordinate_frames)
  {
    int64_t count = 0;
    for (auto &coordinate_frame : coordinate_frames) {
      std::string encoded_name = encoded_coordinate_frame_name(coordinate_frame);
      adios_wrapper.InquireAndPut<double>(encoded_name, coordinate_frame.coordinates());
    }
    return count;
  }

  // Write data defined per block, not per field.
  void DatabaseIO::write_meta_data()
  {
    Ioss::Region                   *region      = get_region();
    const Ioss::NodeBlockContainer &node_blocks = region->get_node_blocks();

    assert(node_blocks.size() == 1);
    spatialDimension = node_blocks[0]->get_property("component_degree").get_int();

    // Region
    write_properties(region, encode_field_name({region->type_string(), region_name}));
    // Node blocks --
    nodeCount = write_meta_data_container<Ioss::NodeBlockContainer>(node_blocks);
    // Edge Blocks --
    const Ioss::EdgeBlockContainer &edge_blocks = region->get_edge_blocks();
    edgeCount = write_meta_data_container<Ioss::EdgeBlockContainer>(edge_blocks);
    // Face Blocks --
    const Ioss::FaceBlockContainer &face_blocks = region->get_face_blocks();
    faceCount = write_meta_data_container<Ioss::FaceBlockContainer>(face_blocks);
    // Element Blocks --
    const Ioss::ElementBlockContainer &element_blocks = region->get_element_blocks();
    elementCount = write_meta_data_container<Ioss::ElementBlockContainer>(element_blocks);
    // Side Blocks ...
    const Ioss::SideSetContainer &ssets = region->get_sidesets();
    write_meta_data_container<Ioss::SideSetContainer>(ssets);
    // Comm Sets --
    const Ioss::CommSetContainer &comm_sets = region->get_commsets();
    write_meta_data_container<Ioss::CommSetContainer>(comm_sets);
    // Coordinate frames --
    const Ioss::CoordinateFrameContainer &coordinate_frames = region->get_coordinate_frames();
    write_meta_data_container<Ioss::CoordinateFrameContainer>(coordinate_frames);
    // Global meta data
    put_data<unsigned long>(static_cast<void *>(&rank), Processor_id_meta);
  }

  void DatabaseIO::check_model()
  {
    Ioss::Region *region = get_region();

    // Add region missing properties.
    if (!region->property_exists("title")) {
      Ioss::Property title_property("title", "IOSS Default Output Title");
      region->property_add(title_property);
    }
  }

  bool DatabaseIO::end_nl(Ioss::State state)
  {
    // Transitioning out of state 'state'
    assert(state == dbState);
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL:
      if (!is_input()) {
        check_model();
        define_model();
        define_global_variables();
        write_meta_data();
      }
      break;
    case Ioss::STATE_DEFINE_TRANSIENT:
      if (!is_input()) {
        if (is_streaming) {
          // If streaming, we need to write the meta data available at
          // every step.
          write_meta_data();
        }
        Ioss::Field::RoleType role = Ioss::Field::RoleType::TRANSIENT;
        define_model(&role);
      }
      break;
    default: // ignore everything else...
      break;
    }

    {
      dbState = Ioss::STATE_UNKNOWN;
    }

    return true;
  }

  bool DatabaseIO::begin_state_nl(int state, double time)
  {
    Ioss::Region *this_region = get_region();
    if (!is_input()) {
      // `BeginStep()` should not be used at the same time as random access. Since at read time,
      // we currently read variables with random access, `BeginStep()` should not be used
      // at read time.
      // Begin  step for transient data
      adios2::StepStatus status = adios_wrapper.BeginStep();
      this_region->property_update("streaming_status", static_cast<int>(status));

      // Add time to adios
      adios_wrapper.PutMetaVariable<double>(Time_meta, time / timeScaleFactor);
    }
    else {
      if (is_streaming) {
        // Begin step for transient data if streaming. Otherwise, data will be accessed with
        // `SetStepSelection()`.
        adios2::StepStatus status = adios_wrapper.BeginStep();
        this_region->property_update("streaming_status", static_cast<int>(status));
      }
    }

    return true;
  }

  // common
  bool DatabaseIO::end_state_nl(int state, double time)
  {
    Ioss::Region *this_region = get_region();

    // if (!is_input()) {
    // End step for transient data
    // End step only for writing. Streaming is closed in `BeginStep()` call if time is different
    // from previous time.
    // if (!is_streaming || !is_input()) {
    adios_wrapper.EndStep();
    //}
    this_region->property_update("streaming_status", -1);

    //}
    return true;
  }

  void DatabaseIO::define_field_meta_variables(const std::string &encoded_name)
  {
    adios_wrapper.DefineMetaVariable<int>(Role_meta, encoded_name);
    adios_wrapper.DefineMetaVariable<std::string>(Var_type_meta, encoded_name);
  }

  template <typename T, typename>
  void DatabaseIO::define_entity_meta_variables(const std::string &encoded_name)
  {
    adios_wrapper.DefineMetaVariable<std::string>(Topology_meta, encoded_name);
  }

  template <typename T, typename, typename>
  void DatabaseIO::define_entity_meta_variables(const std::string &encoded_name)
  {
    // no op
  }

  template <>
  void DatabaseIO::define_entity_meta_variables<Ioss::SideBlock>(const std::string &encoded_name)
  {
    adios_wrapper.DefineMetaVariable<std::string>(Topology_meta, encoded_name);
    adios_wrapper.DefineMetaVariable<std::string>(Parent_topology_meta, encoded_name);
  }

  template <typename T>
  void DatabaseIO::define_model_internal(const Ioss::Field &field, const std::string &encoded_name,
                                         const std::string &entity_type,
                                         const std::string &field_name)
  {
    size_t component_count;
    size_t local_size;

    if (use_transformed_storage(field, entity_type, field_name)) {
      component_count = field.transformed_storage()->component_count();
      local_size      = field.transformed_count();
    }
    else {
      component_count = field.raw_storage()->component_count();
      local_size      = field.raw_count();
    }
    adios_wrapper.DefineVariable<T>(encoded_name, {number_proc, INT_MAX, component_count},
                                    {rank, 0, 0}, {1, local_size, component_count});
  }

  void DatabaseIO::define_properties(const Ioss::GroupingEntity *const entity_block,
                                     const std::string                &encoded_entity_name)
  {
    std::vector<std::string> property_list = properties_to_save(entity_block);
    for (auto &property_name : property_list) {
      Ioss::Property property      = entity_block->get_property(property_name);
      std::string    variable_name = get_property_variable_name(property_name);
      switch (property.get_type()) {
      case Ioss::Property::BasicType::REAL:
        adios_wrapper.DefineMetaVariable<double>(variable_name, encoded_entity_name);
        break;
      case Ioss::Property::BasicType::INTEGER:
        adios_wrapper.DefineMetaVariable<int>(variable_name, encoded_entity_name);
        break;
      case Ioss::Property::BasicType::STRING:
        adios_wrapper.DefineMetaVariable<std::string>(variable_name, encoded_entity_name);
        break;
      case Ioss::Property::BasicType::POINTER:
      case Ioss::Property::BasicType::INVALID:
      default:
        // Do not save properties with an invalid type.
        break;
      }
    }
  }

  template <typename T>
  void DatabaseIO::define_entity_internal(const T &entity_blocks, Ioss::Field::RoleType *role)
  {
    using cv_removed_value_type = typename std::remove_pointer_t<typename T::value_type>;
    for (auto &entity_block : entity_blocks) {
      std::string entity_type = entity_block->type_string();
      std::string entity_name = entity_block->name();
      if (!role) {
        // Meta variables and properties are only defined if state is STATE_DEFINE_MODEL:
        std::string encoded_entity_name = encode_field_name({entity_type, entity_name});
        define_entity_meta_variables<cv_removed_value_type>(encoded_entity_name);
        define_properties(entity_block, encoded_entity_name);
      }
      Ioss::NameList field_names = entity_block->field_describe();
      for (auto &field_name : field_names) {
        // Skip ignored fields
        if (find_field_in_mapset(entity_type, field_name, Ignore_fields)) {
          continue;
        }
        // Define entity block variables
        auto field = entity_block->get_fieldref(field_name);
        if (role && field.get_role() != *role) {
          continue;
        }
        std::string encoded_name = encode_field_name({entity_type, entity_name, field_name});
        switch (field.get_type()) {
        // REAL and DOUBLE have same value. Code is left here for documentation purposes
        // but commented out to compile.
        // case Ioss::Field::BasicType::REAL:
        case Ioss::Field::BasicType::DOUBLE:
          define_model_internal<double>(field, encoded_name, entity_type, field_name);
          break;
        case Ioss::Field::BasicType::INT32:
          define_model_internal<int32_t>(field, encoded_name, entity_type, field_name);
          break;
        case Ioss::Field::BasicType::INT64:
          define_model_internal<int64_t>(field, encoded_name, entity_type, field_name);
          break;
        case Ioss::Field::BasicType::COMPLEX:
          define_model_internal<Complex>(field, encoded_name, entity_type, field_name);
          break;
        case Ioss::Field::BasicType::CHARACTER:
          define_model_internal<char>(field, encoded_name, entity_type, field_name);
          break;
        case Ioss::Field::BasicType::STRING:
          adios_wrapper.DefineVariable<std::string>(
              encoded_name,
              // Global dimensions
              {field.raw_count()},
              // starting offset of the local array in the global space
              {0},
              // local size, could be defined later using SetSelection()
              {field.raw_count()});
          break;
        default:
          std::ostringstream errmsg;
          errmsg
              << "INTERNAL ERROR: Invalid field type. "
              << "Something is wrong in the Ioad::DatabaseIO::define_entity_internal() function. "
              << "Please report.\n";
          IOSS_ERROR(errmsg);
        }
        define_field_meta_variables(encoded_name);
      }
    }
  }

  void DatabaseIO::define_coordinate_frames_internal(
      const Ioss::CoordinateFrameContainer &coordinate_frames)
  {
    for (auto &coordinate_frame : coordinate_frames) {
      std::string encoded_name = encoded_coordinate_frame_name(coordinate_frame);
      adios_wrapper.DefineVariable<double>(encoded_name, {number_proc, 9, 1}, {rank, 0, 0},
                                           {1, 9, 1});
    }
  }

  // Similar to `write_meta_data()` function in other DatabaseIO. This function has been renamed in
  // this database to reflect more precisely what it accomplishes.
  void DatabaseIO::define_model(Ioss::Field::RoleType *role)
  {

    Ioss::Region                   *region      = get_region();
    const Ioss::NodeBlockContainer &node_blocks = region->get_node_blocks();

    // A single nodeblock named "nodeblock_1" will be created for the mesh. It contains information
    // for every node that exists in the model (Ioss-exodus-mapping.pdf).
    assert(node_blocks.size() == 1);

    if (!role) {
      // Schema
      adios_wrapper.DefineAttribute<unsigned int>(Schema_version_string, 1);
      // Region
      define_properties(region, encode_field_name({region->type_string(), region_name}));
    }
    // Node blocks --
    define_entity_internal(node_blocks, role);
    // Edge Blocks --
    const Ioss::EdgeBlockContainer &edge_blocks = region->get_edge_blocks();
    define_entity_internal(edge_blocks, role);
    // Face Blocks --
    const Ioss::FaceBlockContainer &face_blocks = region->get_face_blocks();
    define_entity_internal(face_blocks, role);
    // Element Blocks --
    const Ioss::ElementBlockContainer &element_blocks = region->get_element_blocks();
    define_entity_internal(element_blocks, role);
    // Nodesets ...
    const Ioss::NodeSetContainer &nodesets = region->get_nodesets();
    define_entity_internal(nodesets, role);
    // Edgesets ...
    const Ioss::EdgeSetContainer &edgesets = region->get_edgesets();
    define_entity_internal(edgesets, role);
    // Facesets ...
    const Ioss::FaceSetContainer &facesets = region->get_facesets();
    define_entity_internal(facesets, role);
    // Elementsets ...
    const Ioss::ElementSetContainer &elementsets = region->get_elementsets();
    define_entity_internal(elementsets, role);
    // CommSets ...
    const Ioss::CommSetContainer &csets = region->get_commsets();
    define_entity_internal(csets, role);
    // SideSets ...
    const Ioss::SideSetContainer &ssets = region->get_sidesets();
    define_entity_internal(ssets, role);
    // SideBlocks ...
    for (auto &sset : ssets) {
      std::string encoded_name = encode_sideblock_name(sset->type_string(), sset->name());
      const Ioss::SideBlockContainer &sblocks = sset->get_side_blocks();
      adios_wrapper.DefineVariable<std::string>(encoded_name);
      define_entity_internal(sblocks, role);
    }
    // Coordinate frames ...
    const Ioss::CoordinateFrameContainer &coordinate_frames = region->get_coordinate_frames();
    define_coordinate_frames_internal(coordinate_frames);
  }

  void DatabaseIO::define_global_variables()
  {
    adios_wrapper.DefineAttribute<double>(Time_scale_factor, timeScaleFactor);
    adios_wrapper.DefineMetaVariable<double>(Time_meta);
    adios_wrapper.DefineVariable<unsigned long>(Processor_id_meta, {number_proc}, {rank}, {1});
    adios_wrapper.DefineAttribute<unsigned long>(Processor_number_meta, number_proc);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(reg, field, data, data_size);
  }

  // Returns byte size of integers stored on the database...
  int DatabaseIO::int_byte_size_db() const { return 8; }

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::EDGEBLOCK | Ioss::FACEBLOCK | Ioss::ELEMENTBLOCK |
           Ioss::NODESET | Ioss::EDGESET | Ioss::FACESET | Ioss::ELEMENTSET | Ioss::SIDESET |
           Ioss::SIDEBLOCK | Ioss::REGION;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(nb, field, data, data_size);
  }

  template <typename T> void DatabaseIO::put_data(void *data, const std::string &encoded_name) const
  {
    T *rdata = static_cast<T *>(data);
    adios_wrapper.InquireAndPut<T>(encoded_name, rdata);
  }

  void DatabaseIO::put_meta_variables(const std::string &encoded_name, const Ioss::Field &field,
                                      const std::string &entity_type,
                                      const std::string &field_name) const
  {
    adios_wrapper.PutMetaVariable<int>(Role_meta, field.get_role(), encoded_name);
    bool        transformed_field = use_transformed_storage(field, entity_type, field_name);
    std::string var_type =
        transformed_field ? field.transformed_storage()->name() : field.raw_storage()->name();
    adios_wrapper.PutMetaVariable<std::string>(Var_type_meta, var_type, encoded_name);
  }

  template <typename T>
  int64_t DatabaseIO::put_field_internal_t(T entity, const Ioss::Field &field, void *data,
                                           size_t data_size) const
  {
    if (!data || !data_size) {
      return 0;
    }
    std::string        entity_type = entity->type_string();
    const std::string &field_name  = field.get_name();
    if (find_field_in_mapset(entity_type, field_name, Ignore_fields)) {
      return 0;
    }

    int num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      std::string encoded_name = encode_field_name({entity_type, entity->name(), field_name});
      switch (field.get_type()) {
      case Ioss::Field::BasicType::DOUBLE: put_data<double>(data, encoded_name); break;
      case Ioss::Field::BasicType::INT32: put_data<int32_t>(data, encoded_name); break;
      case Ioss::Field::BasicType::INT64: put_data<int64_t>(data, encoded_name); break;
      case Ioss::Field::BasicType::COMPLEX: put_data<Complex>(data, encoded_name); break;
      case Ioss::Field::BasicType::CHARACTER: put_data<char>(data, encoded_name); break;
      default: {
        std::ostringstream errmsg;
        errmsg << "INTERNAL ERROR: Invalid field type. "
               << "Something is wrong in the Ioad::DatabaseIO::put_field_internal_t() function. "
               << "Please report.\n";
        IOSS_ERROR(errmsg);
      }
      }
      put_meta_variables(encoded_name, field, entity_type, field_name);
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(eb, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(eb, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(fb, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(sb, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(ns, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(es, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(fs, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(es, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal_t(ss, field, data, data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    // TODO: Make sure that Commset should be handled like a "set" and not like a "block".
    return put_field_internal_t(cs, field, data, data_size);
  }

  template <typename T>
  DatabaseIO::FieldInfoType DatabaseIO::get_expected_variable_infos_from_map(
      const EntityMapType &entity_map, const std::string &entity_type,
      const std::string &block_name, const std::string &var_name) const
  {
    if (entity_map.find(block_name) == entity_map.end()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: block name " << block_name << " not found.\n";
      IOSS_ERROR(errmsg);
    }

    if (entity_map.at(block_name).find(var_name) == entity_map.at(block_name).end()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: block name not found or does not contain `" << var_name << "` field.\n";
      IOSS_ERROR(errmsg);
    }

    if (entity_map.at(block_name).at(var_name).second != GetType<T>()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: `" << var_name << "` field type does not match expected type.\n";
      IOSS_ERROR(errmsg);
    }
    std::string variable = encode_field_name({entity_type, block_name, var_name});
    return get_variable_infos<T>(variable);
  }

  DatabaseIO::FieldInfoType DatabaseIO::get_variable_infos_from_map(
      const EntityMapType &entity_map, const std::string &entity_type,
      const std::string &block_name, const std::string &var_name) const
  {
    // No check to verify that the block_name, entity_type, variable_name, and type exist as they
    // have been extracted earlier from the file and everything should still be up to date.
    std::string variable = encode_field_name({entity_type, block_name, var_name});
    std::string type     = entity_map.at(block_name).at(var_name).second;
    // Simply to start the "else if" section with "if".
    if (type == "not supported") {}
#define declare_template_instantiation(T)                                                          \
  else if (type == GetType<T>())                                                                   \
  {                                                                                                \
    return get_variable_infos<T>(variable);                                                        \
  }
    ADIOS2_FOREACH_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    else
    {
      std::ostringstream errmsg;
      errmsg << "INTERNAL ERROR: Invalid variable type. "
             << "Something is wrong in the "
                "Ioad::DatabaseIO::get_variable_infos_from_map() function."
             << "Please report.\n";
      IOSS_ERROR(errmsg);
    }
    return FieldInfoType();
  }

  template <>
  int64_t DatabaseIO::get_entities<Ioss::NodeBlock>(const FieldsMapType &fields_map,
                                                    const FieldsMapType &properties_map)
  {
    const std::string block_type = get_entity_type<Ioss::NodeBlock>();
    // For exodusII, there is only a single node block which contains
    // all of the nodes.
    // The default id assigned is '1' and the name is 'nodeblock_1'
    const EntityMapType &entity_map = fields_map.at(block_type);
    std::string          block_name = "nodeblock_1";

    // `mesh_model_coordinates` field is automatically created in NodeBlock constructor.
    std::string coord_var_name = "mesh_model_coordinates";
    // Get `node_count` and `spatialDimension`.
    FieldInfoType model_coordinates_infos = get_expected_variable_infos_from_map<double>(
        entity_map, block_type, block_name, coord_var_name);
    spatialDimension = model_coordinates_infos.component_count;
    if (!spatialDimension) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Variable `" << coord_var_name
             << "` in BP file without any dimension information.\n";
      IOSS_ERROR(errmsg);
    }
    auto block = new Ioss::NodeBlock(this, block_name, model_coordinates_infos.node_boundaries_size,
                                     spatialDimension);
    block->property_add(Ioss::Property("id", 1));
    block->property_add(Ioss::Property("guid", util().generate_guid(1)));
    for (auto &variable : entity_map.at(block_name)) {
      // Since some fields are created automatically, we need to avoid recreating them when loading
      // the file.
      if (!block->field_exists(variable.first)) {
        FieldInfoType variable_infos =
            get_variable_infos_from_map(entity_map, block_type, block_name, variable.first);
        Ioss::Field field(variable.first, variable_infos.basic_type, variable_infos.variable_type,
                          variable_infos.role, variable_infos.node_boundaries_size);
        block->field_add(field);
      }
    }
    add_entity_properties(block, properties_map);
    bool added = get_region()->add(block);
    if (!added) {
      delete block;
    }
    return model_coordinates_infos.node_boundaries_size;
  }

  template <typename T>
  void DatabaseIO::add_entity_property(Ioss::GroupingEntity *ge, const std::string &encoded_name,
                                       const std::string &property_name)
  {
    T variable;
    // Meta variable read with `Get()` instead of `GetMetaVariable()` because name is already
    // encoded. `GetMetaVariable()` simply encodes the variable name as a meta variable before
    // performing a `Get()` call.
    adios_wrapper.GetSync<T>(encoded_name, variable);
    Ioss::Property property(property_name, variable);
    ge->property_add(property);
  }

  void DatabaseIO::add_entity_properties(Ioss::GroupingEntity *ge,
                                         const FieldsMapType &properties_map, std::string name)
  {
    std::string entity_type = ge->type_string();
    if (properties_map.find(entity_type) == properties_map.end()) {
      return;
    }
    const EntityMapType &entity_property_map = properties_map.at(entity_type);
    GlobalMapType        properties_info;
    if (!name.empty()) {
      // Verify that the given name exists.
      if (entity_property_map.find(name) == entity_property_map.end()) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Property map doesn't contain " << name << " inside " << entity_type
               << ".\n";
        IOSS_ERROR(errmsg);
      }
      properties_info = entity_property_map.at(name);
    }
    else {
      std::string entity_name = ge->name();
      if (entity_property_map.find(entity_name) == entity_property_map.end()) {
        return;
      }
      properties_info = entity_property_map.at(entity_name);
    }
    for (const auto &property_info : properties_info) {
      std::string property_name = property_info.first;
      std::string encoded_name  = property_info.second.first;
      std::string type          = property_info.second.second;
      if (type == GetType<std::string>()) {
        add_entity_property<std::string>(ge, encoded_name, property_name);
      }
      else if (type == GetType<double>()) {
        add_entity_property<double>(ge, encoded_name, property_name);
      }
      else if (type == GetType<int>() || type == "int32_t") {
        add_entity_property<int>(ge, encoded_name, property_name);
      }
      else {
        std::ostringstream errmsg;
        errmsg << "ERROR: Property '" << property_name << "' with type `" << type
               << "` not supported.\n";
        IOSS_ERROR(errmsg);
      }
    }
  }

  template <typename T>
  std::string DatabaseIO::get_property_value(const FieldsMapType &properties_map,
                                             const std::string   &entity_type,
                                             const std::string   &entity_name,
                                             const std::string   &property_name) const
  {
    T property_value;
    if (properties_map.find(entity_type) == properties_map.end()) {
      return property_value;
    }
    const EntityMapType &entity_property_map = properties_map.at(entity_type);
    if (entity_property_map.find(entity_name) == entity_property_map.end()) {
      return property_value;
    }
    const GlobalMapType properties_info = entity_property_map.at(entity_name);
    const auto         &property_info   = properties_info.find(property_name);
    if (property_info != properties_info.end()) {
      const std::string encoded_name  = property_info->second.first;
      const std::string type          = property_info->second.second;
      const std::string expected_type = GetType<T>();
      if (type != expected_type) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Property type `" << type << "` doesn't match expected type `"
               << expected_type << "`.\n";
        IOSS_ERROR(errmsg);
      }
      // Meta variable read with `Get()` instead of `GetMetaVariable()` because name is already
      // encoded. `GetMetaVariable()` simply encodes the variable name as a meta variable before
      // performing a `Get()` call.
      adios_wrapper.GetSync<T>(encoded_name, property_value);
    }
    return property_value;
  }

  template <>
  int64_t DatabaseIO::get_entities<Ioss::SideSet>(const FieldsMapType &fields_map,
                                                  const FieldsMapType &properties_map)
  {
    int64_t     count       = 0;
    std::string entity_type = get_entity_type<Ioss::SideSet>();
    if (fields_map.find(entity_type) == fields_map.end()) {
      return count;
    }
    std::string sideblock_type = get_entity_type<Ioss::SideBlock>();
    // Do SideBlocks exist in the fields_map. This will be used when loading the SideBlock
    // for each SideSet, but we check this here outside of a loop as it only needs to be done once.
    bool no_sideblocks = (fields_map.find(sideblock_type) == fields_map.end());

    const EntityMapType &sidesets_map   = fields_map.at(entity_type);
    const EntityMapType &sideblocks_map = fields_map.at(sideblock_type);
    for (const auto &entity : sidesets_map) {
      std::string    entity_name = entity.first;
      Ioss::SideSet *ss          = new Ioss::SideSet(this, entity_name);
      bool           ss_added    = get_region()->add(ss);
      add_entity_properties(ss, properties_map);
      if (!ss_added) {
        delete ss;
        return count;
      }
      for (const auto &variable_pair : entity.second) {
        // Since some fields are created automatically, we need to avoid recreating them when
        // loading the file.
        // First, check that field is actually a field and not a list of sideblocks
        if (!is_sideblock_name(variable_pair.first)) {
          if (!ss->field_exists(variable_pair.first)) {
            FieldInfoType infos = get_variable_infos_from_map(sidesets_map, entity_type,
                                                              entity_name, variable_pair.first);

            Ioss::Field field(variable_pair.first, infos.basic_type, infos.variable_type,
                              infos.role, infos.node_boundaries_size);
            count += infos.node_boundaries_size;
            ss->field_add(field);
          }
        }
        else {
          // The field in the bp file is actually a list of sideblocks
          // If no sideblock, don't worry about this section and move on to the next
          // sideset.
          if (no_sideblocks) {
            continue;
          }
          adios2::Variable<std::string> sideblock_var =
              adios_wrapper.InquireVariable<std::string>(variable_pair.second.first);

          if (sideblock_var) {
            std::string stringified_names;
            adios_wrapper.GetSync<std::string>(sideblock_var, stringified_names);
            std::vector<std::string> block_names =
                Ioss::tokenize(stringified_names, Sideblock_separator);
            for (std::string block_name : block_names) {
              if (sideblocks_map.find(block_name) != sideblocks_map.end()) {
                bool first = true;
                for (const auto &sideblock_field_pair : sideblocks_map.at(block_name)) {
                  std::string   field_name  = sideblock_field_pair.first;
                  FieldInfoType block_infos = get_variable_infos_from_map(
                      sideblocks_map, sideblock_type, block_name, field_name);
                  // Create sideblock first
                  if (first) {
                    Ioss::SideBlock *sideblock = new Ioss::SideBlock(
                        this, block_name, block_infos.topology, block_infos.parent_topology,
                        block_infos.node_boundaries_size);
                    bool block_added = ss->add(sideblock);
                    if (!block_added) {
                      delete sideblock;
                      std::ostringstream errmsg;
                      errmsg << "ERROR: Could not add sideblock `" << block_name << "` to sideset `"
                             << ss->name() << "`.\n";
                      IOSS_ERROR(errmsg);
                    }
                    first = false;
                  }
                  Ioss::SideBlock *sideblock = ss->get_side_block(block_name);
                  add_entity_properties(sideblock, properties_map);
                  // Add fields to sideblock if it was not automatically added.
                  if (!sideblock->field_exists(sideblock_field_pair.first)) {
                    Ioss::Field field(field_name, block_infos.basic_type, block_infos.variable_type,
                                      block_infos.role, block_infos.node_boundaries_size);
                    sideblock->field_add(field);
                  }
                }
              }
              else {
                // throw an error if sideblock not found? Check if sideblock name matches sideset
                // name: it could
                // have been added automatically at read time if the list was empty (see STK).
              }
            }
          }
        }
      }
    }
    return count;
  }

  template <typename T>
  int64_t DatabaseIO::get_entities(const FieldsMapType &fields_map,
                                   const FieldsMapType &properties_map)
  {
    int64_t     count       = 0;
    std::string entity_type = get_entity_type<T>();
    if (fields_map.find(entity_type) == fields_map.end()) {
      return count;
    }
    const EntityMapType &entity_map = fields_map.at(entity_type);
    for (const auto &variable_pair : entity_map) {
      std::string entity_name = variable_pair.first;
      // Get size and type info for the entity using the first element in the map.
      FieldInfoType infos_to_create_entity = get_variable_infos_from_map(
          entity_map, entity_type, entity_name, variable_pair.second.begin()->first);
      std::string block_type =
          get_property_value<std::string>(properties_map, entity_type, entity_name, "entity_type");
      block_type = block_type.empty() ? infos_to_create_entity.topology : block_type;
      T *entity =
          NewEntity<T>(this, entity_name, block_type, infos_to_create_entity.node_boundaries_size);
      count += infos_to_create_entity.node_boundaries_size;
      add_entity_properties(entity, properties_map);

      bool added = get_region()->add(entity);
      if (!added) {
        delete entity;
      }
      for (const auto &field_pair : variable_pair.second) {
        // Since some fields are created automatically, we need to avoid recreating them when
        // loading the file.
        // Note: We get the information about the first field twice: once before this loop, and
        // once inside the loop. The code could be modified to perform this action only once but
        // most likely the code will be more complex and will not save a lot of computation time.
        if (!entity->field_exists(field_pair.first)) {
          FieldInfoType infos =
              get_variable_infos_from_map(entity_map, entity_type, entity_name, field_pair.first);
          Ioss::Field field(field_pair.first, infos.basic_type, infos.variable_type, infos.role,
                            infos.node_boundaries_size);
          entity->field_add(field);
        }
      }
    }
    return count;
  }

  // TODO: Consolidate this code with generice `get_entities()`  function as there is very little
  // difference between the two.
  template <>
  int64_t DatabaseIO::get_entities<Ioss::CommSet>(const FieldsMapType &fields_map,
                                                  const FieldsMapType &properties_map)
  {
    int64_t     count       = 0;
    std::string entity_type = get_entity_type<Ioss::CommSet>();
    if (fields_map.find(entity_type) == fields_map.end()) {
      return count;
    }
    const EntityMapType &entity_map = fields_map.at(entity_type);
    for (const auto &variable_pair : entity_map) {
      std::string entity_name = variable_pair.first;
      // Get size and type info for the entity using the first element in the map.
      FieldInfoType infos_to_create_entity = get_variable_infos_from_map(
          entity_map, entity_type, entity_name, variable_pair.second.begin()->first);
      std::string block_type =
          get_property_value<std::string>(properties_map, entity_type, entity_name, "entity_type");
      std::string commset_name = "commset_" + block_type;

      block_type            = block_type.empty() ? infos_to_create_entity.topology : block_type;
      Ioss::CommSet *entity = NewEntity<Ioss::CommSet>(this, commset_name, block_type,
                                                       infos_to_create_entity.node_boundaries_size);
      count += infos_to_create_entity.node_boundaries_size;
      add_entity_properties(entity, properties_map);
      // Save original name as entity property
      Ioss::Property property(original_name, entity_name);
      entity->property_add(property);

      bool added = get_region()->add(entity);
      if (!added) {
        delete entity;
      }
      for (const auto &field_pair : variable_pair.second) {
        // Since some fields are created automatically, we need to avoid recreating them when
        // loading the file.
        // Note: We get the information about the first field twice: once before this loop, and
        // once inside the loop. The code could be modified to perform this action only once but
        // most likely the code will be more complex and will not save a lot of computation time.
        if (!entity->field_exists(field_pair.first)) {
          FieldInfoType infos =
              get_variable_infos_from_map(entity_map, entity_type, entity_name, field_pair.first);
          Ioss::Field field(field_pair.first, infos.basic_type, infos.variable_type, infos.role,
                            infos.node_boundaries_size);
          entity->field_add(field);
        }
      }
    }
    return count;
  }

  std::string DatabaseIO::get_optional_string_variable(const std::string &field_name,
                                                       const std::string &string_variable) const
  {
    auto        tokens = Ioss::tokenize(field_name, Name_separator);
    std::string entity = encode_field_name({tokens[0], tokens[1]});
    auto        v      = adios_wrapper.InquireVariable<std::string>(
        adios_wrapper.EncodeMetaVariable(string_variable, entity));
    if (v) {
      return adios_wrapper.GetMetaVariable<std::string>(string_variable, entity);
    }
    else {
      return std::string{};
    }
  }

  template <typename T>
  DatabaseIO::BlockInfoType DatabaseIO::get_block_infos(const adios2::Variable<T> &var) const
  {
    std::map<size_t, std::vector<typename adios2::Variable<T>::Info>> allblocks =
        adios_wrapper.AllStepsBlocksInfo(var);
    if (allblocks.empty()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Empty BP variable\n";
      IOSS_ERROR(errmsg);
    }
    BlockInfoType infos;
    // For each time step.
    for (const auto &blockpair : allblocks) {
      const std::vector<typename adios2::Variable<T>::Info> &blocks = blockpair.second;
      // Find in vector if this variable is defined for the current rank process. This means
      // that there is one block for which the rank encoded as the first value in the `Start` array
      // matches the current rank.
      // Note: one block per rank.
      for (const auto &block : blocks) {
        if (block.Start[0] != rank) {
          infos.global_size += block.Count[1];
          // This is not the block corresponding to the current process (rank).
          continue;
        }
        if (infos.Count.empty()) {
          infos.Count = block.Count;
        }
        else if (infos.Count != block.Count) {
          std::ostringstream errmsg;
          errmsg
              << "ERROR: Variable changes sizes over steps. Not supported by Ioad::DatabaseIO.\n";
          IOSS_ERROR(errmsg);
        }
        // }
        // Steps: will save last step
        infos.steps.push_back(blockpair.first);
        // We found the block we were searching for, no need to continue looping over all
        // the other blocks.
        break;
      }
    }
    return infos;
  }

  template <typename T>
  DatabaseIO::FieldInfoType DatabaseIO::get_variable_infos(const std::string &var_name) const
  {
    FieldInfoType       infos;
    adios2::Variable<T> var = adios_wrapper.InquireVariable<T>(var_name);

    // Dimension is expected to be 3 because:
    // 0 = rank
    // 1 = size of field
    // 2 = number of components
    size_t ndim = var.Shape().size();
    if (ndim != 3) {
      std::ostringstream errmsg;
      errmsg << "ERROR: BP variable dimension should be 3.\n";
      IOSS_ERROR(errmsg);
    }
    BlockInfoType block_infos = get_block_infos<T>(var);
    size_t        laststep    = 0;
    if (!block_infos.steps.empty() && !block_infos.Count.empty()) {
      laststep = block_infos.steps.back();
      // Skipping dimension=0 which is equal to the rank. Dimension=1 encodes beginning and count
      // of node count for this variable. Dimension=2 encodes the number of components for this
      // variables, and should always start at 0;
      infos.component_count      = block_infos.Count[2];
      infos.node_boundaries_size = block_infos.Count[1];
      infos.steps                = block_infos.steps;
    }
    // Get topology: Only relevant for blocks, not for sets.
    infos.topology = get_optional_string_variable(var_name, Topology_meta);
    // Get parent topology: Only relevant for sideblocks.
    infos.parent_topology = get_optional_string_variable(var_name, Parent_topology_meta);
    // Get VariableType
    infos.variable_type = adios_wrapper.GetMetaVariable<std::string>(Var_type_meta, var_name);
    infos.basic_type    = template_to_basic_type<T>();
    // Only transient fields can have steps.
    infos.role =
        static_cast<Ioss::Field::RoleType>(adios_wrapper.GetMetaVariable<int>(Role_meta, var_name));
    if (infos.role != Ioss::Field::RoleType::TRANSIENT && laststep != 0) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Last step should be 0 for non-transient fields. "
             << "Something is wrong in the Ioad::DatabaseIO::get_variable_size() function.\n";
      IOSS_ERROR(errmsg);
    }
    return infos;
  }

  void DatabaseIO::get_globals(const GlobalMapType &globals_map,
                               const FieldsMapType &properties_map)
  {
    // Check "time" attribute and global variable.
    timeScaleFactor = adios_wrapper.GetAttribute<double>(Time_scale_factor, true, 1.0);
    Ioss::SerializeIO serializeIO_(this);
    Ioss::Region     *this_region = get_region();
    if (globals_map.find(Time_meta) != globals_map.end()) {
      // Load time steps
      // 1) Check that the time type is `double` as expected.
      adios2::Variable<double> time_var = adios_wrapper.InquireVariable<double>(Time_meta);
      if (time_var) {
        std::vector<double> tsteps(time_var.Steps());
        if (!is_streaming) {
          time_var.SetStepSelection(std::make_pair(
              time_var.StepsStart(), time_var.Steps())); // Doesn't work with streaming.
        }
        else {
          // For streaming, we probably only want to read the current time as we do not want to
          // use SetStepSelection that would prohibit the usage of `BeginStep()/EndStep()`.
          // std::ostringstream errmsg;
          // errmsg << "ERROR: Streaming is not yet supported for reading.\n";
          // IOSS_ERROR(errmsg);
        }
        adios_wrapper.GetSync(time_var, Data(tsteps));
        for (size_t step = 0; step < time_var.Steps(); step++) {
          // if (tsteps[i] <= last_time) { TODO: Check last time step before writing everything
          this_region->add_state_nl(tsteps[step] * timeScaleFactor);
        }
      }
      else {
        std::ostringstream errmsg;
        errmsg << "ERROR: Timestep global detected in file but cannot read it.\n";
        IOSS_ERROR(errmsg);
      }
    }
    // Add region properties
    add_entity_properties(this_region, properties_map, region_name);
  }

  void DatabaseIO::read_meta_data_nl()
  {
    check_processor_info();
    Ioss::Region *region = get_region();

    // Define properties
    if (is_streaming) {
      region->property_update("streaming", 1);
    }
    // Only get schema version attribute as it is the only one we expect.
    adios_wrapper.GetAttribute<unsigned int>(Schema_version_string);

    // Get all variables
    GlobalMapType globals_map;
    // Entity_type/Entity_name::Property -> name & Type
    FieldsMapType properties_map;
    // Entity_type/Entity_name/Field -> name & Type
    FieldsMapType fields_map;

    const std::map<std::string, std::map<std::string, std::string>> variables =
        adios_wrapper.AvailableVariables();
    for (const auto &vpair : variables) {
      const std::string &name           = vpair.first;
      auto               name_type_pair = std::make_pair(name, vpair.second.at("Type"));
      auto               tokens         = Ioss::tokenize(name, Name_separator);
      switch (tokens.size()) {
      case 1: globals_map[tokens[0]] = name_type_pair; break;
      case 2: {
        // Only save properties.
        std::string meta_var         = "";
        std::string property         = "";
        std::tie(meta_var, property) = adios_wrapper.DecodeMetaName(tokens[1]);
        if (!property.compare(0, property_meta.size(), property_meta)) {
          properties_map[tokens[0]][meta_var][property.substr(property_meta.size())] =
              name_type_pair;
        }
        break;
      }
      case 3: {
        std::string is_meta            = "";
        std::tie(std::ignore, is_meta) = adios_wrapper.DecodeMetaName(tokens[2]);
        if (is_meta.empty()) {
          fields_map[tokens[0]][tokens[1]][tokens[2]] = name_type_pair;
        }
        break;
      }
      default: {
        std::ostringstream errmsg;
        errmsg << "ERROR: Invalid encoded entity name.\n";
        IOSS_ERROR(errmsg);
      }
      }
    }

    get_globals(globals_map, properties_map);
    nodeCount    = get_entities<Ioss::NodeBlock>(fields_map, properties_map);
    edgeCount    = get_entities<Ioss::EdgeBlock>(fields_map, properties_map);
    faceCount    = get_entities<Ioss::FaceBlock>(fields_map, properties_map);
    elementCount = get_entities<Ioss::ElementBlock>(fields_map, properties_map);
    check_side_topology();
    get_entities<Ioss::SideSet>(fields_map, properties_map);
    get_entities<Ioss::NodeSet>(fields_map, properties_map);
    get_entities<Ioss::EdgeSet>(fields_map, properties_map);
    get_entities<Ioss::FaceSet>(fields_map, properties_map);
    get_entities<Ioss::ElementSet>(fields_map, properties_map);
    get_entities<Ioss::CommSet>(fields_map, properties_map);

    read_region(fields_map);
    read_communication_metadata();
  }

  void DatabaseIO::read_region(const FieldsMapType &fields_map)
  {
    // // Add properties and fields to the 'owning' region.
    // // Also defines member variables of this class...

    if (nodeCount == 0) {
      Ioss::WarnOut() << "No nodes were found in the model, file '" << decoded_filename() << "'\n";
    }
    else if (nodeCount < 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "ERROR: Negative node count was found in the model\n"
             << "       File: '" << decoded_filename() << "'.\n";
      IOSS_ERROR(errmsg);
    }

    if (elementCount == 0) {
      Ioss::WarnOut() << "No elements were found in the model, file: '" << decoded_filename()
                      << "'\n";
    }

    if (elementCount < 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "ERROR: Negative element count was found in the model, file: '"
             << decoded_filename() << "'";
      IOSS_ERROR(errmsg);
    }
    Ioss::Region *region = get_region();

    // See if any coordinate frames exist on mesh.  If so, define them on region.
    // add_coordinate_frames();
    if (fields_map.find(coordinate_frame_name) != fields_map.end()) {
      const EntityMapType &entity_map = fields_map.at(coordinate_frame_name);
      for (const auto &variable_pair : entity_map) {
        std::vector<double> coord(9, 0);

        int64_t     id       = std::stoll(variable_pair.first);
        char        tag      = variable_pair.second.begin()->first[0];
        std::string var_name = variable_pair.second.begin()->second.first;
        adios_wrapper.GetSync(var_name, Data(coord));
        Ioss::CoordinateFrame coordinate_frame(id, tag, Data(coord));
        region->add(coordinate_frame);
      }
    }
  }

  void DatabaseIO::read_communication_metadata()
  {
    // Number of processors file was decomposed for
    int num_proc = adios_wrapper.GetAttribute<unsigned long>(Processor_number_meta);

    // Get global data (over all processors)
    int64_t global_nodes    = nodeCount;
    int64_t global_elements = elementCount;
    int64_t global_eblocks  = 0; // unused
    int64_t global_nsets    = 0; // unused
    int64_t global_ssets    = 0; // unused

    int64_t num_internal_nodes = nodeCount;
    int64_t num_border_nodes   = 0;
    int64_t num_internal_elems = elementCount;
    int64_t num_border_elems   = 0;

    Ioss::Region *region = get_region();

    region->property_add(Ioss::Property("processor_count", num_proc));

    region->property_add(Ioss::Property("internal_node_count", num_internal_nodes));
    region->property_add(Ioss::Property("border_node_count", num_border_nodes));
    region->property_add(Ioss::Property("internal_element_count", num_internal_elems));
    region->property_add(Ioss::Property("border_element_count", num_border_elems));
    region->property_add(Ioss::Property("global_node_count", global_nodes));
    region->property_add(Ioss::Property("global_element_count", global_elements));
    region->property_add(Ioss::Property("global_element_block_count", global_eblocks));
    region->property_add(Ioss::Property("global_node_set_count", global_nsets));
    region->property_add(Ioss::Property("global_side_set_count", global_ssets));

    // Possibly, the following 4 fields should be nodesets and element
    // sets instead of fields on the region...
    region->field_add(Ioss::Field("internal_nodes", region->field_int_type(), IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_internal_nodes));
    region->field_add(Ioss::Field("border_nodes", region->field_int_type(), IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_border_nodes));
    region->field_add(Ioss::Field("internal_elements", region->field_int_type(), IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_internal_elems));
    region->field_add(Ioss::Field("border_elements", region->field_int_type(), IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_border_elems));

    assert(nodeCount == num_internal_nodes + num_border_nodes);
    assert(elementCount == num_internal_elems + num_border_elems);
  }
  void DatabaseIO::check_processor_info()
  {
    std::ostringstream errmsg;
    // Get Processor information
    unsigned long number_proc_read =
        adios_wrapper.GetAttribute<unsigned long>(Processor_number_meta);

    if (number_proc < number_proc_read) {
      errmsg << "Processor decomposition count in file (" << number_proc_read
             << ") is larger than current processor count (" << number_proc
             << "). Configuration not supported.\n";
      IOSS_ERROR(errmsg);
    }
    else if (number_proc > number_proc_read) {
      Ioss::WarnOut() << "This file was originally written on " << number_proc_read
                      << " processors, but is now being read using " << number_proc
                      << " processors.\n";
    }
    if (rank < number_proc_read) {
      // Only get info for processors that actually have an id.
      unsigned long processor_id;
      get_data<unsigned long>(static_cast<void *>(&processor_id), Processor_id_meta);
      if (rank != processor_id) {
        Ioss::WarnOut() << "This file was originally written on processor " << processor_id
                        << ", but is now being read on processor " << rank
                        << ". This may cause problems if there is any processor-dependent data on "
                           "the file.\n";
      }
    }
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(reg, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      if (field.get_name() == "mesh_model_coordinates_x" ||
          field.get_name() == "mesh_model_coordinates_y" ||
          field.get_name() == "mesh_model_coordinates_z") {
        Ioss::Field         coord_field = nb->get_field("mesh_model_coordinates");
        std::vector<double> coord(num_to_get * spatialDimension);
        get_field_internal_t(nb, coord_field, Data(coord), data_size * spatialDimension);

        // Cast 'data' to correct size -- double
        double *rdata  = static_cast<double *>(data);
        int     offset = 0;
        if (field.get_name() == "mesh_model_coordinates_x") {
          offset = 0;
        }
        else if (field.get_name() == "mesh_model_coordinates_y") {
          offset = 1;
        }
        else if (field.get_name() == "mesh_model_coordinates_z") {
          offset = 2;
        }
        for (size_t i = 0; i < num_to_get; i++) {
          rdata[i] = coord[i * spatialDimension + offset];
        }
      }
      else {
        return get_field_internal_t(nb, field, data, data_size);
      }
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(eb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(fb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(eb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(sb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(ns, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(es, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(fs, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(es, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(ss, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal_t(cs, field, data, data_size);
  }

  int64_t DatabaseIO::get_field_internal_t(const Ioss::GroupingEntity *entity,
                                           const Ioss::Field &field, void *data,
                                           size_t data_size) const
  {
    if (!data || !data_size) {
      return 0;
    }
    int num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      const std::string entity_type = entity->type_string();

      // Check if field name has changed. Rely on property `original_name` if
      // it exists.
      const std::string  entity_name = entity->property_exists(original_name)
                                           ? entity->get_property(original_name).get_string()
                                           : entity->name();
      const std::string &field_name  = field.get_name();

      if (find_field_in_mapset(entity_type, field_name, Ignore_fields)) {
        return num_to_get;
      }

      std::string encoded_name       = encode_field_name({entity_type, entity_name, field_name});
      bool        use_step_selection = false;
      if (field.get_role() == Ioss::Field::RoleType::TRANSIENT && !is_streaming) {
        use_step_selection = true;
      }
      switch (field.get_type()) {
      case Ioss::Field::BasicType::DOUBLE:
        get_data<double>(data, encoded_name, use_step_selection);
        break;
      case Ioss::Field::BasicType::INT32:
        get_data<int32_t>(data, encoded_name, use_step_selection);
        break;
      case Ioss::Field::BasicType::INT64:
        get_data<int64_t>(data, encoded_name, use_step_selection);
        break;
      case Ioss::Field::BasicType::COMPLEX:
        get_data<Complex>(data, encoded_name, use_step_selection);
        break;
      case Ioss::Field::BasicType::CHARACTER:
        get_data<char>(data, encoded_name, use_step_selection);
        break;
      default:
        std::ostringstream errmsg;
        errmsg << "INTERNAL ERROR: Invalid field type. "
               << "Something is wrong in the Ioad::DatabaseIO::get_field_internal_t() function. "
               << "Please report.\n";
        IOSS_ERROR(errmsg);
      }
    }
    return num_to_get;
  }

  template <typename T>
  void DatabaseIO::get_data(void *data, const std::string &encoded_name,
                            bool use_step_selection) const
  {
    adios2::Variable<T> entities = adios_wrapper.InquireVariable<T>(encoded_name);
    if (entities) {
      T           *rdata  = static_cast<T *>(data);
      adios2::Dims size   = entities.Shape();
      size[0]             = 1;
      adios2::Dims offset = entities.Start();
      offset[0]           = rank;

      entities.SetSelection(adios2::Box<adios2::Dims>(offset, size));
      // if transient, set step that should be read if not streaming.
      if (use_step_selection) {
        size_t step = get_current_state();
        entities.SetStepSelection(std::make_pair(step, 1));
      }
      else {
        // if streaming, we are reading step selected by `BeginStep()/EndStep()`.
      }

      // TODO: Set selection per rank. Support written by N processes, and loaded by M processes.
      adios_wrapper.GetSync<T>(entities, rdata);
    }
    else {
      if (!is_streaming) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Required `" << encoded_name << "` variable not found in file.\n";
        IOSS_ERROR(errmsg);
      }
      else {
        Ioss::WarnOut() << "The variable `" << encoded_name << "` was not found.\n";
      }
    }
  }

  void DatabaseIO::compute_block_membership_nl(Ioss::SideBlock          *efblock,
                                               std::vector<std::string> &block_membership) const
  {
    const Ioss::ElementBlockContainer &element_blocks = get_region()->get_element_blocks();
    assert(Ioss::Utils::check_block_order(element_blocks));

    Ioss::Int64Vector block_ids(element_blocks.size());
    if (block_ids.size() == 1) {
      block_ids[0] = 1;
    }
    else {
      Ioss::Int64Vector element_side;
      efblock->get_field_data("element_side", element_side);

      size_t              number_sides = element_side.size() / 2;
      Ioss::ElementBlock *block        = nullptr;
      for (size_t iel = 0; iel < number_sides; iel++) {
        int64_t elem_id = element_side[2 * iel]; // Vector contains both element and side.
        // elem_id         = elemMap.global_to_local(elem_id);
        if (block == nullptr || !block->contains(elem_id)) {
          block = get_region()->get_element_block(elem_id);
          assert(block != nullptr);
          size_t block_order = block->get_property("original_block_order").get_int();
          assert(block_order < block_ids.size());
          block_ids[block_order] = 1;
        }
      }
    }

    // Synchronize among all processors....
    if (isParallel) {
      util().global_array_minmax(block_ids, Ioss::ParallelUtils::DO_MAX);
    }

    for (const auto block : element_blocks) {
      size_t block_order = block->get_property("original_block_order").get_int();
      assert(block_order < block_ids.size());
      if (block_ids[block_order] == 1) {
        if (!Ioss::Utils::block_is_omitted(block)) {
          block_membership.push_back(block->name());
        }
      }
    }
  }

  int DatabaseIO::get_current_state() const
  {
    // value returned is 1-based, whereas ADIOS expect 0-based values
    int step = get_region()->get_current_state() - 1;

    if (step < 0) {
      std::ostringstream errmsg;
      errmsg << "ERROR: No currently active state.  The calling code must call "
                "Ioss::Region::begin_state(int step)\n"
             << "       to set the database timestep from which to read the transient data.\n"
             << "       [" << get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    }
    return step;
  }

} // namespace Ioad
