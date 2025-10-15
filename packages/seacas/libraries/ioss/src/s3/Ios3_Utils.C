// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "Ios3_Utils.h"

#include "Ioss_CommSet.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h" // for Region
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h" // for Region
#include "Ioss_Region.h"
#include "Ioss_Region.h" // for Region
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_State.h" // for State
#include "Ioss_StructuredBlock.h"

#include "cereal/archives/portable_binary.hpp"
#include "cereal/types/array.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/utility.hpp"
#include "cereal/types/vector.hpp"
#include <sstream>

namespace Ios3 {

  state_entry_t::state_entry_t(const Ioss::Region &r)
      : count(r.get_property("state_count").get_int()),
        value{0, count * sizeof(state_entry_t::basic_type)}
  {
  }

  sideblock_entry_t::sideblock_entry_t(const Ioss::SideBlock &sb) : entity_count(sb.entity_count())
  {
  }

  /* The constructor for SideBlocks requires attributes that are
     neither properties nor fields.  This function collects and
     serializes these attributes so that the SideBlock can be
     reconstructed from the data stored in S3. */
  PackedBytes pack_sideblock(const Ioss::SideBlock &sb)
  {
    // meta, entry, data, data_size
    auto v = PackedBytes(sizeof(meta_entry_t) + sizeof(sideblock_entry_t));

    auto meta          = reinterpret_cast<meta_entry_t *>(v.data());
    meta->ioss_type    = meta_entry_t::IossType::Ios3SideBlock;
    meta->value.offset = 0;
    meta->value.size   = sizeof(sideblock_entry_t);

    auto entry = reinterpret_cast<sideblock_entry_t *>(
        reinterpret_cast<void *>(v.data() + sizeof(meta_entry_t) + meta->value.offset));

    entry->entity_count = sb.entity_count();

    return v;
  }

  int64_t unpack_sideblocks(PackedBytes &v)
  {
    auto meta = reinterpret_cast<meta_entry_t *>(v.data());
    assert(meta->ioss_type == meta_entry_t::IossType::Ios3SideBlock);

    auto entry =
        reinterpret_cast<sideblock_entry_t *>(v.data() + sizeof(meta_entry_t) + meta->value.offset);

    return entry->entity_count;
  }

  /* The constructor for StructuredBlocks requires attributes that are
     neither properties nor fields. This function collects and
     serializes these attributes so that the StructuredBlock can be
     reconstructed from the data stored in S3. */
  PackedBytes pack_structuredblock(const Ioss::StructuredBlock &sb)
  {
    std::ostringstream                  outstream;
    cereal::PortableBinaryOutputArchive oarchive(outstream);
    oarchive(sb);
    std::streamoff length = outstream.tellp();

    // meta, entry, data, data_size
    auto        v         = PackedBytes(length);
    std::string outstring = outstream.str();
    outstring.copy(reinterpret_cast<char *>(v.data()), length);

    return v;
  }

  void unpack_structuredblock(PackedBytes &v, Ioss::StructuredBlock &sb)
  {
    std::stringstream instream;
    char             *p = (char *)v.data();
    for (uint64_t i = 0; i < v.size(); i++) {
      instream << p[i];
    }
    instream.seekg(0, instream.end);
    instream.seekg(0, instream.beg);
    cereal::PortableBinaryInputArchive iarchive(instream);
    iarchive(sb);
  }

  PackedBytes pack_states(const Ioss::Region &r)
  {
    // meta, entry, data, data_size
    auto state_count = r.get_property("state_count").get_int();

    auto data_size = state_count * sizeof(state_entry_t::basic_type);

    auto v = PackedBytes(sizeof(meta_entry_t) + sizeof(state_entry_t) + data_size);

    auto meta          = reinterpret_cast<meta_entry_t *>(v.data());
    meta->ioss_type    = meta_entry_t::IossType::Ios3States;
    meta->value.offset = 0;
    meta->value.size   = sizeof(state_entry_t) + data_size;

    auto entry = reinterpret_cast<state_entry_t *>(
        reinterpret_cast<void *>(v.data() + sizeof(meta_entry_t)));
    entry->count        = state_count;
    entry->value.offset = 0;
    entry->value.size   = data_size;

    auto data = reinterpret_cast<Ios3::state_entry_t::basic_type *>(
        reinterpret_cast<void *>(entry->data + entry->value.offset));

    for (uint64_t state(1); state <= entry->count; state++) {
      data[state - 1] = r.get_state_time(state);
    }

    return v;
  };

  key_t make_map_search_key(int rank, const std::string &type, const std::string &name)
  {
    std::string key = "::Map::" + type;
    if (not name.empty()) {
      key += "::Name::" + name;
    }
    return key_t(std::to_string(rank), key);
  }

  key_t make_map_key(int rank, const std::string &type, const std::string &name)
  {
    return key_t(std::to_string(rank), "::Map::" + type + "::Name::" + name);
  }

  key_t make_node_map_search_key(int rank, const std::string &name)
  {
    return make_map_search_key(rank, "NodeMap", name);
  }

  key_t make_node_map_key(int rank, const std::string &name)
  {
    return make_map_key(rank, "NodeMap", name);
  }

  key_t make_edge_map_search_key(int rank, const std::string &name)
  {
    return make_map_search_key(rank, "EdgeMap", name);
  }

  key_t make_edge_map_key(int rank, const std::string &name)
  {
    return make_map_key(rank, "EdgeMap", name);
  }

  key_t make_face_map_search_key(int rank, const std::string &name)
  {
    return make_map_search_key(rank, "FaceMap", name);
  }

  key_t make_face_map_key(int rank, const std::string &name)
  {
    return make_map_key(rank, "FaceMap", name);
  }

  key_t make_elem_map_search_key(int rank, const std::string &name)
  {
    return make_map_search_key(rank, "ElemMap", name);
  }

  key_t make_elem_map_key(int rank, const std::string &name)
  {
    return make_map_key(rank, "ElemMap", name);
  }

  key_t make_states_search_key(int rank, const Ioss::Region &region)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return key_t(std::to_string(rank), "::TimeSteps");
  }

  key_t make_states_key(int rank, const Ioss::Region &region)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return key_t(std::to_string(rank), "::TimeSteps");
  }

  key_t sideblocks_search_key(int rank, const Ioss::Region &region, const Ioss::SideSet &sideset)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    std::string value = fmt::format("::SideSet::{}::SideBlock::", sideset.name());
    return key_t(std::to_string(rank), value);
  }

  key_t make_sideblock_key(int rank, const Ioss::Region &region, const Ioss::SideSet &sideset,
                           const Ioss::SideBlock &sideblock)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    std::string value =
        fmt::format("::SideSet::{}::SideBlock::{}", sideset.name(), sideblock.name());
    return key_t(std::to_string(rank), value);
  }

  key_t structuredblock_search_key(int rank, const Ioss::Region &region,
                                   const Ioss::StructuredBlock &structuredblock)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    std::string value = fmt::format("::StructuredBlock::{}::Attributes", structuredblock.name());
    return key_t(std::to_string(rank), value);
  }

  key_t make_structuredblock_key(int rank, const Ioss::Region &region,
                                 const Ioss::StructuredBlock &structuredblock)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    std::string value = fmt::format("::StructuredBlock::{}::Attributes", structuredblock.name());
    return key_t(std::to_string(rank), value);
  }

  key_t make_key(int rank, const Ioss::Region &region, const Ioss::GroupingEntity &grouping_entity,
                 const Ioss::Field &field, const std::string &name)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value =
        fmt::format("::State::{}::Entity::{}{}::Field::RoleType::{}::BasicType::{}::Name::{}",
                    region.get_current_state(), grouping_entity.type_string(), grouping_entity_name,
                    to_string(field.get_role()), to_string(field.get_type()), name);
    return key_t(std::to_string(rank), value);
  }

  key_t make_key(int rank, const Ioss::Region &region, const Ioss::GroupingEntity &grouping_entity,
                 const Ioss::Field &field)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value =
        fmt::format("::State::{}::Entity::{}{}::Field::RoleType::{}::BasicType::{}::Name::{}",
                    region.get_current_state(), grouping_entity.type_string(), grouping_entity_name,
                    to_string(field.get_role()), to_string(field.get_type()), field.get_name());
    return key_t(std::to_string(rank), value);
  }

  key_t make_key(int rank, const Ioss::Region &region, const Ioss::GroupingEntity &grouping_entity,
                 const Ioss::Property &property)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value =
        fmt::format("::State::{}::Entity::{}{}::Property::BasicType::{}::Name::{}",
                    region.get_current_state(), grouping_entity.type_string(), grouping_entity_name,
                    to_string(property.get_type()), property.get_name());
    return key_t(std::to_string(rank), value);
  }

  key_t make_key(int rank, const Ioss::Region &region, const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value = fmt::format("::State::{}::Entity::{}{}", region.get_current_state(),
                                    grouping_entity.type_string(), grouping_entity_name);
    return key_t(std::to_string(rank), value);
  }

  key_t entity_search_key(int rank, const Ioss::Region &region, const std::string &entity)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    std::string value =
        fmt::format("::State::{}::Entity::{}::", region.get_current_state(), entity);
    return key_t(std::to_string(rank), value);
  }

  key_t entity_search_key(int rank, const Ioss::Region &region,
                          const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value = fmt::format("::State::{}::Entity::{}{}", region.get_current_state(),
                                    grouping_entity.type_string(), grouping_entity_name);
    return key_t(std::to_string(rank), value);
  }

  key_t property_search_key(int rank, const Ioss::Region &region,
                            const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value =
        fmt::format("::State::{}::Entity::{}{}::Property::", region.get_current_state(),
                    grouping_entity.type_string(), grouping_entity_name);
    return key_t(std::to_string(rank), value);
  }

  key_t make_property_key(int rank, const Ioss::Region &region, const std::string &entity_type,
                          const std::string &entity_name, const std::string &property_type,
                          const std::string &property_name)
  {
    std::string grouping_entity_name;
    if (entity_type == "Region") {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + entity_name;
    }

    std::string value = fmt::format("::State::{}::Entity::{}{}::Property::BasicType::{}::Name::{}",
                                    region.get_current_state(), entity_type, grouping_entity_name,
                                    property_type, property_name);
    return key_t(std::to_string(rank), value);
  }

  key_t field_search_key(int rank, const Ioss::Region &region,
                         const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value =
        fmt::format("::State::{}::Entity::{}{}::Field::", region.get_current_state(),
                    grouping_entity.type_string(), grouping_entity_name);
    return key_t(std::to_string(rank), value);
  }

  key_t field_search_key(int rank, int state, const Ioss::Region &,
                         const Ioss::GroupingEntity &grouping_entity)
  {
    // The default Region state is -1. When the mesh data is being
    // read from S3, the state has not yet been updated.  This version
    // of the function allows the caller to specify which state it
    // should be reading keys for.  This is required to allow
    // Ios3::DatabaseIO to identify the TRANSIENT fields in the S3 key
    // structure.

    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "::Name::" + grouping_entity.name();
    }

    std::string value = fmt::format("::State::{}::Entity::{}{}::Field::", state,
                                    grouping_entity.type_string(), grouping_entity_name);
    return key_t(std::to_string(rank), value);
  }

  std::string to_string(const Ioss::Property::BasicType &t)
  {
    switch (t) {
    case Ioss::Property::BasicType::INVALID: return std::string("INVALID");
    case Ioss::Property::BasicType::REAL: return std::string("REAL");
    case Ioss::Property::BasicType::INTEGER: return std::string("INTEGER");
    case Ioss::Property::BasicType::POINTER: return std::string("POINTER");
    case Ioss::Property::BasicType::STRING: return std::string("STRING");
    default: return std::string("INVALID");
    }
  }

  std::string to_string(const Ioss::Field::BasicType &t)
  {
    // INTEGER == INT32
    // DOUBLE == REAL
    switch (t) {
    case Ioss::Field::BasicType::INVALID: return std::string("INVALID");
    case Ioss::Field::BasicType::REAL: return std::string("REAL");
    case Ioss::Field::BasicType::INTEGER: return std::string("INTEGER");
    case Ioss::Field::BasicType::INT64: return std::string("INT64");
    case Ioss::Field::BasicType::COMPLEX: return std::string("COMPLEX,");
    case Ioss::Field::BasicType::STRING: return std::string("STRING,");
    case Ioss::Field::BasicType::CHARACTER: return std::string("CHARACTER");
    default: return std::string("INVALID");
    }
  }

  std::string to_string(const Ioss::Field::RoleType &t)
  {
    // INTEGER == INT32
    // DOUBLE == REAL
    switch (t) {
    case Ioss::Field::RoleType::INTERNAL: return std::string("INTERNAL");
    case Ioss::Field::RoleType::MESH: return std::string("MESH");
    case Ioss::Field::RoleType::ATTRIBUTE: return std::string("ATTRIBUTE");
    case Ioss::Field::RoleType::COMMUNICATION: return std::string("COMMUNICATION");
    case Ioss::Field::RoleType::INFORMATION: return std::string("INFORMATION");
    case Ioss::Field::RoleType::REDUCTION: return std::string("REDUCTION");
    case Ioss::Field::RoleType::TRANSIENT: return std::string("TRANSIENT");
    default: return std::string("INVALID");
    }
  }

  std::string to_string(const Ioss::EntityType &t) { return Ioss::Utils::entity_type_to_string(t); }

  std::string get_entity_name(const std::string &k, const std::string &target)
  {
    std::string name;

    std::string front = fmt::format("{}::", target);
    auto        begin = k.find(front);
    if (begin != std::string::npos) {
      name = k.substr(begin + front.size());
    }

    return name;
  }

  std::set<std::string> get_entity_names(const std::vector<std::string> &keys,
                                         const std::string              &target)
  {
    std::set<std::string> names;
    for (auto k : keys) {
      std::string front = fmt::format("Entity::{}::Name::", target);
      auto        begin = k.find(front);
      if (begin != std::string::npos) {
        {
          std::string back = std::string("::Property::BasicType");
          auto        end  = k.find(back);
          if (end != std::string::npos) {
            auto name = k.substr(begin + front.size(), end - begin - front.size());
            names.insert(name);
          }
        }
        {
          std::string back = std::string("::Field::RoleType");
          auto        end  = k.find(back);
          if (end != std::string::npos) {
            auto name = k.substr(begin + front.size(), end - begin - front.size());
            names.insert(name);
          }
        }
      }
    }

    return names;
  }

} // namespace Ios3
