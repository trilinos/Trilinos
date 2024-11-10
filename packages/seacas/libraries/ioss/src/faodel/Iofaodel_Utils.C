// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Iofaodel_Utils.h"

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

namespace Iofaodel {

  state_entry_t::state_entry_t(const Ioss::Region &r)
      : count(r.get_property("state_count").get_int()),
        value{0, count * sizeof(state_entry_t::basic_type)}
  {
  }

  sideblock_entry_t::sideblock_entry_t(const Ioss::SideBlock &sb) : entity_count(sb.entity_count())
  {
  }

  /* The constructor for SideBlocks requires attributes that are neither properties nor fields.
     This function collects and serializes these attributes so that the SideBlock can be
     reconstructed from the data stored in Faodel. */
  lunasa::DataObject pack_sideblock(const Ioss::SideBlock &sb)
  {
    // meta, entry, data, data_size
    auto ldo = lunasa::DataObject(sizeof(meta_entry_t), sizeof(sideblock_entry_t),
                                  lunasa::DataObject::AllocatorType::eager);

    auto meta          = static_cast<meta_entry_t *>(ldo.GetMetaPtr());
    meta->ioss_type    = meta_entry_t::IossType::IofaodelSideBlock;
    meta->value.offset = 0;
    meta->value.size   = sizeof(sideblock_entry_t);

    auto entry = static_cast<sideblock_entry_t *>(
        static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));

    entry->entity_count = sb.entity_count();

    return ldo;
  }

  int64_t unpack_sideblocks(lunasa::DataObject ldo)
  {
    auto meta = static_cast<meta_entry_t *>(ldo.GetMetaPtr());
    assert(meta->ioss_type == meta_entry_t::IossType::IofaodelSideBlock);

    auto entry = static_cast<sideblock_entry_t *>(
        static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));

    return entry->entity_count;
  }

  /* The constructor for StructuredBlocks requires attributes that are neither properties nor
     fields. This function collects and serializes these attributes so that the StructuredBlock can
     be reconstructed from the data stored in Faodel. */
  lunasa::DataObject pack_structuredblock(const Ioss::StructuredBlock &sb)
  {
    std::ostringstream                  outstream;
    cereal::PortableBinaryOutputArchive oarchive(outstream);
    oarchive(sb);
    std::streamoff length = outstream.tellp();

    // meta, entry, data, data_size
    auto        ldo       = lunasa::DataObject(length);
    std::string outstring = outstream.str();
    outstring.copy(static_cast<char *>(ldo.GetDataPtr()), length);
    char *p = (char *)ldo.GetDataPtr();

    return ldo;
  }

  void unpack_structuredblock(lunasa::DataObject &ldo, Ioss::StructuredBlock &sb)
  {
    std::stringstream instream;
    char             *p = (char *)ldo.GetDataPtr();
    for (int i = 0; i < ldo.GetDataSize(); i++) {
      instream << p[i];
    }
    instream.seekg(0, instream.end);
    std::streamoff length = instream.tellg();
    instream.seekg(0, instream.beg);
    cereal::PortableBinaryInputArchive iarchive(instream);
    iarchive(sb);
  }

  lunasa::DataObject pack_states(const Ioss::Region &r)
  {
    // meta, entry, data, data_size
    auto state_count = r.get_property("state_count").get_int();

    // std::cerr << "state count: " << state_count << std::endl;
    auto data_size = state_count * sizeof(state_entry_t::basic_type);

    auto ldo = lunasa::DataObject(sizeof(meta_entry_t), sizeof(state_entry_t) + data_size,
                                  lunasa::DataObject::AllocatorType::eager);

    auto meta          = static_cast<meta_entry_t *>(ldo.GetMetaPtr());
    meta->ioss_type    = meta_entry_t::IossType::IofaodelStates;
    meta->value.offset = 0;
    meta->value.size   = sizeof(state_entry_t) + data_size;

    auto entry = static_cast<state_entry_t *>(
        static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));
    entry->count        = state_count;
    entry->value.offset = 0;
    entry->value.size   = data_size;

    auto data = static_cast<Iofaodel::state_entry_t::basic_type *>(
        static_cast<void *>(entry->data + entry->value.offset));

    for (auto state(1); state <= entry->count; state++)
      data[state - 1] = r.get_state_time(state);

    return ldo;
  };

  kelpie::Key make_states_search_key(int rank, const Ioss::Region &region)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return kelpie::Key(std::to_string(rank), "/TimeSteps*");
  }

  kelpie::Key make_states_key(int rank, const Ioss::Region &region)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return kelpie::Key(std::to_string(rank), "/TimeSteps");
  }

  kelpie::Key sideblocks_search_key(int rank, const Ioss::Region &region,
                                    const Ioss::SideSet &sideset)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    return kelpie::Key(std::to_string(rank), "/SideSet/" + sideset.name() + "/SideBlock/*");
  }

  kelpie::Key make_sideblock_key(int rank, const Ioss::Region &region, const Ioss::SideSet &sideset,
                                 const Ioss::SideBlock &sideblock)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }

    return kelpie::Key(std::to_string(rank),
                       "/SideSet/" + sideset.name() + "/SideBlock/" + sideblock.name());
  }

  kelpie::Key structuredblock_search_key(int rank, const Ioss::Region &region,
                                         const Ioss::StructuredBlock &structuredblock)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return kelpie::Key(std::to_string(rank),
                       "/StructuredBlock/" + structuredblock.name() + "/Attributes*");
  }

  kelpie::Key make_structuredblock_key(int rank, const Ioss::Region &region,
                                       const Ioss::StructuredBlock &structuredblock)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return kelpie::Key(std::to_string(rank),
                       "/StructuredBlock/" + structuredblock.name() + "/Attributes");
  }

  kelpie::Key make_key(int rank, const Ioss::Region &region,
                       const Ioss::GroupingEntity &grouping_entity, const Ioss::Field &field)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           grouping_entity.type_string() + grouping_entity_name +
                           "/Field/RoleType/" + to_string(field.get_role()) + "/BasicType/" +
                           to_string(field.get_type()) + "/Name/" + field.get_name());
  }

  kelpie::Key make_key(int rank, const Ioss::Region &region,
                       const Ioss::GroupingEntity &grouping_entity, const Ioss::Property &property)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           grouping_entity.type_string() + grouping_entity_name +
                           "/Property/BasicType/" + to_string(property.get_type()) + "/Name/" +
                           property.get_name());
  }

  kelpie::Key make_key(int rank, const Ioss::Region &region,
                       const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           grouping_entity.type_string() + grouping_entity_name);
  }

  kelpie::Key entity_search_key(int rank, const Ioss::Region &region, const std::string &entity)
  {
    auto region_name = region.name();
    if (region_name.empty()) {
      region_name = "UNNAMED";
    }
    return kelpie::Key(std::to_string(rank), "/State/" +
                                                 std::to_string(region.get_current_state()) +
                                                 "/Entity/" + entity + "/*");
  }

  kelpie::Key entity_search_key(int rank, const Ioss::Region &region,
                                const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           grouping_entity.type_string() + grouping_entity_name + "*");
  }

  kelpie::Key property_search_key(int rank, const Ioss::Region &region,
                                  const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           grouping_entity.type_string() + grouping_entity_name + "/Property/*");
  }

  kelpie::Key make_property_key(int rank, const Ioss::Region &region,
                                const std::string &entity_type, const std::string &entity_name,
                                const std::string &property_type, const std::string &property_name)
  {
    std::string grouping_entity_name;
    if (entity_type == "Region") {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + entity_name;
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           entity_type + grouping_entity_name + "/Property/BasicType/" +
                           property_type + "/Name/" + property_name);
  }

  kelpie::Key field_search_key(int rank, const Ioss::Region &region,
                               const Ioss::GroupingEntity &grouping_entity)
  {
    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank),
                       "/State/" + std::to_string(region.get_current_state()) + "/Entity/" +
                           grouping_entity.type_string() + grouping_entity_name + "/Field/*");
  }

  kelpie::Key field_search_key(int rank, int state, const Ioss::Region &region,
                               const Ioss::GroupingEntity &grouping_entity)
  {
    // The default Region state is -1. When the mesh data is being read from Faodel, the state
    // has not yet been updated.  This version of the function allows the caller to specify which
    // state it should be reading keys for.  This is required to allow Iofaodel::DatabaseIO to
    // identify the TRANSIENT fields in the Faodel key structure.

    std::string grouping_entity_name;
    if (grouping_entity.type() == Ioss::EntityType::REGION) {
      grouping_entity_name = "";
    }
    else {
      grouping_entity_name = "/Name/" + grouping_entity.name();
    }

    return kelpie::Key(std::to_string(rank), "/State/" + std::to_string(state) + "/Entity/" +
                                                 grouping_entity.type_string() +
                                                 grouping_entity_name + "/Field/*");
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

  std::string get_entity_name(const kelpie::Key &k, const std::string &target)
  {
    std::string name;

    std::string front(target + "/");
    auto        begin = k.K2().find(front);
    if (begin != std::string::npos) {
      name = k.K2().substr(begin + front.size());
    }

    return name;
  }

  std::set<std::string> get_entity_names(const std::vector<kelpie::Key> &keys,
                                         const std::string              &target)
  {
    std::set<std::string> names;
    for (auto k : keys) {
      std::string front("Entity/" + target + "/Name/");
      auto        begin = k.K2().find(front);
      if (begin != std::string::npos) {

        {
          std::string back = std::string("/Property/BasicType");
          auto        end  = k.K2().find(back);
          if (end != std::string::npos) {
            auto name = k.K2().substr(begin + front.size(), end - begin - front.size());
            names.insert(name);
          }
        }

        {
          std::string back = std::string("/Field/RoleType");
          auto        end  = k.K2().find(back);
          if (end != std::string::npos) {
            auto name = k.K2().substr(begin + front.size(), end - begin - front.size());
            names.insert(name);
          }
        }
      }
    }

    return names;
  }

} // namespace Iofaodel
