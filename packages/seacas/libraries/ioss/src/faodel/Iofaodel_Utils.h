// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <Ioss_GroupingEntity.h>
#include <Ioss_Property.h>
#include <Ioss_Region.h>

#include <kelpie/Key.hh>
#include <lunasa/DataObject.hh>

#include <set>

namespace Iofaodel {

  // Keys and LDOs go together, the pair below makes publishing easy
  using DataPair = std::pair<kelpie::Key, lunasa::DataObject>;

  struct value_entry_t
  {
    size_t offset, size;
  };

  struct meta_entry_t
  {
    enum class IossType {
      IossProperty,
      IossField,
      IofaodelStates,
      IofaodelSideBlock,
      IofaodelStructuredBlock
    };
    IossType      ioss_type;
    value_entry_t value; // offset from LDO::GetDataPtr and size
    // NOTE an added char data[0] would point to the next meta_entry_T
  };

  struct state_entry_t
  {
    using basic_type = double;
    size_t        count;
    value_entry_t value;
    char          data[0];

    explicit state_entry_t(const Ioss::Region &r);
  };

  struct sideblock_entry_t
  {
    size_t entity_count;

    explicit sideblock_entry_t(const Ioss::SideBlock &sb);
  };

  lunasa::DataObject pack_states(const Ioss::Region &r);

  lunasa::DataObject pack_sideblock(const Ioss::SideBlock &sb);
  int64_t            unpack_sideblocks(lunasa::DataObject ldo);

  lunasa::DataObject pack_structuredblock(const Ioss::StructuredBlock &sb);
  void               unpack_structuredblock(lunasa::DataObject &ldo, Ioss::StructuredBlock &sb);

  kelpie::Key make_states_search_key(int parallel_rank, const Ioss::Region &region);

  kelpie::Key make_states_key(int parallel_rank, const Ioss::Region &region);

  kelpie::Key sideblocks_search_key(int rank, const Ioss::Region &region,
                                    const Ioss::SideSet &sideset);

  kelpie::Key make_sideblock_key(int rank, const Ioss::Region &region, const Ioss::SideSet &sideset,
                                 const Ioss::SideBlock &sideblock);

  kelpie::Key structuredblock_search_key(int parallel_rank, const Ioss::Region &region,
                                         const Ioss::StructuredBlock &structuredblock);

  kelpie::Key make_structuredblock_key(int parallel_rank, const Ioss::Region &region,
                                       const Ioss::StructuredBlock &structuredblock);

  kelpie::Key make_key(int parallel_rank, const Ioss::Region &region,
                       const Ioss::GroupingEntity &grouping_entity, const Ioss::Field &field);

  kelpie::Key make_key(int parallel_rank, const Ioss::Region &region,
                       const Ioss::GroupingEntity &grouping_entity, const Ioss::Property &property);

  kelpie::Key make_key(int parallel_rank, const Ioss::Region &region,
                       const Ioss::GroupingEntity &grouping_entity);

  kelpie::Key entity_search_key(int rank, const Ioss::Region &region,
                                const std::string &entity_name);

  kelpie::Key entity_search_key(int rank, const Ioss::Region &region,
                                const Ioss::GroupingEntity &entity);

  kelpie::Key property_search_key(int parallel_rank, const Ioss::Region &region,
                                  const Ioss::GroupingEntity &grouping_entity);

  kelpie::Key make_property_key(int rank, const Ioss::Region &region,
                                const std::string &entity_type, const std::string &entity_name,
                                const std::string &property_type, const std::string &property_name);

  kelpie::Key field_search_key(int parallel_rank, const Ioss::Region &region,
                               const Ioss::GroupingEntity &grouping_entity);

  kelpie::Key field_search_key(int parallel_rank, int state, const Ioss::Region &region,
                               const Ioss::GroupingEntity &grouping_entity);

  std::string to_string(const Ioss::Property::BasicType &t);
  std::string to_string(const Ioss::Field::BasicType &t);
  std::string to_string(const Ioss::Field::RoleType &t);
  std::string to_string(const Ioss::EntityType &t);

  std::string           get_entity_name(const kelpie::Key &k, const std::string &target);
  std::set<std::string> get_entity_names(const std::vector<kelpie::Key> &keys,
                                         const std::string              &target);

} // namespace Iofaodel
