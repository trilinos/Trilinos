// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "ios3_export.h"

#include "Ioss_GroupingEntity.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"

#include <set>

namespace Ios3 {

  using PackedBytes = std::vector<unsigned char>;

  using key_t = std::pair<std::string, std::string>;

  struct IOS3_EXPORT value_entry_t
  {
    size_t offset{0};
    size_t size{0};
  };

  struct IOS3_EXPORT meta_entry_t
  {
    enum class IossType { IossProperty, IossField, Ios3States, Ios3SideBlock, Ios3StructuredBlock };
    IossType      ioss_type;
    value_entry_t value;
    // NOTE an added char data[0] would point to the next meta_entry_T
  };

  struct IOS3_EXPORT state_entry_t
  {
    using basic_type = double;
    size_t        count{0};
    value_entry_t value;
    char          data[0];

    explicit state_entry_t(const Ioss::Region &r);
  };

  struct IOS3_EXPORT sideblock_entry_t
  {
    size_t entity_count{0};

    explicit sideblock_entry_t(const Ioss::SideBlock &sb);
  };

  IOS3_EXPORT PackedBytes pack_states(const Ioss::Region &r);

  IOS3_EXPORT PackedBytes pack_sideblock(const Ioss::SideBlock &sb);
  IOS3_EXPORT int64_t     unpack_sideblocks(PackedBytes &v);

  IOS3_EXPORT PackedBytes pack_structuredblock(const Ioss::StructuredBlock &sb);
  IOS3_EXPORT void        unpack_structuredblock(PackedBytes &v, Ioss::StructuredBlock &sb);

  IOS3_EXPORT key_t make_node_map_search_key(int rank, const std::string &name = "");
  IOS3_EXPORT key_t make_node_map_key(int rank, const std::string &name);
  IOS3_EXPORT key_t make_edge_map_search_key(int rank, const std::string &name = "");
  IOS3_EXPORT key_t make_edge_map_key(int rank, const std::string &name);
  IOS3_EXPORT key_t make_face_map_search_key(int rank, const std::string &name = "");
  IOS3_EXPORT key_t make_face_map_key(int rank, const std::string &name);
  IOS3_EXPORT key_t make_elem_map_search_key(int rank, const std::string &name = "");
  IOS3_EXPORT key_t make_elem_map_key(int rank, const std::string &name);

  IOS3_EXPORT key_t make_states_search_key(int parallel_rank, const Ioss::Region &region);
  IOS3_EXPORT key_t make_states_key(int parallel_rank, const Ioss::Region &region);

  IOS3_EXPORT key_t sideblocks_search_key(int rank, const Ioss::Region &region,
                                          const Ioss::SideSet &sideset);
  IOS3_EXPORT key_t make_sideblock_key(int rank, const Ioss::Region &region,
                                       const Ioss::SideSet   &sideset,
                                       const Ioss::SideBlock &sideblock);

  IOS3_EXPORT key_t structuredblock_search_key(int parallel_rank, const Ioss::Region &region,
                                               const Ioss::StructuredBlock &structuredblock);
  IOS3_EXPORT key_t make_structuredblock_key(int parallel_rank, const Ioss::Region &region,
                                             const Ioss::StructuredBlock &structuredblock);

  IOS3_EXPORT key_t make_key(int parallel_rank, const Ioss::Region &region,
                             const Ioss::GroupingEntity &grouping_entity, const Ioss::Field &field,
                             const std::string &name);
  IOS3_EXPORT key_t make_key(int parallel_rank, const Ioss::Region &region,
                             const Ioss::GroupingEntity &grouping_entity, const Ioss::Field &field);

  IOS3_EXPORT key_t make_key(int parallel_rank, const Ioss::Region &region,
                             const Ioss::GroupingEntity &grouping_entity,
                             const Ioss::Property       &property);

  IOS3_EXPORT key_t make_key(int parallel_rank, const Ioss::Region &region,
                             const Ioss::GroupingEntity &grouping_entity);

  IOS3_EXPORT key_t entity_search_key(int rank, const Ioss::Region &region,
                                      const std::string &entity_name);
  IOS3_EXPORT key_t entity_search_key(int rank, const Ioss::Region &region,
                                      const Ioss::GroupingEntity &entity);

  IOS3_EXPORT key_t property_search_key(int parallel_rank, const Ioss::Region &region,
                                        const Ioss::GroupingEntity &grouping_entity);
  IOS3_EXPORT key_t make_property_key(int rank, const Ioss::Region &region,
                                      const std::string &entity_type,
                                      const std::string &entity_name,
                                      const std::string &property_type,
                                      const std::string &property_name);

  IOS3_EXPORT key_t field_search_key(int parallel_rank, const Ioss::Region &region,
                                     const Ioss::GroupingEntity &grouping_entity);
  IOS3_EXPORT key_t field_search_key(int parallel_rank, int state, const Ioss::Region &region,
                                     const Ioss::GroupingEntity &grouping_entity);

  IOS3_EXPORT std::string to_string(const Ioss::Property::BasicType &t);
  IOS3_EXPORT std::string to_string(const Ioss::Field::BasicType &t);
  IOS3_EXPORT std::string to_string(const Ioss::Field::RoleType &t);
  IOS3_EXPORT std::string to_string(const Ioss::EntityType &t);

  IOS3_EXPORT std::string get_entity_name(const std::string &k, const std::string &target);
  IOS3_EXPORT std::set<std::string> get_entity_names(const std::vector<std::string> &keys,
                                                     const std::string              &target);

} // namespace Ios3
