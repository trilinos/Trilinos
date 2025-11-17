// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ios3_FieldSerialization.h"
#include "Ios3_Utils.h"

namespace Ios3 {

  size_t data_size(const Ioss::Field &f) { return f.get_size(); }

  int map_fields(const Ioss::Region &region, const Ioss::GroupingEntity &entity, FieldFunction op)
  {
    int num_failed = 0;

    std::vector<std::string> description;
    entity.field_describe(&description);
    for (auto name : description)
      num_failed += op(region, entity, entity.get_field(name));

    return num_failed;
  }

  int map_fields(const Ioss::Region &region, FieldFunction op)
  {
    int num_failed = 0;

    std::vector<std::string> description;
    region.field_describe(&description);
    for (auto name : description) {
      auto role = region.get_field(name).get_role();
      if (role == Ioss::Field::RoleType::TRANSIENT || role == Ioss::Field::RoleType::REDUCTION) {
        op(region, region, region.get_field(name));
      }
    }

    for (auto entity : region.get_edge_blocks())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_element_blocks())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_face_blocks())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_node_blocks())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_edgesets())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_elementsets())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_facesets())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_nodesets())
      num_failed += map_fields(region, *entity, op);

    for (auto entity : region.get_sidesets())
      num_failed += map_fields(region, *entity, op);

    return num_failed;
  }

  field_entry_t::field_entry_t(const Ioss::Field &field, const size_t start)
      : basic_type(field.get_type()), role_type(field.get_role()), is_valid(field.is_valid()),
        raw_count(field.raw_count()), name{start, field.get_name().size()},
        value{name.offset + name.size, Ios3::data_size(field)},
        storage{value.offset + value.size, field.raw_storage()->name().size()},
        data_size(name.size + value.size + storage.size)
  {
  }

  PackedBytes pack_field(const Ioss::Region &, const Ioss::GroupingEntity &entity,
                         const Ioss::Field &field)
  {
    field_entry_t field_entry(field);

    PackedBytes v(sizeof(field_entry_t) + field_entry.data_size);

    // copy field_entry_t to meta section
    std::memcpy(Data(v), &field_entry, sizeof(field_entry_t));

    auto entry       = reinterpret_cast<field_entry_t *>(Data(v));
    auto name_ptr    = reinterpret_cast<char *>(entry->data) + entry->name.offset;
    auto value_ptr   = reinterpret_cast<void *>(entry->data + entry->value.offset);
    auto storage_ptr = reinterpret_cast<char *>(entry->data) + entry->storage.offset;

    // copy name to data section
    std::memcpy(name_ptr, field.get_name().data(), entry->name.size);

    std::memcpy(storage_ptr, field.raw_storage()->name().data(), entry->storage.size);

    std::string field_name(name_ptr, entry->name.size);
    entity.get_field_data(field_name, value_ptr, entry->value.size);

    return v;
  }

  PackedBytes pack_field(const Ioss::Region &, const Ioss::GroupingEntity &, const Ioss::Field &f,
                         void *data, size_t)
  {
    field_entry_t field_entry(f);

    PackedBytes v(sizeof(field_entry_t) + field_entry.data_size);

    // copy field_entry_t to meta section
    std::memcpy(Data(v), &field_entry, sizeof(field_entry_t));

    auto entry    = reinterpret_cast<field_entry_t *>(Data(v));
    auto name_ptr = reinterpret_cast<char *>(entry->data) + entry->name.offset;
    auto value_ptr =
        reinterpret_cast<void *>(reinterpret_cast<char *>(entry->data) + entry->value.offset);
    auto storage_ptr = reinterpret_cast<char *>(entry->data) + entry->storage.offset;

    // copy name to data section
    std::memcpy(name_ptr, f.get_name().data(), entry->name.size);

    std::memcpy(storage_ptr, f.raw_storage()->name().data(), entry->storage.size);

    std::string field_name(name_ptr, entry->name.size);
    std::memcpy(value_ptr, data, entry->value.size);

    return v;
  }

} // namespace Ios3
