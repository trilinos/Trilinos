// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Iofaodel_FieldSerialization.h"
#include "Iofaodel_Utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

namespace Iofaodel {

  size_t data_size(const Ioss::Field &f) { return f.get_size(); }

  void map_fields(const Ioss::Region &region, const Ioss::GroupingEntity &entity, FieldFunction op)
  {
    std::vector<std::string> description;
    entity.field_describe(&description);
    for (auto name : description)
      op(region, entity, entity.get_field(name));
  }

  void map_fields(const Ioss::Region &region, FieldFunction op)
  {

    std::vector<std::string> description;
    region.field_describe(&description);
    for (auto name : description) {
      auto role = region.get_field(name).get_role();
      if (role == Ioss::Field::RoleType::TRANSIENT || role == Ioss::Field::RoleType::REDUCTION) {
        op(region, region, region.get_field(name));
      }
    }

    for (auto entity : region.get_edge_blocks())
      map_fields(region, *entity, op);

    for (auto entity : region.get_element_blocks())
      map_fields(region, *entity, op);

    for (auto entity : region.get_face_blocks())
      map_fields(region, *entity, op);

    for (auto entity : region.get_node_blocks())
      map_fields(region, *entity, op);

    for (auto entity : region.get_edgesets())
      map_fields(region, *entity, op);

    for (auto entity : region.get_elementsets())
      map_fields(region, *entity, op);

    for (auto entity : region.get_facesets())
      map_fields(region, *entity, op);

    for (auto entity : region.get_nodesets())
      map_fields(region, *entity, op);

    for (auto entity : region.get_sidesets())
      map_fields(region, *entity, op);
  }

  // Put this in the meta data section of the LDO
  field_entry_t::field_entry_t(const Ioss::Field &field, const size_t start)
      : basic_type(field.get_type()), role_type(field.get_role()), is_valid(field.is_valid()),
        raw_count(field.raw_count()), name{start, field.get_name().size()},
        value{name.offset + name.size, Iofaodel::data_size(field)},
        storage{value.offset + value.size, field.raw_storage()->name().size()},
        data_size(name.size + value.size + storage.size)
  {
  }

  lunasa::DataObject pack_field(const Ioss::Region &region, const Ioss::GroupingEntity &entity,
                                const Ioss::Field &field)
  {
    field_entry_t field_entry(field);

    meta_entry_t meta_entry{meta_entry_t::IossType::IossField, 0, field_entry.data_size};

    auto ldo =
        lunasa::DataObject(sizeof(meta_entry_t), sizeof(field_entry_t) + field_entry.data_size,
                           lunasa::DataObject::AllocatorType::eager);

    // copy index to meta section
    std::memcpy(static_cast<char *>(ldo.GetMetaPtr()), &meta_entry, sizeof(meta_entry_t));

    // copy property_entry_t to meta section
    std::memcpy(static_cast<char *>(ldo.GetDataPtr()), &field_entry, sizeof(field_entry_t));

    auto entry       = static_cast<field_entry_t *>(ldo.GetDataPtr());
    auto name_ptr    = static_cast<char *>(entry->data) + entry->name.offset;
    auto value_ptr   = static_cast<void *>(entry->data + entry->value.offset);
    auto storage_ptr = static_cast<char *>(entry->data) + entry->storage.offset;

    // copy name to data section
    std::memcpy(name_ptr, field.get_name().data(), entry->name.size);

    std::memcpy(storage_ptr, field.raw_storage()->name().data(), entry->storage.size);

    // Copy field data from GroupingEntity to the LDO
    std::string field_name(name_ptr, entry->name.size);
    entity.get_field_data(field_name, value_ptr, entry->value.size);

    return ldo;
  }

  lunasa::DataObject pack_field(const Ioss::Region &r, const Ioss::GroupingEntity &e,
                                const Ioss::Field &f, void *data, size_t data_size)
  {
    field_entry_t field_entry(f);

    meta_entry_t meta_entry{meta_entry_t::IossType::IossField, 0, field_entry.data_size};

    auto ldo =
        lunasa::DataObject(sizeof(field_entry_t), sizeof(field_entry_t) + field_entry.data_size,
                           lunasa::DataObject::AllocatorType::eager);

    // copy index to meta section
    std::memcpy(static_cast<char *>(ldo.GetMetaPtr()), &meta_entry, sizeof(meta_entry_t));

    // copy property_entry_t to meta section
    std::memcpy(static_cast<char *>(ldo.GetDataPtr()), &field_entry, sizeof(field_entry_t));

    auto entry       = static_cast<field_entry_t *>(ldo.GetDataPtr());
    auto name_ptr    = static_cast<char *>(entry->data) + entry->name.offset;
    auto value_ptr   = static_cast<void *>(static_cast<char *>(entry->data) + entry->value.offset);
    auto storage_ptr = static_cast<char *>(entry->data) + entry->storage.offset;

    // copy name to data section
    std::memcpy(name_ptr, f.get_name().data(), entry->name.size);

    std::memcpy(storage_ptr, f.raw_storage()->name().data(), entry->storage.size);

    // Copy field data from GroupingEntity to the LDO
    std::string field_name(name_ptr, entry->name.size);
    // entity.get_field_data(field_name, value_ptr, entry->value.size);
    std::memcpy(value_ptr, data, entry->value.size);

    return ldo;
  }

} // namespace Iofaodel
