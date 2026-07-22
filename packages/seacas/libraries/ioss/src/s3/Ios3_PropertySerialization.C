// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "Ios3_PropertySerialization.h"
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

namespace Ios3 {

  size_t data_size(const Ioss::Property &p)
  {
    auto type = p.get_type();
    if (type == Ioss::Property::BasicType::REAL)
      return sizeof(double);
    else if (type == Ioss::Property::BasicType::INTEGER)
      return sizeof(int64_t);
    else if (type == Ioss::Property::BasicType::POINTER)
      return sizeof(int64_t);
    else if (type == Ioss::Property::BasicType::STRING)
      return p.get_string().size();
    else
      return 0;
  }

  int map_properties(const Ioss::Region &region, const Ioss::GroupingEntity &entity,
                     PropertyFunction op)
  {
    int num_failed = 0;

    std::vector<std::string> description;
    entity.property_describe(&description);
    for (auto name : description)
      num_failed += op(region, entity, entity.get_property(name));

    return num_failed;
  }

  int map_properties(const Ioss::Region &region, PropertyFunction op)
  {
    int num_failed = 0;

    num_failed += map_properties(region, region, op);

    for (auto entity : region.get_edge_blocks())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_element_blocks())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_face_blocks())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_node_blocks())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_structured_blocks())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_edgesets())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_elementsets())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_facesets())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_nodesets())
      num_failed += map_properties(region, *entity, op);

    for (auto entity : region.get_commsets())
      num_failed += map_properties(region, *entity, op);

    for (auto sideset : region.get_sidesets()) {
      num_failed += map_properties(region, *sideset, op);
      for (auto sideblock : sideset->get_side_blocks()) {
        num_failed += map_properties(region, *sideblock, op);
      }
    }

    return num_failed;
  }

  property_entry_t::property_entry_t(const Ioss::Property &property, const size_t start)
      : basic_type(property.get_type()), is_implicit(property.is_implicit()),
        is_valid(property.is_valid()), name{start, property.get_name().size()},
        value{name.offset + name.size, Ios3::data_size(property)}, data_size(name.size + value.size)
  {
  }

  PackedBytes pack_property(const Ioss::Region &, const Ioss::GroupingEntity &,
                            const Ioss::Property &property)
  {
    property_entry_t property_entry(property);

    PackedBytes v(sizeof(property_entry_t) + property_entry.data_size);

    // copy property_entry_t to meta section
    std::memcpy(reinterpret_cast<char *>(Data(v)), &property_entry, sizeof(property_entry_t));

    auto entry    = reinterpret_cast<property_entry_t *>(Data(v));
    auto name_ptr = reinterpret_cast<char *>(entry->data) + entry->name.offset;
    auto value_ptr =
        reinterpret_cast<void *>(reinterpret_cast<char *>(entry->data) + entry->value.offset);

    // copy name to data section
    std::memcpy(name_ptr, property.get_name().data(), entry->name.size);

    // copy value to data section
    if (property.get_type() == Ioss::Property::BasicType::INTEGER) {
      auto value = reinterpret_cast<int64_t *>(value_ptr);
      *value     = property.get_int();
    }
    else if (property.get_type() == Ioss::Property::BasicType::REAL) {
      auto value = reinterpret_cast<double *>(value_ptr);
      *value     = property.get_int();
    }
    else if (property.get_type() == Ioss::Property::BasicType::STRING) {
      std::memcpy(reinterpret_cast<char *>(value_ptr), property.get_string().data(),
                  property.get_string().size());
    }

    if (entry->value.size != Ios3::data_size(property)) {
      fmt::print(stderr, "pack_property: value.size mismatch: {} != {}", entry->value.size,
                 Ios3::data_size(property));
    }
    return v;
  }

  int64_t property_get_int(PackedBytes &p)
  {
    auto prop(reinterpret_cast<Ios3::property_entry_t *>(Data(p)));
    auto value_ptr = reinterpret_cast<void *>(prop->data + prop->value.offset);
    return *(reinterpret_cast<int64_t *>(value_ptr));
  }

  std::string property_get_string(PackedBytes &p)
  {
    auto prop(reinterpret_cast<Ios3::property_entry_t *>(Data(p)));
    return std::string(prop->data + prop->value.offset, prop->value.size);
  }
} // namespace Ios3
