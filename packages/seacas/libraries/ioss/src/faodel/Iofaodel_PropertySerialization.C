// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Iofaodel_PropertySerialization.h"
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

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

namespace Iofaodel {

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

  void map_properties(const Ioss::Region &region, const Ioss::GroupingEntity &entity,
                      PropertyFunction op)
  {
    std::vector<std::string> description;
    entity.property_describe(&description);
    for (auto name : description)
      op(region, entity, entity.get_property(name));
  }

  void map_properties(const Ioss::Region &region, PropertyFunction op)
  {
    map_properties(region, region, op);

    for (auto entity : region.get_edge_blocks())
      map_properties(region, *entity, op);

    for (auto entity : region.get_element_blocks())
      map_properties(region, *entity, op);

    for (auto entity : region.get_face_blocks())
      map_properties(region, *entity, op);

    for (auto entity : region.get_node_blocks())
      map_properties(region, *entity, op);

    for (auto entity : region.get_structured_blocks())
      map_properties(region, *entity, op);

    for (auto entity : region.get_edgesets())
      map_properties(region, *entity, op);

    for (auto entity : region.get_elementsets())
      map_properties(region, *entity, op);

    for (auto entity : region.get_facesets())
      map_properties(region, *entity, op);

    for (auto entity : region.get_nodesets())
      map_properties(region, *entity, op);

    for (auto entity : region.get_commsets())
      map_properties(region, *entity, op);

    for (auto sideset : region.get_sidesets()) {
      map_properties(region, *sideset, op);
      for (auto sideblock : sideset->get_side_blocks()) {
        map_properties(region, *sideblock, op);
      }
    }
  }

  // Put this in the meta data section of the LDO
  property_entry_t::property_entry_t(const Ioss::Property &property, const size_t start)
      : basic_type(property.get_type()),

        is_implicit(property.is_implicit()), is_valid(property.is_valid()),

        name{start, property.get_name().size()},
        value{name.offset + name.size, Iofaodel::data_size(property)},
        data_size(name.size + value.size)
  {
  }

  lunasa::DataObject pack_property(const Ioss::Region &region, const Ioss::GroupingEntity &entity,
                                   const Ioss::Property &property)
  {
    property_entry_t property_entry(property);

    meta_entry_t meta_entry{meta_entry_t::IossType::IossProperty, 0, property_entry.data_size};

    auto ldo = lunasa::DataObject(sizeof(meta_entry_t),
                                  sizeof(property_entry_t) + property_entry.data_size,
                                  lunasa::DataObject::AllocatorType::eager);

    // copy meta_entry_t to meta section
    std::memcpy(static_cast<char *>(ldo.GetMetaPtr()), &meta_entry, sizeof(meta_entry_t));

    // copy property_entry_t to meta section
    std::memcpy(static_cast<char *>(ldo.GetDataPtr()), &property_entry, sizeof(property_entry_t));

    auto entry     = static_cast<property_entry_t *>(ldo.GetDataPtr());
    auto name_ptr  = static_cast<char *>(entry->data) + entry->name.offset;
    auto value_ptr = static_cast<void *>(static_cast<char *>(entry->data) + entry->value.offset);

    // copy name to data section
    std::memcpy(name_ptr, property.get_name().data(), entry->name.size);

    // copy value to data section
    if (property.get_type() == Ioss::Property::BasicType::INTEGER) {
      auto value = static_cast<int64_t *>(value_ptr);
      *value     = property.get_int();
    }
    else if (property.get_type() == Ioss::Property::BasicType::REAL) {
      auto value = static_cast<double *>(value_ptr);
      *value     = property.get_int();
    }
    else if (property.get_type() == Ioss::Property::BasicType::STRING) {
      std::memcpy(static_cast<char *>(value_ptr), property.get_string().data(),
                  property.get_string().size());
    }

    if (entry->value.size != Iofaodel::data_size(property))
      std::cerr << "value.size mismatch: " << entry->value.size
                << " ?= " << Iofaodel::data_size(property);
    return ldo;
  }

  int64_t property_get_int(lunasa::DataObject ldo)
  {
    auto prop(static_cast<Iofaodel::property_entry_t *>(ldo.GetDataPtr()));
    auto value_ptr = static_cast<void *>(prop->data + prop->value.offset);
    return *(static_cast<int64_t *>(value_ptr));
  }

  std::string property_get_string(lunasa::DataObject ldo)
  {
    auto prop(static_cast<Iofaodel::property_entry_t *>(ldo.GetDataPtr()));
    return std::string(prop->data + prop->value.offset, prop->value.size);
  }
} // namespace Iofaodel
