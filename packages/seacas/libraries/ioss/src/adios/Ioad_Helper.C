// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_SideBlock.h" // for SideBlockContainer
#include "adios/Ioad_Constants.h"
#include "adios/Ioad_Helper.h"

namespace Ioad {
  int find_field_in_mapset(const std::string &entity_type, const std::string &field_name,
                           const std::map<std::string, std::set<std::string>> &mapset)
  {
    if (mapset.find(entity_type) == mapset.end()) {
      // No field for this entity_type in the map.
      return 0;
    }
    const std::set<std::string> &entity_set = mapset.at(entity_type);
    if (entity_set.find(field_name) != entity_set.end()) {
      return 1;
    }
    return 0;
  }

  std::string get_property_variable_name(const std::string &property_name)
  {
    return property_meta + property_name;
  }

  std::vector<std::string> properties_to_save(const Ioss::GroupingEntity *const entity_block)
  {
    std::vector<std::string> property_list;
    entity_block->property_describe(&property_list);

    for (auto ignore_property : Ignore_properties) {
      property_list.erase(std::remove(property_list.begin(), property_list.end(), ignore_property),
                          property_list.end());
    }
    property_list.erase(std::remove_if(property_list.begin(), property_list.end(),
                                       [&](std::string property_name) -> bool {
                                         Ioss::Property property =
                                             entity_block->get_property(property_name);
                                         return (property.is_invalid() || property.is_implicit());
                                       }),
                        property_list.end());
    return property_list;
  }

  std::string stringify_side_block_names(const Ioss::SideBlockContainer &sblocks)
  {
    std::string stringified_sblock_names;
    for (auto sblock : sblocks) {
      stringified_sblock_names += Sideblock_separator + sblock->name();
    }
    return stringified_sblock_names;
  }

  std::string encode_field_name(std::vector<std::string> names)
  {
    std::string encoded_name;
    size_t      count = 0;
    for (std::vector<std::string>::iterator it = names.begin(); it != names.end() - 1;
         ++it, ++count) {
      encoded_name += *it + Name_separator;
    }
    // sentinel value to avoid having to test if the string is empty or not.
    names.push_back("");
    return encoded_name + names[count];
  }

  std::string encode_sideblock_name(const std::string &type_string, const std::string &name)
  {
    return encode_field_name({type_string, name, sideblock_names});
  }

  bool is_sideblock_name(const std::string &name)
  {
    size_t pos = name.find(sideblock_names);
    return (pos == 0);
  }

  bool use_transformed_storage(const Ioss::Field &field, const std::string &entity_type,
                               const std::string &field_name)
  {
    return (field.get_role() == Ioss::Field::RoleType::TRANSIENT ||
            find_field_in_mapset(entity_type, field_name, Use_transformed_storage_map));
  }
} // namespace Ioad
