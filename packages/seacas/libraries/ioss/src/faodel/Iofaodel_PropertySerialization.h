// Copyright(C) 1999-2022, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iofaodel_export.h"

#include "Iofaodel_Utils.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Property.h"

#include <string>
#include <vector>

#include <lunasa/DataObject.hh>

namespace Iofaodel {

  IOFAODEL_EXPORT size_t data_size(const Ioss::Property &p);

  // Caller should write their own version of this
  // PropertyFunction should return a function or a lambda that matches the
  // signature below. The function it returns may or may not capture variables
  // that are given to the user-defined function.
  // Some examples are given in this file and are also useful
  //
  using PropertyFunction = std::function<void(const Ioss::Region &, const Ioss::GroupingEntity &,
                                              const Ioss::Property &)>;

  // Applies PropertyFunction 'op' to all properties encountered in the
  // Ioss::Region and it's various Ioss::GroupingEntities
  IOFAODEL_EXPORT void map_properties(const Ioss::Region &region, PropertyFunction op);

  // Applies PropertyFunction 'op' to all properties encountered in the
  // Ioss::GroupingEntity
  IOFAODEL_EXPORT void map_properties(const Ioss::Region         &region,
                                      const Ioss::GroupingEntity &grouping_entity,
                                      PropertyFunction            op);

  IOFAODEL_EXPORT lunasa::DataObject pack_property(const Ioss::Region         &region,
                                                   const Ioss::GroupingEntity &entity,
                                                   const Ioss::Property       &property);

  // Put this in the meta data section of the LDO
  struct IOFAODEL_EXPORT property_entry_t
  {
    Ioss::Property::BasicType basic_type;

    bool is_implicit;
    bool is_valid;

    value_entry_t name;  // offsets from data[0]
    value_entry_t value; // offsets from data[0]

    size_t data_size; // Total size of data stored in data[0] ptr
    char   data[0];

    explicit property_entry_t(const Ioss::Property &property, const size_t start = 0);
  };

  IOFAODEL_EXPORT int64_t property_get_int(lunasa::DataObject ldo);
  std::string             property_get_string(lunasa::DataObject ldo);

} // namespace Iofaodel
