// Copyright(C) 1999-2022, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iofaodel_export.h"

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

#include <string>
#include <vector>

namespace Iofaodel {

  IOFAODEL_EXPORT size_t data_size(const Ioss::Field &f);

  // Caller should write their own version of this
  // FieldFunction should return a function or a lambda that matches the
  // signature below. The function it returns may or may not capture variables
  // that are given to the user-defined function.
  // Some examples are given in this file and are also useful
  //
  using FieldFunction =
      std::function<void(const Ioss::Region &, const Ioss::GroupingEntity &, const Ioss::Field &)>;

  // Applies FieldFunction 'op' to all fields encountered in the
  // Ioss::Region and it's various Ioss::GroupingEntities
  IOFAODEL_EXPORT void map_fields(const Ioss::Region &region, FieldFunction op);

  // Applies FieldFunction 'op' to all fields encountered in the
  // Ioss::GroupingEntity
  IOFAODEL_EXPORT void map_fields(const Ioss::Region         &region,
                                  const Ioss::GroupingEntity &grouping_entity, FieldFunction op);

  IOFAODEL_EXPORT lunasa::DataObject pack_field(const Ioss::Region         &region,
                                                const Ioss::GroupingEntity &entity,
                                                const Ioss::Field          &field);

  IOFAODEL_EXPORT lunasa::DataObject pack_field(const Ioss::Region         &r,
                                                const Ioss::GroupingEntity &e, const Ioss::Field &f,
                                                void *data, size_t data_size);

  // Put this in the meta data section of the LDO
  struct IOFAODEL_EXPORT field_entry_t
  {
    Ioss::Field::BasicType basic_type;
    Ioss::Field::RoleType  role_type;
    bool                   is_implicit;
    bool                   is_valid;
    size_t                 raw_count;

    // value_entry_t storage;
    value_entry_t name;
    value_entry_t value;
    value_entry_t storage;
    size_t        data_size; // Total size of data stored in LDO data section

    char data[0];

    explicit field_entry_t(const Ioss::Field &field, const size_t start = 0);
  };

} // namespace Iofaodel
