// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "ios3_export.h"

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

#include <string>
#include <vector>

namespace Ios3 {

  IOS3_EXPORT size_t data_size(const Ioss::Field &f);

  // Caller should write their own version of this FieldFunction
  // should return a function or a lamba that matches the signature
  // below. The function it returns may or may not capture variables
  // that are given to the user-defined function.  Some examples are
  // given in this file and are also useful
  //
  using FieldFunction =
      std::function<int(const Ioss::Region &, const Ioss::GroupingEntity &, const Ioss::Field &)>;

  using PackedBytes = std::vector<unsigned char>;

  // Applies FieldFunction 'op' to all fields encountered in the
  // Ioss::Region and it's various Ioss::GroupingEntities
  //
  IOS3_EXPORT int map_fields(const Ioss::Region &region, FieldFunction op);

  // Applies FieldFunction 'op' to all fields encountered in the
  // Ioss::GroupingEntity
  //
  IOS3_EXPORT int map_fields(const Ioss::Region         &region,
                             const Ioss::GroupingEntity &grouping_entity, FieldFunction op);

  IOS3_EXPORT PackedBytes pack_field(const Ioss::Region &region, const Ioss::GroupingEntity &entity,
                                     const Ioss::Field &field);

  IOS3_EXPORT PackedBytes pack_field(const Ioss::Region &r, const Ioss::GroupingEntity &e,
                                     const Ioss::Field &f, void *data, size_t data_size);

  struct IOS3_EXPORT field_entry_t
  {
    Ioss::Field::BasicType basic_type;
    Ioss::Field::RoleType  role_type;
    bool                   is_implicit{false};
    bool                   is_valid{false};
    size_t                 raw_count{0};

    value_entry_t name;    // offset from data[0]
    value_entry_t value;   // offset from data[0]
    value_entry_t storage; // offset from data[0]

    size_t data_size{0}; // Total size of data stored in the value vector

    char data[0];

    explicit field_entry_t(const Ioss::Field &field, const size_t start = 0);
  };

} // namespace Ios3
