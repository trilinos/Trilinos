// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "ioss_export.h"

#include "Ioss_EntityType.h" // for EntityType, etc
#include <Ioss_EntityBlock.h>
#include <Ioss_Property.h> // for Property
#include <cstddef>         // for size_t
#include <cstdint>         // for int64_t
#include <string>          // for string
namespace Ioss {
  class DatabaseIO;
} // namespace Ioss
namespace Ioss {
  class Field;
} // namespace Ioss

namespace Ioss {

  /** \brief A collection of element faces with the same topology.
   */
  class IOSS_EXPORT FaceBlock : public EntityBlock
  {
  public:
    FaceBlock(DatabaseIO *io_database, const std::string &my_name, const std::string &face_type,
              int64_t number_faces);

    FaceBlock(const FaceBlock &) = default;
    ~FaceBlock() override;

    std::string type_string() const override { return "FaceBlock"; }
    std::string short_type_string() const override { return "faceblock"; }
    std::string contains_string() const override { return "Face"; }
    EntityType  type() const override { return FACEBLOCK; }

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'face_block_count' for a region.
    Property get_implicit_property(const std::string &my_name) const override;

  protected:
    int64_t internal_get_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

    int64_t internal_put_field_data(const Field &field, void *data,
                                    size_t data_size) const override;
  };
} // namespace Ioss
