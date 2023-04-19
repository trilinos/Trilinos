// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_CodeTypes.h>
#include <null/Ionull_DatabaseIO.h>

#include "Ioss_Assembly.h"
#include "Ioss_Blob.h"
#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_EntitySet.h"
#include "Ioss_EntityType.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_State.h"
#include "Ioss_VariableType.h"

namespace Ionull {
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, Ioss_MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props)
  {
  }

  DatabaseIO::~DatabaseIO() {}

  void DatabaseIO::read_meta_data__() {}

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::EDGEBLOCK | Ioss::FACEBLOCK | Ioss::ELEMENTBLOCK |
           Ioss::NODESET | Ioss::EDGESET | Ioss::FACESET | Ioss::ELEMENTSET | Ioss::SIDESET |
           Ioss::SIDEBLOCK | Ioss::REGION | Ioss::SUPERELEMENT | Ioss::ASSEMBLY | Ioss::BLOB |
           Ioss::STRUCTUREDBLOCK;
  }

  bool DatabaseIO::begin__(Ioss::State /* state */) { return true; }

  bool DatabaseIO::end__(Ioss::State /* state */) { return true; }

  bool DatabaseIO::begin_state__(int /* state */, double) { return true; }

  bool DatabaseIO::end_state__(int /* state */, double) { return true; }

  int64_t DatabaseIO::put_field_internal(const Ioss::Region *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Assembly *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *, const Ioss::Field &field,
                                         void *, size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *, const Ioss::Field &field,
                                         void *, size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Blob *, const Ioss::Field &field, void *,
                                         size_t data_size) const
  {
    return field.verify(data_size);
  }

} // namespace Ionull
