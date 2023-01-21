// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef Iovs_cgns_DatabaseIO_h
#define Iovs_cgns_DatabaseIO_h

#include "iovs_export.h"

#include <Ioss_DatabaseIO.h>

namespace Iovs_cgns {
  class CatalystCGNSMeshBase;

  class IOVS_EXPORT DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               Ioss_MPI_Comm communicator, const Ioss::PropertyManager &props);

    ~DatabaseIO() override;

    const std::string get_format() const override { return "Embedded CGNS Visualization"; }

    unsigned entity_field_support() const override { return Ioss::REGION; }

    int  int_byte_size_db() const override { return int_byte_size_api(); }
    void write_meta_data();

  private:
    bool begin__(Ioss::State state) override;
    bool end__(Ioss::State state) override;

    bool begin_state__(int state, double time) override;
    bool end_state__(int state, double time) override;

    void read_meta_data__() override;

    IOSS_NOOP_GFI(Ioss::Region)
    IOSS_NOOP_GFI(Ioss::NodeBlock)
    IOSS_NOOP_GFI(Ioss::EdgeBlock)
    IOSS_NOOP_GFI(Ioss::FaceBlock)
    IOSS_NOOP_GFI(Ioss::ElementBlock)
    IOSS_NOOP_GFI(Ioss::StructuredBlock)
    IOSS_NOOP_GFI(Ioss::SideBlock)
    IOSS_NOOP_GFI(Ioss::NodeSet)
    IOSS_NOOP_GFI(Ioss::EdgeSet)
    IOSS_NOOP_GFI(Ioss::FaceSet)
    IOSS_NOOP_GFI(Ioss::ElementSet)
    IOSS_NOOP_GFI(Ioss::SideSet)
    IOSS_NOOP_GFI(Ioss::CommSet)
    IOSS_NOOP_GFI(Ioss::Assembly)
    IOSS_NOOP_GFI(Ioss::Blob)

    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override;
    IOSS_NOOP_PFI(Ioss::Region)
    IOSS_NOOP_PFI(Ioss::NodeBlock)
    IOSS_NOOP_PFI(Ioss::EdgeBlock)
    IOSS_NOOP_PFI(Ioss::FaceBlock)
    IOSS_NOOP_PFI(Ioss::SideBlock)
    IOSS_NOOP_PFI(Ioss::NodeSet)
    IOSS_NOOP_PFI(Ioss::EdgeSet)
    IOSS_NOOP_PFI(Ioss::FaceSet)
    IOSS_NOOP_PFI(Ioss::ElementSet)
    IOSS_NOOP_PFI(Ioss::SideSet)
    IOSS_NOOP_PFI(Ioss::CommSet)
    IOSS_NOOP_PFI(Ioss::Assembly)
    IOSS_NOOP_PFI(Ioss::Blob)

    std::unique_ptr<CatalystCGNSMeshBase> catCGNSMesh;
  };
} // namespace Iovs_cgns

#endif
