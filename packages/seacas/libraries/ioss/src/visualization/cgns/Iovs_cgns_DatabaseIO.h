// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef Iovs_cgns_DatabaseIO_h
#define Iovs_cgns_DatabaseIO_h

#include "Ioss_DatabaseIO.h"
#include <visualization/cgns/CatalystCGNSMeshBase.h>

namespace Iovs_cgns {

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               Ioss_MPI_Comm communicator, const Ioss::PropertyManager &props);

    ~DatabaseIO() override;

    std::string get_format() const override { return "Embedded CGNS Visualization"; }

    unsigned entity_field_support() const override { return Ioss::REGION; }

    int  int_byte_size_db() const override { return int_byte_size_api(); }
    void write_meta_data();

  private:
    bool begin_nl(Ioss::State state) override;
    bool end_nl(Ioss::State state) override;

    bool begin_state_nl(int state, double time) override;
    bool end_state_nl(int state, double time) override;

    void read_meta_data_nl() override;

    int64_t get_field_internal(const Ioss::Region * /*reg*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceBlock * /*fb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideBlock * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeSet * /*es*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceSet * /*fs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementSet * /*es*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideSet * /*ss*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::CommSet * /*cs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::StructuredBlock * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::Assembly * /*as*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::Blob * /*bl*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }

    int64_t put_field_internal(const Ioss::Region * /*reg*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::NodeBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::EdgeBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceBlock * /*fb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::SideBlock * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::NodeSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::EdgeSet * /*es*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceSet * /*fs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementSet * /*es*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::SideSet * /*ss*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::CommSet * /*cs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override;
    int64_t put_field_internal(const Ioss::Assembly * /*as*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::Blob * /*bl*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }

    std::unique_ptr<CatalystCGNSMeshBase> catCGNSMesh;
    bool                                  isIdOutputCreated;
    void initZoneDataFromStructuredBlock(CatalystCGNSMeshBase::ZoneData &zoneData,
                                         const Ioss::StructuredBlock    *sb) const;
    void createIdOutput();
    void outputId(const std::string idName, std::vector<int64_t> &ids, bool isCellField,
                  const Ioss::StructuredBlock *sb);
  };
} // namespace Iovs_cgns

#endif
