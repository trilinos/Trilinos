// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "ios3_export.h"

#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_IOFactory.h"  // for IOFactory

#include "Ios3_AwsHelpers.h"

namespace Ios3 {

  class IOS3_EXPORT IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              Ioss_MPI_Comm                communicator,
                              const Ioss::PropertyManager &properties) const;
  };

  class IOS3_EXPORT DatabaseIO : public Ioss::DatabaseIO
  {
  private:
    bool put_properties() const;
    void finalize_database() const override;

  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               Ioss_MPI_Comm communicator, const Ioss::PropertyManager &props);

    ~DatabaseIO() override;

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    IOSS_NODISCARD virtual unsigned entity_field_support() const override
    {
      return Ioss::NODEBLOCK | Ioss::EDGEBLOCK | Ioss::FACEBLOCK | Ioss::ELEMENTBLOCK |
             Ioss::NODESET | Ioss::EDGESET | Ioss::FACESET | Ioss::ELEMENTSET | Ioss::SIDESET |
             Ioss::SIDEBLOCK | Ioss::REGION | Ioss::SUPERELEMENT;
    }

    int spatial_dimension() const { return spatialDimension; }
    int node_block_count() const { return nodeBlockCount; }
    int element_block_count() const { return elementBlockCount; }
    int nodeset_count() const { return nodesetCount; }
    int sideset_count() const { return sidesetCount; }

    /** Return a string specifying underlying format of database (exodus, cgns, ...) */
    IOSS_NODISCARD virtual std::string get_format() const override { return "s3"; };

    IOSS_NODISCARD virtual int int_byte_size_db() const override { return sizeof(int); }

    virtual int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    IOSS_NOOP_GFI(Ioss::Assembly)
    IOSS_NOOP_GFI(Ioss::Blob)

    virtual int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    virtual int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                       void *data, size_t data_size) const override;
    IOSS_NOOP_PFI(Ioss::Assembly)
    IOSS_NOOP_PFI(Ioss::Blob)

  private:
    int64_t get_field_internal(const Ioss::GroupingEntity &e, const Ioss::Field &f, void *data,
                               size_t data_size) const;
    int64_t put_field_internal(const Ioss::GroupingEntity &e, const Ioss::Field &f, void *data,
                               size_t data_size) const;

    virtual bool begin_nl(Ioss::State state) override;
    virtual bool end_nl(Ioss::State state) override;

    virtual void read_meta_data_nl() override;

    void put_qa();
    void put_info();
    void write_meta_data(Ioss::IfDatabaseExistsBehavior behavior);

    /*
     * TODO identify all the get_*{blocks|sets} needed here
     */
    void get_step_times_nl() override;

    void read_region();
    void read_entity_properties(std::vector<std::string> &keys, Ioss::GroupingEntity &entity);
    Ioss::Property read_property(std::vector<unsigned char> &value);
    void           read_entity_fields(std::vector<std::string> &keys, Ioss::GroupingEntity &entity);

    Ioss::Map &get_node_map() const;

    void get_edgeblocks();
    void get_elemblocks();
    void get_faceblocks();
    void get_nodeblocks();
    void get_structuredblocks();

    void get_edgesets();
    void get_elemsets();
    void get_facesets();
    void get_nodesets();
    void get_sidesets();
    void get_commsets();

    int spatialDimension{3};

    int nodeBlockCount{0};
    int elementBlockCount{0};
    int nodesetCount{0};
    int sidesetCount{0};

    mutable Ioss::Map nodeMap;

    Ios3::helpers::HelperParameters               helper_params;
    std::shared_ptr<Ios3::helpers::HelperContext> helper_context;

    std::string bucket_name;
  };

} // namespace Ios3
