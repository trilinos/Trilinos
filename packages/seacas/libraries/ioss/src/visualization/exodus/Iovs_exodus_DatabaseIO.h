// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// -*- Mode: c++ -*-
#ifndef IOSS_Iovs_exodus_DatabaseIO_h
#define IOSS_Iovs_exodus_DatabaseIO_h

#include "iovs_export.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EntityType.h"
#include "Ioss_Field.h"
#include "Ioss_Map.h"
#include "Ioss_Utils.h"

#include <algorithm>
#include <ctime>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

/** \brief A namespace for the visualization database format.
 */
namespace Iovs_exodus {
  class CatalystExodusMeshBase;

  typedef std::set<std::pair<int64_t, int64_t>> EntityIdSet;

  class IOVS_EXPORT DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               Ioss_MPI_Comm communicator, const Ioss::PropertyManager &props);

    ~DatabaseIO() override;

    std::string get_format() const override { return "Embedded Visualization"; }

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const override { return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK; }

    // static int parseCatalystFile(const std::string &filepath, std::string &json_result);

    int int_byte_size_db() const override { return int_byte_size_api(); }

  private:
    bool begin_nl(Ioss::State state) override;
    bool end_nl(Ioss::State state) override;

    bool begin_state_nl(int state, double time) override;
    bool end_state_nl(int state, double time) override;

    void read_meta_data_nl() override;

    // For the time being, treat vis as write only. Consider glue pipelines.
    int64_t get_field_internal(const Ioss::Region * /*reg*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceBlock * /*fb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideBlock * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeSet * /*es*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceSet * /*fs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementSet * /*es*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideSet * /*ss*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::CommSet * /*cs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
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

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override;
    int64_t put_field_internal(const Ioss::Assembly *as, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::Blob *bl, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    void write_meta_data();

    int64_t handle_node_ids(void *ids, int64_t num_to_get);
    int64_t handle_element_ids(const Ioss::ElementBlock *eb, void *ids, size_t num_to_get);

    const Ioss::Map &get_node_map() const;
    const Ioss::Map &get_element_map() const;

    DatabaseIO();                              // Do not implement
    DatabaseIO(const DatabaseIO &);            // Do not implement
    DatabaseIO &operator=(const DatabaseIO &); // Do not implement

    bool isInput;
    bool singleProcOnly; // True if history or heartbeat which is only written from proc 0...
    bool doLogging;      // True if logging field input/output
    int  nodeBlockCount;
    int  elementBlockCount;

    std::unique_ptr<CatalystExodusMeshBase> catExoMesh;

    mutable bool        globalNodeAndElementIDsCreated;
    void                create_global_node_and_element_ids() const;
    mutable EntityIdSet ids_{};
  };
} // namespace Iovs_exodus

#endif
