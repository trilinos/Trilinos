/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*--------------------------------------------------------------------*/
/*    Copyright 2000-2010 NTESS.                         */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// -*- Mode: c++ -*-
#ifndef IOSS_Iovs_DatabaseIO_h
#define IOSS_Iovs_DatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_EntityType.h>
#include <Ioss_Field.h>
#include <Ioss_Map.h>
#include <Ioss_Utils.h>

#include <algorithm>
#include <ctime>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

class ParaViewCatalystIossAdapterBase;

/** \brief A namespace for the visualization database format.
 */
namespace Iovs {

  typedef std::set<std::pair<int64_t, int64_t>> EntityIdSet;

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               MPI_Comm communicator, const Ioss::PropertyManager &props);
    ~DatabaseIO() override;

    const std::string get_format() const override { return "Embedded Visualization"; }

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const override
    {
      return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::NODESET | Ioss::SIDESET | Ioss::SIDEBLOCK;
    }

    static int parseCatalystFile(const std::string &filepath, std::string &json_result);

    int int_byte_size_db() const override { return int_byte_size_api(); }

  private:
    bool begin__(Ioss::State state) override;
    bool end__(Ioss::State state) override;

    bool begin_state__(int state, double time) override;
    bool end_state__(int state, double time) override;

    void read_meta_data__() override;

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
    int64_t get_field_internal(const Ioss::EdgeBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementBlock * /*eb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideBlock * /*fb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideSet * /*fs*/, const Ioss::Field & /*field*/,
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

    int64_t get_field_internal(const Ioss::Assembly * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }

    int64_t get_field_internal(const Ioss::Blob * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceBlock * /*nb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementSet * /*ns*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::CommSet * /*cs*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::StructuredBlock * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /* data_size */) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::Assembly * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }

    int64_t put_field_internal(const Ioss::Blob * /*sb*/, const Ioss::Field & /*field*/,
                               void * /*data*/, size_t /*data_size*/) const override
    {
      return 0;
    }

    void write_meta_data();

    static ParaViewCatalystIossAdapterBase *
    load_plugin_library(const std::string &plugin_name, const std::string &plugin_library_name);

    static std::string create_output_file_path(const std::string &          input_deck_name,
                                               const Ioss::PropertyManager &properties);
    // static bool plugin_library_exists(const std::string &plugin_name);

    int64_t handle_node_ids(void *ids, int64_t num_to_get);
    int64_t handle_element_ids(const Ioss::ElementBlock *eb, void *ids, size_t num_to_get);

    const Ioss::Map &get_node_map() const;
    const Ioss::Map &get_element_map() const;

    DatabaseIO();                              // Do not implement
    DatabaseIO(const DatabaseIO &);            // Do not implement
    DatabaseIO &operator=(const DatabaseIO &); // Do not implement

    std::string        databaseTitle{};
    static std::string paraview_script_filename;
    std::string        catalyst_block_file_name{};
    std::string        paraview_json_parse{};
    std::string        sierra_input_deck_name{};
    std::string        catalyst_output_directory{};
    std::string        paraview_script_extra_filename{};
    int                enableLogging;
    int                debugLevel;
    int                underscoreVectors;
    int                applyDisplacements;
    int                createSideSets;
    int                createNodeSets;
    static int         useCount;
    static int         uniqueID;

    int nodeBlockCount;
    int elementBlockCount;

    // Handle to the ParaView Catalyst dynamic library
    // that is loaded via Ioss user plugin at runtime.
    ParaViewCatalystIossAdapterBase *pvcsa;
    mutable bool                     globalNodeAndElementIDsCreated;
    void                             create_global_node_and_element_ids() const;
    mutable EntityIdSet              ids_{};
  };
} // namespace Iovs

#endif
