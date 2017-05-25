/*--------------------------------------------------------------------*/
/*    Copyright 2000-2010 Sandia Corporation.                         */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// -*- Mode: c++ -*-
#ifndef SIERRA_Iovs_DatabaseIO_h
#define SIERRA_Iovs_DatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_EntityType.h>
#include <Ioss_Field.h>
#include <Ioss_Map.h>
#include <Ioss_Utils.h>

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
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
               MPI_Comm communicator, const Ioss::PropertyManager &properties);
    ~DatabaseIO();

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const
    {
      return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::NODESET | Ioss::SIDESET | Ioss::SIDEBLOCK;
    }

    // Eliminate as much memory as possible, but still retain meta data information
    // Typically, eliminate the maps...
    void release_memory();

    /*!
     * Determine the local position of the node with the global id
     * 'global'.  If 'must_exist' is false, then the global id possibly
     * does not exist in the map; otherwise, it must exist and will
     * throw an exception if not found.
     */
    int64_t node_global_to_local(int64_t global, bool must_exist) const
    {
      return nodeMap.global_to_local(global, must_exist);
    }

    int64_t element_global_to_local(int64_t global) const
    {
      return elemMap.global_to_local(global);
    }

    bool begin(Ioss::State state);
    bool end(Ioss::State state);

    bool begin_state(Ioss::Region *region, int state, double time);
    bool end_state(Ioss::Region *region, int state, double time);

    void read_meta_data();

    static int parseCatalystFile(const std::string &filepath, std::string &json_result);

  private:
    // For the time being, treat vis as write only. Consider glue pipelines.
    int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }

    virtual int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                       void *data, size_t data_size) const
    {
      return 0;
    }

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const;

    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const;
    int64_t put_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const;
    int64_t put_field_internal(const Ioss::SideBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const;

    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const;
    int64_t put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const
    {
      return 0;
    }
    virtual int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                       void *data, size_t data_size) const
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

    bool isInput;
    bool singleProcOnly; // True if history or heartbeat which is only written from proc 0...
    bool doLogging;      // True if logging field input/output

    std::string        databaseTitle;
    static std::string paraview_script_filename;
    std::string        catalyst_block_file_name;
    std::string        paraview_json_parse;
    std::string        sierra_input_deck_name;
    std::string        catalyst_output_directory;
    std::string        paraview_script_extra_filename;
    int                enableLogging;
    int                debugLevel;
    int                underscoreVectors;
    int                applyDisplacements;
    int                createSideSets;
    int                createNodeSets;
    static int         useCount;
    static int         uniqueID;

    int64_t nodeCount;
    int64_t elementCount;

    int nodeBlockCount;
    int elementBlockCount;

    // Handle to the ParaView Catalyst dynamic library
    // that is loaded via Ioss user plugin at runtime.
    ParaViewCatalystIossAdapterBase *pvcsa;
    mutable bool                     globalNodeAndElementIDsCreated;
    void                             create_global_node_and_element_ids() const;
    mutable EntityIdSet              ids_;

    // Bulk Data

    // MAPS -- Used to convert from local exodusII ids/names to Ioss
    // database global ids/names
    // Maps internal (1..num_entity) ids to global ids used on the
    //               sierra side.   global = XXXMap.map[local]
    // XXXMap.map[0] contains: -1 if sequential, 0 if ordering unknown, 1
    // if nonsequential
    mutable Ioss::Map nodeMap;
    mutable Ioss::Map elemMap;
  };
}

#endif
