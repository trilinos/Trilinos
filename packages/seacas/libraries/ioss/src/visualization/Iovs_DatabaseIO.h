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

#include <Ioss_DatabaseIO.h>
#include <Ioss_EntityType.h>
#include <Ioss_Field.h>
#include <Ioss_DBUsage.h>
#include <Ioss_Map.h>
#include <Ioss_Utils.h>

// with introduction of paraview sierra catalyst plugin, the Iovs stuff is
// always included and NO_PARAVIEWMESH_SUPPORT is never defined.  With the
// plugin architecture, there is no overhead for sierra when the plugin is
// not loaded.  The #define test is left here for now in case developers
// need to use it.
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
#include <iBase.h>
#include <iMesh.h>
#include <iField.h>
#endif

#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include <time.h>

namespace Iovs {

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string& filename,
	       Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
	       const Ioss::PropertyManager &properties);
    ~DatabaseIO();

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const
    { return Ioss::NODEBLOCK|Ioss::ELEMENTBLOCK|Ioss::NODESET|Ioss::SIDESET|Ioss::SIDEBLOCK;}

    /*!
     * Determine the local position of the node with the global id
     * 'global'.  If 'must_exist' is false, then the global id possibly
     * does not exist in the map; otherwise, it must exist and will
     * throw an exception if not found.
     */
    int64_t node_global_to_local(int64_t global, bool must_exist) const
    {return nodeMap.global_to_local(global, must_exist);}
    
    int64_t element_global_to_local(int64_t global) const
    {return elemMap.global_to_local(global);}

    bool begin(Ioss::State state);
    bool   end(Ioss::State state);

    bool begin_state(Ioss::Region *region, int state, double time);
    bool   end_state(Ioss::Region *region, int state, double time);

    void read_meta_data ();

  private:
    // For the time being, treat vis as write only. Consider glue pipelines.
    int64_t get_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::SideBlock* fb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t get_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}

    int64_t put_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			   void *data, size_t data_size) const;

    int64_t put_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::SideBlock* eb, const Ioss::Field& field,
			   void *data, size_t data_size) const { return 0; }

    int64_t put_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int64_t put_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}

    void write_meta_data();

    int64_t handle_node_ids(void* ids, int64_t num_to_get);
    int64_t handle_element_ids(const Ioss::ElementBlock *eb, void* ids, size_t num_to_get);

    const Ioss::Map& get_node_map()         const;
    const Ioss::Map& get_element_map()         const;

    DatabaseIO(); // Do not implement
    DatabaseIO(const DatabaseIO&); // Do not implement
    DatabaseIO& operator=(const DatabaseIO&); // Do not implement

    bool isInput;
    bool singleProcOnly; // True if history or heartbeat which is only written from proc 0...
    bool doLogging; // True if logging field input/output

    // Private member data...
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
    iMesh_Instance mesh_instance; // interface to the vis component
    iBase_EntitySetHandle rootset;
#endif
    // mutable EntityIdSet ids_;

    std::string databaseTitle;
    int spatialDimension;

    int64_t nodeCount;
    int64_t elementCount;

    int nodeBlockCount;
    int elementBlockCount;

    // Bulk Data

    // MAPS -- Used to convert from local exodusII ids/names to Sierra
    // database global ids/names
    // Maps internal (1..num_entity) ids to global ids used on the
    //               sierra side.   global = XXXMap.map[local]
    // XXXMap.map[0] contains: -1 if sequential, 0 if ordering unknown, 1
    // if nonsequential
    mutable Ioss::Map nodeMap;
    mutable Ioss::Map elemMap;
  };
  
};

#endif
