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
    int64_t iovs_internal_node_global_to_local(int64_t global, bool must_exist) const;
    int64_t    node_global_to_local(int64_t global, bool must_exist) const;
    // this is unnecessary for vis
    int64_t element_global_to_local(int64_t global) const { return global; }

    bool begin(Ioss::State state);
    bool   end(Ioss::State state);

    bool begin_state(Ioss::Region *region, int state, double time);
    bool   end_state(Ioss::Region *region, int state, double time);

    void read_meta_data ();

    void compute_block_membership(int id, std::vector<std::string> &block_membership) const;
    void compute_block_membership(Ioss::EntityBlock *efblock,
				  std::vector<std::string> &block_membership) const;

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

    void build_element_reorder_map(int64_t start, int64_t count);
    void build_node_reorder_map(void *new_ids, int64_t count);

    const Ioss::MapContainer& get_node_map()         const;
    const Ioss::MapContainer& get_element_map()         const;

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
    int nodesetCount;
    int sidesetCount;

    // Communication Set Data
    int *nodeCmapIds;
    int *nodeCmapNodeCnts;
    int *elemCmapIds;
    int *elemCmapElemCnts;
    int commsetNodeCount;
    int commsetElemCount;

    // Bulk Data

    // MAPS -- Used to convert from local exodusII ids/names to Sierra
    // database global ids/names

    //---Node Map -- Maps internal (1..NUMNP) ids to global ids used on the
    //               sierra side.   global = nodeMap[local]
    // nodeMap[0] contains: -1 if sequential, 0 if ordering unknown, 1
    // if nonsequential
    mutable Ioss::MapContainer        nodeMap;
    mutable Ioss::MapContainer        reorderNodeMap;
    mutable Ioss::ReverseMapContainer reverseNodeMap;
    // (local==global)

    //---Element Map -- Maps internal (1..NUMEL) ids to global ids used on the
    //               sierra side.   global = elementMap[local]
    // elementMap[0] contains: -1 if sequential, 0 if ordering unknown,
    // 1 if nonsequential
     mutable Ioss::MapContainer        elementMap;
     mutable Ioss::MapContainer        reorderElementMap;
     mutable Ioss::ReverseMapContainer reverseElementMap;
    mutable bool sequentialEG2L; // true if reverse element map is
    // sequential (local==global)

    // --- Nodal/Element/Attribute Variable Names -- Maps from sierra
    // field names to index of nodal/element/attribute variable in
    // exodusII. Note that the component suffix of the field is added on
    // prior to searching the map for the index.  For example, given the
    // Sierra field 'displ' which is a VECTOR_3D, the names stored in
    // 'elementMap' would be 'displ_x', 'displ_y' and 'displ_z'.  All
    // names are converted to lowercase.

    // VariableNameMap nodalVariables;
    // VariableNameMap elementVariables;
    // VariableNameMap elementAttributes;
    // VariableNameMap nodesetVariables;
    // VariableNameMap sidesetVariables;
    // VariableNameMap globalVariables;

    // Similar for entity attributes (currently only element blocks); however, the
    // names are not mapped to the individual components.
    // The other difference is that that name stored in the map is a
    // combination of the entity-block name and the attribute name
    // so we can use a single VariableNameMap for the region instead
    // of one-per entity block.
    // VariableNameMap attributeNames;

    // mutable ValueContainer  globalValues;

    mutable bool fileExists; // False if file has never been opened/created
    mutable bool minimizeOpenFiles;

    mutable bool blockAdjacenciesCalculated; // True if the lazy creation of
    // block adjacencies has been calculated.
    mutable std::vector<std::vector<bool> > blockAdjacency;

    int exodusMode;
    time_t timeLastFlush;
  };
  
  // ------------------------------------------------------------------------
  // Node and Element mapping functions.  The ExodusII database
  // stores ids in a local-id system (1..NUMNP), (1..NUMEL) but
  // Sierra wants entities in a global system. These routines
  // take care of the mapping from local <-> global

  typedef std::vector<Ioss::IdPair>::iterator RMapI;
  inline int64_t DatabaseIO::node_global_to_local(int64_t global, bool must_exist) const
    {
      return iovs_internal_node_global_to_local(global, must_exist);
      /*
      if (nodeMap.empty()) {
	get_node_map();
      }
      int64_t local = global;
      if (nodeMap[0] != -1) {
	std::pair<RMapI, RMapI> iter = std::equal_range(reverseNodeMap.begin(),
							reverseNodeMap.end(),
							global,
							Ioss::IdPairCompare());
	if (iter.first != iter.second)
	  local = (iter.first)->second;
	else
	  local = 0;
	if (must_exist && iter.first == iter.second) {
	  std::ostringstream errmsg;
	  errmsg << "Node with global id equal to " << global
		 << " does not exist in this mesh on this processor\n";
	  IOSS_ERROR(errmsg);
	}
      } else if (!must_exist && global > nodeCount) {
	local = 0;
      }
      if (local > nodeCount || (local <= 0 && must_exist)) {
	std::ostringstream errmsg;
	errmsg << "Node with global id equal to " << global
	       << " returns a local id of " << local
	       << " which is invalid. This should not happen, please report.\n";
	IOSS_ERROR(errmsg);
      }
      return local;
      */
    }

  /*
    inline int DatabaseIO::element_global_to_local(int global) const
  */
};

#endif
