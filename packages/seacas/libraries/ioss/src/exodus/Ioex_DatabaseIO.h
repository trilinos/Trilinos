// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// -*- Mode: c++ -*-
#ifndef IOSS_Ioex_DatabaseIO_h
#define IOSS_Ioex_DatabaseIO_h

#include <Ioss_DatabaseIO.h>
#include <Ioss_Field.h>
#include <Ioss_DBUsage.h>
#include <Ioss_Map.h>
#include <Ioss_Utils.h>

#include <exodusII.h>

#include <stdint.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include <time.h>

namespace Ioss {
  class GroupingEntity;
  class Region;
  class EntityBlock;
  class NodeBlock;
  class EdgeBlock;
  class FaceBlock;
  class ElementBlock;
  class EntitySet;
  class NodeSet;
  class EdgeSet;
  class FaceSet;
  class ElementSet;
  class SideBlock;
  class SideSet;
  class CommSet;
  class ElementTopology;
}

namespace Ioex {
  struct CommunicationMetaData;

  // Used for variable name index mapping
  typedef std::map<std::string, int, std::less<std::string> > VariableNameMap;
  typedef VariableNameMap::value_type VNMValuePair;

  // Used to store reduction variables
  typedef std::vector<double> ValueContainer;

  // Used for persistent entity IDs
  // The set contains a pair of <ex_entity_type, int>.
  // The ex_entity_type is the exodus entity type defined in
  // exodus's exodusII.h. A couple examples are:
  // EX_ELEM_BLOCK element block and EX_NODE_SET nodeset.
  //
  // The 'int' is the entity id.  The set is used for output databases
  // to ensure that there are no id collisions.
  typedef std::set<std::pair<int64_t, int64_t> > EntityIdSet;

  class DatabaseIO : public Ioss::DatabaseIO
  {
    public:
      DatabaseIO(Ioss::Region *region, const std::string& filename,
                 Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                 const Ioss::PropertyManager &properties);
      virtual ~DatabaseIO();

      // Check to see if database state is ok...
      // If 'write_message' true, then output a warning message indicating the problem.
      // If 'error_message' non-null, then put the warning message into the string and return it.
      // If 'bad_count' non-null, it counts the number of processors where the file does not exist.
      //    if ok returns false, but *bad_count==0, then the routine does not support this argument.
      virtual bool ok(bool write_message = false, std::string *error_message=NULL, int *bad_count=NULL) const = 0;

      // Eliminate as much memory as possible, but still retain meta data information
      // Typically, eliminate the maps...
      void release_memory();

      // Check capabilities of input/output database...  Returns an
      // unsigned int with the supported Ioss::EntityTypes or'ed
      // together. If "return_value & Ioss::EntityType" is set, then the
      // database supports that type (e.g. return_value & Ioss::FACESET)
      unsigned entity_field_support() const;

      bool open_group(const std::string &group_name);
      bool create_subgroup(const std::string &group_name);

      bool begin(Ioss::State state);
      bool   end(Ioss::State state);

      bool begin_state(Ioss::Region *region, int state, double time);
      bool   end_state(Ioss::Region *region, int state, double time);
      virtual void get_step_times() = 0;
      
      int    spatial_dimension()   const     {return spatialDimension;}
      int64_t    node_count()          const     {return nodeCount;}
      int64_t    side_count()          const     {return 0;}
      int64_t    element_count()       const     {return elementCount;}
      int    node_block_count()    const     {return m_groupCount[EX_NODE_BLOCK];}
      int    element_block_count() const     {return m_groupCount[EX_ELEM_BLOCK];}
      int    sideset_count()       const     {return m_groupCount[EX_SIDE_SET];}
      int    nodeset_count()       const     {return m_groupCount[EX_NODE_SET];}
      int    maximum_symbol_length() const   {return maximumNameLength;}

      void get_block_adjacencies(const Ioss::ElementBlock *eb,
                                 std::vector<std::string> &block_adjacency) const;

      virtual void compute_block_membership(int64_t id, std::vector<std::string> &block_membership) const = 0;
      void compute_block_membership(Ioss::SideBlock *efblock,
                                    std::vector<std::string> &block_membership) const;

  protected:
      int64_t get_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
                                 void *data, size_t data_size) const;
      virtual int64_t get_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::SideBlock* fb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t get_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;

      int64_t put_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
                                 void *data, size_t data_size) const;
      virtual int64_t put_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::SideBlock* fb, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;
      virtual int64_t put_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
					 void *data, size_t data_size) const = 0;

      virtual void write_meta_data() = 0;
      void write_results_metadata();


      // Private member functions
      DatabaseIO(const DatabaseIO& from); // do not implement
      DatabaseIO& operator=(const DatabaseIO& from); // do not implement

      virtual void openDatabase() const {
        get_file_pointer();
      }

      virtual void closeDatabase() const {
        free_file_pointer();
      }

  public:
      // Temporarily made public for use during Salinas transition
      // to using Ioss
      virtual int get_file_pointer() const = 0; // Open file and set exodusFilePtr.
  protected:
      int free_file_pointer() const; // Close file and set exodusFilePtr.

      int get_current_state() const; // Get current state with error checks and usage message.
      void put_qa();
      void put_info();

      virtual void compute_block_adjacencies() const = 0;

      template <typename T>
	void internal_write_results_metadata(ex_entity_type type,
					     std::vector<T*> entities,
					     int &glob_index);

      void generate_sideset_truth_table();

      void output_results_names(ex_entity_type type, VariableNameMap &variables) const;
      int  gather_names(ex_entity_type type,
                        VariableNameMap &variables,
                        const Ioss::GroupingEntity *ge,
                        int index, bool reduction);

      void get_nodeblocks();

      void add_attribute_fields(ex_entity_type ent_type, Ioss::GroupingEntity *block,
                                int attribute_count,  const std::string& type);

      void output_other_meta_data();

      int64_t internal_add_results_fields(ex_entity_type type,
                                          Ioss::GroupingEntity *entity,
                                          int64_t position, int64_t block_count,
                                          Ioss::IntVector &truth_table,
                                          Ioex::VariableNameMap &variables);
      int64_t add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity, int64_t position=0);

      void add_region_fields();
      void store_reduction_field(ex_entity_type type,
                                 const Ioss::Field& field,
                                 const Ioss::GroupingEntity *ge,
                                 void *variables) const;

      void get_reduction_field(ex_entity_type type,
                               const Ioss::Field& field,
                               const Ioss::GroupingEntity *ge,
                               void *variables) const;
      void write_reduction_fields() const;
      void read_reduction_fields() const;

      // Handle special output time requests -- primarily restart (cycle, keep, overwrite)
      // Given the global region step, return the step on the database...
      int get_database_step(int global_step) const;

      void finalize_write(double sim_time);


      // Private member data...
  protected:
      mutable int exodusFilePtr;
      mutable std::string m_groupName;
      
      mutable EntityIdSet ids_;

      mutable int exodusMode;
      mutable int dbRealWordSize;

      mutable int maximumNameLength;
      int spatialDimension;

      int64_t nodeCount;
      int64_t edgeCount;
      int64_t faceCount;
      int64_t elementCount;

      mutable std::map<ex_entity_type,int> m_groupCount;

      // Communication Set Data
      Ioss::Int64Vector nodeCmapIds;
      Ioss::Int64Vector nodeCmapNodeCnts;
      Ioss::Int64Vector elemCmapIds;
      Ioss::Int64Vector elemCmapElemCnts;
      int64_t commsetNodeCount;
      int64_t commsetElemCount;

      // Bulk Data

      // MAPS -- Used to convert from local exodusII ids/names to Sierra
      // database global ids/names

      //---Node Map -- Maps internal (1..NUMNP) ids to global ids used on the
      //               sierra side.   global = nodeMap[local]
      // nodeMap[0] contains: -1 if sequential, 0 if ordering unknown, 1
      // if nonsequential

      mutable Ioss::Map nodeMap;
      mutable Ioss::Map edgeMap;
      mutable Ioss::Map faceMap;
      mutable Ioss::Map elemMap;

      // --- Nodal/Element/Attribute Variable Names -- Maps from sierra
      // field names to index of nodal/element/attribute variable in
      // exodusII. Note that the component suffix of the field is added on
      // prior to searching the map for the index.  For example, given the
      // Sierra field 'displ' which is a VECTOR_3D, the names stored in
      // 'elementMap' would be 'displ_x', 'displ_y' and 'displ_z'.  All
      // names are converted to lowercase.

      mutable std::map<ex_entity_type,Ioss::IntVector> m_truthTable;
      mutable std::map<ex_entity_type,VariableNameMap> m_variables;

      mutable ValueContainer  globalValues;

      mutable std::vector<std::vector<bool> > blockAdjacency;
      mutable std::vector<unsigned char> nodeConnectivityStatus;

      // For a database with omitted blocks, this map contains the indices of the
      // active nodes for each nodeset.  If the nodeset is not reduced in size,
      // the map's vector will be empty for that nodeset. If the vector is not
      // empty, then some nodes on that nodeset are only connected to omitted elements.
      mutable std::map<std::string, Ioss::Int64Vector> activeNodesetNodesIndex;
	
      time_t timeLastFlush;

      mutable bool fileExists; // False if file has never been opened/created
      mutable bool minimizeOpenFiles;

      mutable bool blockAdjacenciesCalculated; // True if the lazy creation of
      // block adjacencies has been calculated.
      mutable bool nodeConnectivityStatusCalculated; // True if the lazy creation of
      // nodeConnectivityStatus has been calculated.
  };
}
#endif
