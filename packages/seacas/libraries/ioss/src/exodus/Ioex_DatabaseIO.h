// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#include <Ioss_DBUsage.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Field.h>
#include <Ioss_Map.h>
#include <Ioss_Utils.h>

#include <exodusII.h>

#include <algorithm>
#include <cstdint>
#include <ctime>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

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
  class StructuredBlock;
  class CommSet;
  class ElementTopology;
} // namespace Ioss

/** \brief A namespace for the exodus database format.
 */
namespace Ioex {
  struct CommunicationMetaData;

  // Used for variable name index mapping
  using VariableNameMap = std::map<std::string, int, std::less<std::string>>;
  using VNMValuePair    = VariableNameMap::value_type;

  // Used to store reduction variables
  using ValueContainer = std::vector<double>;

  // Used for persistent entity IDs
  // The set contains a pair of <ex_entity_type, int>.
  // The ex_entity_type is the exodus entity type defined in
  // exodus's exodusII.h. A couple examples are:
  // EX_ELEM_BLOCK element block and EX_NODE_SET nodeset.
  //
  // The 'int' is the entity id.  The set is used for output databases
  // to ensure that there are no id collisions.
  using EntityIdSet = std::set<std::pair<int64_t, int64_t>>;

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               MPI_Comm communicator, const Ioss::PropertyManager &props);
    DatabaseIO(const DatabaseIO &from) = delete;
    DatabaseIO &operator=(const DatabaseIO &from) = delete;

    ~DatabaseIO() override;

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const override;

  protected:
    // Check to see if database state is ok...
    // If 'write_message' true, then output a warning message indicating the problem.
    // If 'error_message' non-null, then put the warning message into the string and return it.
    // If 'bad_count' non-null, it counts the number of processors where the file does not exist.
    //    if ok returns false, but *bad_count==0, then the routine does not support this argument.
    bool ok__(bool write_message = false, std::string *error_message = nullptr,
              int *bad_count = nullptr) const override = 0;

    // Eliminate as much memory as possible, but still retain meta data information
    // Typically, eliminate the maps...
    void release_memory__() override;

    bool open_group__(const std::string &group_name) override;
    bool create_subgroup__(const std::string &group_name) override;

    bool begin__(Ioss::State state) override;
    bool end__(Ioss::State state) override;

    bool begin_state__(Ioss::Region *region, int state, double time) override;
    bool end_state__(Ioss::Region *region, int state, double time) override;
    void get_step_times__() override = 0;

    int maximum_symbol_length() const override { return maximumNameLength; }

    // NOTE: If this is called after write_meta_data, it will have no affect.
    //       Also, it only affects output databases, not input.
    void set_maximum_symbol_length(int requested_symbol_size) override
    {
      if (!is_input()) {
        maximumNameLength = requested_symbol_size;
      }
    }

    void get_block_adjacencies__(const Ioss::ElementBlock *eb,
                                 std::vector<std::string> &block_adjacency) const override;

    void compute_block_membership__(Ioss::SideBlock *         efblock,
                                    std::vector<std::string> &block_membership) const override;

    int  int_byte_size_db() const override;
    void set_int_byte_size_api(Ioss::DataSize size) const override;

  protected:
    int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override = 0;

    virtual void write_meta_data() = 0;
    void         write_results_metadata();

    void openDatabase__() const override { get_file_pointer(); }

    void closeDatabase__() const override { free_file_pointer(); }

  public:
    // Temporarily made public for use during Salinas transition
    // to using Ioss
    virtual int get_file_pointer() const = 0; // Open file and set exodusFilePtr.
  protected:
    virtual int free_file_pointer() const; // Close file and set exodusFilePtr.

    int  get_current_state() const; // Get current state with error checks and usage message.
    void put_qa();
    void put_info();

    virtual void compute_block_adjacencies() const = 0;

    template <typename T>
    void internal_write_results_metadata(ex_entity_type type, std::vector<T *> entities,
                                         int &glob_index);

    void generate_sideset_truth_table();

    void output_results_names(ex_entity_type type, VariableNameMap &variables) const;
    int gather_names(ex_entity_type type, VariableNameMap &variables,
                     const Ioss::GroupingEntity *ge, int index, bool reduction);

    void get_nodeblocks();

    void add_attribute_fields(ex_entity_type entity_type, Ioss::GroupingEntity *block,
                              int attribute_count, const std::string &type);

    void output_other_meta_data();

    int64_t internal_add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity,
                                        int64_t position, int64_t block_count,
                                        Ioss::IntVector &      truth_table,
                                        Ioex::VariableNameMap &variables);
    int64_t add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity,
                               int64_t position = 0);

    void add_region_fields();
    void store_reduction_field(ex_entity_type type, const Ioss::Field &field,
                               const Ioss::GroupingEntity *ge, void *variables) const;

    void get_reduction_field(ex_entity_type type, const Ioss::Field &field,
                             const Ioss::GroupingEntity *ge, void *variables) const;
    void write_reduction_fields() const;
    void read_reduction_fields() const;

    // Handle special output time requests -- primarily restart (cycle, keep, overwrite)
    // Given the global region step, return the step on the database...
    int get_database_step(int global_step) const;

    void flush_database__() const override;
    void finalize_write(double sim_time);

    // Private member data...
  protected:
    mutable int         exodusFilePtr;
    mutable std::string m_groupName;

    mutable EntityIdSet ids_;

    mutable int exodusMode;
    mutable int dbRealWordSize;

    mutable int maximumNameLength;
    int         spatialDimension;

    int64_t nodeCount;
    int64_t edgeCount;
    int64_t faceCount;
    int64_t elementCount;

    mutable std::map<ex_entity_type, int> m_groupCount;

    // Communication Set Data
    Ioss::Int64Vector nodeCmapIds;
    Ioss::Int64Vector nodeCmapNodeCnts;
    Ioss::Int64Vector elemCmapIds;
    Ioss::Int64Vector elemCmapElemCnts;
    int64_t           commsetNodeCount;
    int64_t           commsetElemCount;

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

    mutable std::map<ex_entity_type, Ioss::IntVector> m_truthTable;
    mutable std::map<ex_entity_type, VariableNameMap> m_variables;

    mutable ValueContainer globalValues;

    mutable std::vector<std::vector<bool>> blockAdjacency;
    mutable std::vector<unsigned char>     nodeConnectivityStatus;

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
} // namespace Ioex
#endif
