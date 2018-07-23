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
#ifndef IOSS_Iofx_DatabaseIO_h
#define IOSS_Iofx_DatabaseIO_h

#include <Ioss_DBUsage.h>
#include <Ioss_Field.h>
#include <Ioss_Map.h>
#include <Ioss_Utils.h>
#include <exodus/Ioex_DatabaseIO.h>

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

namespace Ioex {
  struct CommunicationMetaData;
} // namespace Ioex

/** \brief A namespace for the file-per-process version of the
 *  parallel exodus database format.
 */
namespace Iofx {
  class DatabaseIO : public Ioex::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               MPI_Comm communicator, const Ioss::PropertyManager &props);
    DatabaseIO(const DatabaseIO &from) = delete;
    DatabaseIO &operator=(const DatabaseIO &from) = delete;
    ~DatabaseIO() override{};

  private:
    // Check to see if database state is ok...
    // If 'write_message' true, then output a warning message indicating the problem.
    // If 'error_message' non-null, then put the warning message into the string and return it.
    // If 'bad_count' non-null, it counts the number of processors where the file does not exist.
    //    if ok returns false, but *bad_count==0, then the routine does not support this argument.
    bool ok__(bool write_message = false, std::string *error_msg = nullptr,
              int *bad_count = nullptr) const override;

    void get_step_times__() override;

    bool open_input_file(bool write_message, std::string *error_msg, int *bad_count,
                         bool abort_if_error) const;
    bool handle_output_file(bool write_message, std::string *error_msg, int *bad_count,
                            bool overwrite, bool abort_if_error) const;
    bool check_valid_file_ptr(bool write_message, std::string *error_msg, int *bad_count,
                              bool abort_if_error) const;

    int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override
    {
      return -1;
    }

    int64_t get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override
    {
      return -1;
    }
    int64_t put_Xset_field_internal(ex_entity_type type, const Ioss::EntitySet *ns,
                                    const Ioss::Field &field, void *data, size_t data_size) const;
    int64_t get_Xset_field_internal(ex_entity_type type, const Ioss::EntitySet *ns,
                                    const Ioss::Field &field, void *data, size_t data_size) const;

  public:
    // Temporarily made public for use during Salinas transition
    // to using Ioss
    int get_file_pointer() const override; // Open file and set exodusFilePtr.
  private:
    int64_t read_nodal_coordinates();
    void    read_elements(const Ioss::ElementBlock &block);

    void compute_node_status() const;

    // Metadata-related functions.
    void read_meta_data__() override;
    void read_communication_metadata();

    int64_t read_transient_field(ex_entity_type type, const Ioex::VariableNameMap &variables,
                                 const Ioss::Field &field, const Ioss::GroupingEntity *ge,
                                 void *data) const;

    int64_t read_attribute_field(ex_entity_type type, const Ioss::Field &field,
                                 const Ioss::GroupingEntity *ge, void *data) const;

    int64_t write_attribute_field(ex_entity_type type, const Ioss::Field &field,
                                  const Ioss::GroupingEntity *ge, void *data) const;

    // Handles subsetting of side blocks.
    int64_t read_ss_transient_field(const Ioss::Field &field, int64_t id, void *variables,
                                    std::vector<int> &is_valid_side) const;

    // Should be made more generic again so can rejoin with write_element_transient field
    void write_nodal_transient_field(ex_entity_type type, const Ioss::Field &field,
                                     const Ioss::NodeBlock *ge, int64_t count,
                                     void *variables) const;
    // Should be made more generic again so can rejoin with write_nodal_transient field
    void write_entity_transient_field(ex_entity_type type, const Ioss::Field &field,
                                      const Ioss::GroupingEntity *ge, int64_t count,
                                      void *variables) const;
    void write_meta_data() override;
    void gather_communication_metadata(Ioex::CommunicationMetaData *meta);

    // Read related metadata and store it in the region...
    void read_region();
    void get_edgeblocks();
    void get_faceblocks();
    void get_elemblocks();
    void get_blocks(ex_entity_type entity_type, int rank_offset, const std::string &basename);

    void get_sidesets();

    template <typename T>
    void get_sets(ex_entity_type type, int64_t count, const std::string &base,
                  const T * /*unused*/);
    void get_nodesets();
    void get_edgesets();
    void get_facesets();
    void get_elemsets();

    void get_commsets();

    // ID Mapping functions.
    const Ioss::Map &get_map(ex_entity_type type) const;
    const Ioss::Map &get_map(Ioss::Map &entity_map, int64_t entityCount, ex_entity_type entity_type,
                             ex_inquiry inquiry_type) const;

    // Internal data handling
    int64_t handle_node_ids(void *ids, int64_t num_to_get) const;
    int64_t handle_element_ids(const Ioss::ElementBlock *eb, void *ids, size_t num_to_get) const;
    int64_t handle_face_ids(const Ioss::FaceBlock *eb, void *ids, size_t num_to_get) const;
    int64_t handle_edge_ids(const Ioss::EdgeBlock *eb, void *ids, size_t num_to_get) const;

    int64_t get_side_connectivity(const Ioss::SideBlock *fb, int64_t id, int64_t my_side_count,
                                  void *fconnect, bool map_ids) const;
    template <typename INT>
    int64_t get_side_connectivity_internal(const Ioss::SideBlock *fb, int64_t id,
                                           int64_t side_count, INT *fconnect, bool map_ids) const;
    int64_t get_side_distributions(const Ioss::SideBlock *fb, int64_t id, int64_t my_side_count,
                                   double *dist_fact, size_t data_size) const;

    int64_t get_side_field(const Ioss::SideBlock *ef_blk, const Ioss::Field &field, void *data,
                           size_t data_size) const;
    int64_t put_side_field(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                           size_t data_size) const;

    mutable bool isSerialParallel{
        false}; //!< true if application code is controlling the processor id.
  };
} // namespace Iofx
#endif
