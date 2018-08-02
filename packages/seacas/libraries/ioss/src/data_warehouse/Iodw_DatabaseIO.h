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

#ifndef Iodw_DatabaseIO_h
#define Iodw_DatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>    // for DatabaseUsage
#include <Ioss_DatabaseIO.h> // for DatabaseIO
#include <Ioss_IOFactory.h>  // for IOFactory
#include <Ioss_Map.h>        // for Map
#include <Ioss_State.h>      // for State
#include <cstddef>           // for size_t
#include <cstdint>           // for int64_t
#include <string>            // for string
#include <vector>            // for vector

namespace Ioss {
  class CommSet;
  class EdgeBlock;
  class EdgeSet;
  class ElementBlock;
  class ElementSet;
  class EntityBlock;
  class FaceBlock;
  class FaceSet;
  class Field;
  class GroupingEntity;
  class NodeBlock;
  class NodeSet;
  class PropertyManager;
  class Region;
  class SideBlock;
  class SideSet;
  class StructuredBlock;
} // namespace Ioss

/** \brief A namespace for the pamgen database format.
 */
namespace Iodw {
  class IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              MPI_Comm communicator, const Ioss::PropertyManager &properties) const;
  };

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               MPI_Comm communicator, const Ioss::PropertyManager &properties);
    DatabaseIO(const DatabaseIO &from) = delete;
    DatabaseIO &operator=(const DatabaseIO &from) = delete;
    ~DatabaseIO();

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const override { return 0; }

    std::string title() const { return databaseTitle; }
    int         spatial_dimension() const { return spatialDimension; }
    int         node_count() const { return nodeCount; }
    int         side_count() const { return 0; }
    int         element_count() const { return elementCount; }
    int         node_block_count() const { return nodeBlockCount; }
    int         element_block_count() const { return elementBlockCount; }
    int         sideset_count() const { return sidesetCount; }
    int         nodeset_count() const { return nodesetCount; }
    int         maximum_symbol_length() const override { return 32; }

    void compute_block_membership(Ioss::SideBlock *         efblock,
                                  std::vector<std::string> &block_membership) const;

  private:
    void read_meta_data__() override;

    bool begin__(Ioss::State state) override { return false; };
    bool end__(Ioss::State state) override { return false; };

    bool begin_state__(Ioss::Region *region, int state, double time) override { return false; };
    bool end_state__(Ioss::Region *region, int state, double time) override { return false; };

    void read_region();
    void read_communication_metadata();

    /*
     * TODO identify all the get_*{blocks|sets} needed here
     */
    void get_edgeblocks();
    void get_elemblocks();
    void get_faceblocks();
    void get_nodeblocks();

    void get_edgesets();
    void get_elemsets();
    void get_facesets();
    void get_nodesets();

    int get_side_connectivity(const Ioss::SideBlock *fb, int id, int side_count, int *fconnect,
                              size_t data_size) const;
    int get_side_distributions(const Ioss::SideBlock *fb, int id, int side_count, double *dist_fact,
                               size_t data_size) const;

    const Ioss::Map &get_node_map() const;
    const Ioss::Map &get_element_map() const;

    int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
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
    int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override;
    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
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
                               void *data, size_t data_size) const override;

    std::string databaseTitle;

    int spatialDimension;

    int nodeBlockCount;
    int elementBlockCount;
    int nodesetCount;
    int sidesetCount;

    // Communication Set Data
    Ioss::IntVector nodeCmapIds;
    Ioss::IntVector nodeCmapNodeCnts;
    Ioss::IntVector elemCmapIds;
    Ioss::IntVector elemCmapElemCnts;
    int             commsetNodeCount;
    int             commsetElemCount;
  };
} // namespace Iodw
#endif // Iodw_DatabaseIO_h
