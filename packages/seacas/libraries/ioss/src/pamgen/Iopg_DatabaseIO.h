// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Iopg_DatabaseIO_h
#define IOSS_Iopg_DatabaseIO_h

#include "Ioss_State.h" // for State
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>    // for DatabaseUsage
#include <Ioss_DatabaseIO.h> // for DatabaseIO
#include <Ioss_IOFactory.h>  // for IOFactory
#include <Ioss_Map.h>        // for Map
#include <stddef.h>          // for size_t
#include <stdint.h>          // for int64_t
#include <string>            // for string
#include <vector>            // for vector

namespace Ioss {
  class Assembly;
  class Blob;
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
namespace Iopg {
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

    const std::string get_format() const override { return "PamGen"; }

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const override
    {
      return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::NODESET | Ioss::SIDESET | Ioss::REGION;
    }

    int int_byte_size_db() const override { return 4; }

    std::string title() const { return databaseTitle; }
    int         maximum_symbol_length() const override { return 32; }

    void compute_block_membership__(Ioss::SideBlock *         efblock,
                                    std::vector<std::string> &block_membership) const override;

  private:
    void read_meta_data__() override;

    bool begin__(Ioss::State state) override;
    bool end__(Ioss::State state) override;

    void read_region();
    void read_communication_metadata();

    void get_nodeblocks();
    void get_elemblocks();
    void get_nodesets();
    void get_sidesets();
    void get_commsets();

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
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override
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
    int64_t put_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return 0;
    }
    int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override
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

    std::string databaseTitle;

    int spatialDimension{3};

    int nodeBlockCount{0};
    int elementBlockCount{0};
    int nodesetCount{0};
    int sidesetCount{0};

    // Communication Set Data
    Ioss::IntVector nodeCmapIds;
    Ioss::IntVector nodeCmapNodeCnts;
    Ioss::IntVector elemCmapIds;
    Ioss::IntVector elemCmapElemCnts;
    int             commsetNodeCount{0};
    int             commsetElemCount{0};
  };
} // namespace Iopg
#endif // IOSS_Iopg_DatabaseIO_h
