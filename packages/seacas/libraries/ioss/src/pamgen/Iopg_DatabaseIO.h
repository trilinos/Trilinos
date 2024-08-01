// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iopg_export.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_DBUsage.h"    // for DatabaseUsage
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_IOFactory.h"  // for IOFactory
#include "Ioss_Map.h"        // for Map
#include "Ioss_State.h"      // for State
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
  class IOPG_EXPORT IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();

  private:
    IOFactory();
    IOSS_NODISCARD Ioss::DatabaseIO *make_IO(const std::string           &filename,
                                             Ioss::DatabaseUsage          db_usage,
                                             Ioss_MPI_Comm                communicator,
                                             const Ioss::PropertyManager &properties) const;
  };

  class IOPG_EXPORT DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               Ioss_MPI_Comm communicator, const Ioss::PropertyManager &properties);
    ~DatabaseIO();

    IOSS_NODISCARD std::string get_format() const override { return "PamGen"; }

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    IOSS_NODISCARD unsigned entity_field_support() const override
    {
      return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::NODESET | Ioss::SIDESET | Ioss::REGION;
    }

    IOSS_NODISCARD int int_byte_size_db() const override { return 4; }

    IOSS_NODISCARD std::string title() const { return databaseTitle; }
    IOSS_NODISCARD int         maximum_symbol_length() const override { return 32; }

    void compute_block_membership_nl(Ioss::SideBlock *efblock,
                                     Ioss::NameList  &block_membership) const override;

  private:
    void read_meta_data_nl() override;

    bool begin_nl(Ioss::State state) override;
    bool end_nl(Ioss::State state) override;

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

    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    IOSS_NOOP_GFI(Ioss::Region)
    IOSS_NOOP_GFI(Ioss::EdgeBlock)
    IOSS_NOOP_GFI(Ioss::FaceBlock)
    IOSS_NOOP_GFI(Ioss::StructuredBlock)
    IOSS_NOOP_GFI(Ioss::EdgeSet)
    IOSS_NOOP_GFI(Ioss::FaceSet)
    IOSS_NOOP_GFI(Ioss::ElementSet)
    IOSS_NOOP_GFI(Ioss::SideSet)
    IOSS_NOOP_GFI(Ioss::Blob)
    IOSS_NOOP_GFI(Ioss::Assembly)

    // Input only database -- these will never be called...
    IOSS_NOOP_PFI(Ioss::Region)
    IOSS_NOOP_PFI(Ioss::NodeBlock)
    IOSS_NOOP_PFI(Ioss::EdgeBlock)
    IOSS_NOOP_PFI(Ioss::FaceBlock)
    IOSS_NOOP_PFI(Ioss::ElementBlock)
    IOSS_NOOP_PFI(Ioss::StructuredBlock)
    IOSS_NOOP_PFI(Ioss::SideBlock)
    IOSS_NOOP_PFI(Ioss::NodeSet)
    IOSS_NOOP_PFI(Ioss::EdgeSet)
    IOSS_NOOP_PFI(Ioss::FaceSet)
    IOSS_NOOP_PFI(Ioss::ElementSet)
    IOSS_NOOP_PFI(Ioss::SideSet)
    IOSS_NOOP_PFI(Ioss::CommSet)
    IOSS_NOOP_PFI(Ioss::Assembly)
    IOSS_NOOP_PFI(Ioss::Blob)

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
