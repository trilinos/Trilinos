// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioss_SideBlock_h
#define IOSS_Ioss_SideBlock_h

#include <Ioss_ElementBlock.h>
#include <Ioss_EntityBlock.h> // for EntityBlock
#include <Ioss_EntityType.h>  // for EntityType, etc
#include <Ioss_Property.h>    // for Property
#include <Ioss_SideSet.h>
#include <cstddef> // for size_t
#include <cstdint> // for int64_t
#include <string>  // for string
#include <vector>  // for vector
namespace Ioss {
  class DatabaseIO;
} // namespace Ioss
namespace Ioss {
  class ElementTopology;
} // namespace Ioss
namespace Ioss {
  class Field;
} // namespace Ioss

namespace Ioss {

  /** \brief A collection of element sides having the same topology.
   */
  class SideBlock : public EntityBlock
  {
  public:
    friend class SideSet;

    SideBlock(DatabaseIO *io_database, const std::string &my_name, const std::string &side_type,
              const std::string &element_type, size_t side_count);

    SideBlock(const SideBlock &other);

    std::string type_string() const override { return "SideBlock"; }
    std::string short_type_string() const override { return "sideblock"; }
    std::string contains_string() const override { return "Element/Side pair"; }
    EntityType  type() const override { return SIDEBLOCK; }

    const SideSet *             owner() const { return owner_; }
    const Ioss::GroupingEntity *contained_in() const override { return owner_; }

    void block_membership(std::vector<std::string> &block_members) override;

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string &my_name) const override;

    // For faceblock, edgeblock, if they are split by element block, then this
    // will be non-nullptr and is a pointer to the parent element block for this
    // faceblock or edgeblock. Has no meaning for other EntityBlock types or split
    // types.
    const ElementBlock *parent_element_block() const
    {
      return dynamic_cast<const ElementBlock *>(parentBlock_);
    }

    void set_parent_element_block(const ElementBlock *element_block)
    {
      parentBlock_ = element_block;
    }

    const EntityBlock *parent_block() const { return parentBlock_; }
    void               set_parent_block(const EntityBlock *block) { parentBlock_ = block; }

    // Describes the contained entities element block topology
    const ElementTopology *parent_element_topology() const { return parentTopology_; }

    // For faceblock, edgeblock, return whether the surface is applied
    // to the same face/edge for all elements in the surface. If not,
    // return 0; otherwise return the consistent face number.
    int  get_consistent_side_number() const;
    void set_consistent_side_number(int side) { consistentSideNumber = side; }

    bool operator==(const SideBlock &) const;
    bool operator!=(const SideBlock &) const;
    bool equal(const SideBlock &) const;

  protected:
    int64_t internal_get_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

    int64_t internal_put_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

  private:
    bool equal_(const SideBlock &, bool quiet) const;

    const SideSet *    owner_{nullptr};
    ElementTopology *  parentTopology_{nullptr}; // Topology of parent element (if any)
    const EntityBlock *parentBlock_{nullptr};

    // Pointer to the SideSet (if any) that contains this side block.
    std::vector<std::string> blockMembership{}; // What element blocks do the
                                                // elements in this sideset belong to.
    mutable int consistentSideNumber{-1};
  };
} // namespace Ioss
#endif
