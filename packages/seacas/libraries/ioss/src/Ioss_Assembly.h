// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioss_Assembly_h
#define IOSS_Ioss_Assembly_h

#include "Ioss_EntityType.h" // for EntityType, etc
#include "Ioss_Property.h"   // for Property
#include <Ioss_GroupingEntity.h>
#include <cstddef> // for size_t
#include <cstdint> // for int64_t
#include <string>  // for string
namespace Ioss {
  class DatabaseIO;
  class Field;
} // namespace Ioss

namespace Ioss {

  using EntityContainer = std::vector<const Ioss::GroupingEntity *>;

  /** \brief A homogeneous collection of other GroupingEntities.
   */
  class Assembly : public GroupingEntity
  {
  public:
    Assembly()           = default; // Used for template typing only
    ~Assembly() override = default;
    Assembly(const Assembly &);

    Assembly(DatabaseIO *io_database, const std::string &my_name);

    std::string type_string() const override { return "Assembly"; }
    std::string short_type_string() const override { return "assembly"; }
    std::string contains_string() const override
    {
      return m_members.empty() ? "<EMPTY>" : m_members[0]->type_string();
    }
    EntityType type() const override { return ASSEMBLY; }

    EntityType get_member_type() const { return m_type; }

    bool                   add(const GroupingEntity *member);
    bool                   remove(const GroupingEntity *member);
    const EntityContainer &get_members() const;
    const GroupingEntity * get_member(const std::string &my_name) const;
    void                   remove_members();
    size_t                 member_count() const { return m_members.size(); }

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string &my_name) const override;

  protected:
    int64_t internal_get_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

    int64_t internal_put_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

  private:
    EntityContainer m_members;
    EntityType      m_type{INVALID_TYPE};
  };
} // namespace Ioss
#endif
