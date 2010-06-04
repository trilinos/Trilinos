/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_EdgeBlock_h
#define SIERRA_Ioss_EdgeBlock_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_EntityBlock.h>
#include <string>

namespace Ioss {
  class DatabaseIO;
  class EdgeSet;

  class EdgeBlock : public EntityBlock {
  public:
    friend class EdgeSet;

    EdgeBlock(const DatabaseIO *io_database, const std::string& name,
	      const std::string& edge_type, const std::string& parent_type,
	      size_t edge_count);

    std::string type_string() const {return "EdgeBlock";}
    EntityType type() const {return EDGEBLOCK;}

    const EdgeSet* owner() const {return owner_;}

    void block_membership(std::vector<std::string> &block_members);

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string& name) const;

  protected:
    int internal_get_field_data(const Field& field,
				void *data, size_t data_size) const;

    int internal_put_field_data(const Field& field,
				void *data, size_t data_size) const;
  private:
    // Pointer to the EdgeSet (if any) that contains this edge block.
    const EdgeSet *owner_;
    std::vector<std::string> blockMembership; // What element blocks do the
                                             // elements in this faceset belong to.
  };
}
#endif
