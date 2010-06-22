/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_EdgeSet_h
#define IOSS_Ioss_EdgeSet_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_GroupingEntity.h>
#include <string>

#include <vector>

namespace Ioss {
  class DatabaseIO;
  class EdgeBlock;
  class EntityBlock;
  
  typedef std::vector<EdgeBlock*> EdgeBlockContainer;

  class EdgeSet : public GroupingEntity {
  public:
    EdgeSet(const DatabaseIO *io_database, const std::string& name);
    ~EdgeSet();

    std::string type_string() const {return "EdgeSet";}
    EntityType type() const {return EDGESET;}

    bool add(EdgeBlock    *edge_block);
    const EdgeBlockContainer&    get_edge_blocks() const;
    EdgeBlock*    get_edge_block(const std::string& name) const;
    size_t edge_block_count() const {return edgeBlocks.size();}

    size_t block_count() const {return edgeBlocks.size();}
    EntityBlock* get_block(size_t which) const;

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
    EdgeBlockContainer      edgeBlocks;
    std::vector<std::string> blockMembership; // What element blocks do the
						 // elements in this edgeset belong to.
  };
}
#endif
