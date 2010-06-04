/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_EntityBlock_h
#define SIERRA_Ioss_EntityBlock_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_GroupingEntity.h>
#include <string>

namespace Ioss {
  class ElementTopology;
  class ElementBlock;
  class DatabaseIO;

  class EntityBlock : public GroupingEntity {
  public:
    virtual Property
      get_implicit_property(const std::string& name) const = 0;

    // Describes the contained entities topology
    const ElementTopology *topology() const {return topology_;}

    // Describes the contained entities element block topology
    const ElementTopology *parent_element_topology() const
      { return parentTopology_;}

    // For faceblock, edgeblock, if they are split by element block, then this
    // will be non-NULL and is a pointer to the parent element block for this
    // faceblock or edgeblock. Has no meaning for other EntityBlock types or split
    // types.
    const ElementBlock *parent_element_block() const
      { return parentElementBlock_;}
    void set_parent_element_block(ElementBlock *element_block)
    { parentElementBlock_ = element_block; }
    
    // For faceblock, edgeblock, return whether the surface is applied
    // to the same face/edge for all elements in the surface. If not,
    // return 0; otherwise return the consistent face number.
    int  get_consistent_side_number() const;
    void set_consistent_side_number(int side) {consistentSideNumber = side;}

  protected:
    EntityBlock(const DatabaseIO *io_database,
		const std::string& name, const std::string& entity_type,
		const std::string& parent_topology_type,
		size_t entity_count);

  private:
    EntityBlock(const EntityBlock&); // do not implement
    EntityBlock& operator=(const EntityBlock&); // do not implement
    ElementTopology *topology_;
    ElementTopology *parentTopology_; // Topology of parent element (if any)
    ElementBlock    *parentElementBlock_;
  protected:
    mutable int consistentSideNumber;
  };
}
#endif
