/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_EntityBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>

#include <string>
#include <assert.h>

Ioss::EntityBlock::EntityBlock(const Ioss::DatabaseIO *io_database,
			       const std::string& my_name,
			       const std::string& entity_type,
			       const std::string& parent_topology_type,
			       size_t entity_count)
  : Ioss::GroupingEntity(io_database, my_name), parentTopology_(NULL),
    parentElementBlock_(NULL), consistentSideNumber(-1)
{
  topology_ = ElementTopology::factory(entity_type);
  assert(topology_ != NULL);

  parentTopology_ = ElementTopology::factory(parent_topology_type);
  assert(parentTopology_ != NULL);

  properties.add(Ioss::Property("entity_count", static_cast<int>(entity_count)));
  properties.add(Ioss::Property(this, "topology_node_count",
			       Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this, "topology_type",
			       Ioss::Property::STRING));
  properties.add(Ioss::Property(this, "parent_topology_type",
			       Ioss::Property::STRING));

  fields.add(Ioss::Field("ids",
			Ioss::Field::INTEGER, "scalar",
			Ioss::Field::MESH, entity_count));

  fields.add(Ioss::Field("connectivity", Ioss::Field::INTEGER,
			topology_->name(),
			Ioss::Field::MESH, entity_count));
}

Ioss::Property Ioss::EntityBlock::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "topology_node_count")
    return Ioss::Property(my_name, topology()->number_nodes());
  else if (my_name == "topology_type")
    return Ioss::Property(my_name, topology()->name());
  else if (my_name == "parent_topology_type")
    return Ioss::Property(my_name, parent_element_topology()->name());
  else
    return Ioss::GroupingEntity::get_implicit_property(my_name);
}

int  Ioss::EntityBlock::get_consistent_side_number() const
{
  if (consistentSideNumber == -1) {
    // It wasn't calculated during the metadata reading of the surfaces.
    // Determine it now...
    int ecount  = get_property("entity_count").get_int();
    std::vector<int> element_side(2*ecount);
    if (field_exists("element_side")) {
      get_field_data("element_side", element_side);
      
      int side = ecount > 0 ? element_side[1] : 0;
      for (int i=3; i < 2*ecount; i+=2) {
	int this_side = element_side[i];
	if (this_side != side) {
	  side = 999; // Indicates the sides are not consistent ;
	  break;
	}
      }
      int side_max = get_database()->util().global_minmax(side, Ioss::ParallelUtils::DO_MAX);
      if (side_max != 999)
	consistentSideNumber = side_max;
      else
	consistentSideNumber = 0;
    } else {
      consistentSideNumber = 0;
    }
  }
  return consistentSideNumber;
}
