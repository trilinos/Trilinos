/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_EdgeBlock.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_ElementTopology.h>

#include <string>

Ioss::EdgeBlock::EdgeBlock(const Ioss::DatabaseIO *io_database,
			   const std::string& my_name,
			   const std::string& edge_type,
			   const std::string& parent_type,
			   size_t edge_count)
  : Ioss::EntityBlock(io_database, my_name, edge_type, parent_type, edge_count), owner_(NULL)
{
  properties.add(Ioss::Property(this, "distribution_factor_count",
			       Ioss::Property::INTEGER));

  fields.add(Ioss::Field("element_side",
			Ioss::Field::INTEGER, "pair",
			Ioss::Field::MESH, edge_count));

  // Distribution factors are optional...
}

int Ioss::EdgeBlock::internal_get_field_data(const Ioss::Field& field,
					     void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::EdgeBlock::internal_put_field_data(const Ioss::Field& field,
					     void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::EdgeBlock::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "distribution_factor_count") {
    if (field_exists("distribution_factors")) {
      int nnodes = topology()->number_nodes();
      int nedge  = get_property("entity_count").get_int();
      return Ioss::Property(my_name, nnodes*nedge);
    } else {
      return Ioss::Property(my_name, 0);
    }
  }
  else
    return Ioss::EntityBlock::get_implicit_property(my_name);
}

void Ioss::EdgeBlock::block_membership(std::vector<std::string> &block_members)
{
  // Simplest case.  If the surfaces are split by element block, then this will return non-null
  // and we are done.
  const Ioss::ElementBlock *eb = parent_element_block();
  if (eb != NULL) {
    block_members.push_back(eb->name());
    return;
  }

  if (blockMembership.empty()) {
    get_database()->compute_block_membership(this, blockMembership);
  } 
  block_members = blockMembership;
}
