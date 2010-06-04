/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_NodeSet.h>

#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_Utils.h>
#include <string>

namespace {
  const std::string SCALAR() { return std::string("scalar");}
}

Ioss::NodeSet::NodeSet(const Ioss::DatabaseIO *io_database, const std::string& my_name,
		       size_t number_nodes)
  : Ioss::GroupingEntity(io_database, my_name)
{
  properties.add(Ioss::Property("entity_count", static_cast<int>(number_nodes)));
  properties.add(Ioss::Property("distribution_factor_count",
				static_cast<int>(number_nodes)));
      // Add the standard fields...
  fields.add(Ioss::Field("ids", Ioss::Field::INTEGER, SCALAR(),
			 Ioss::Field::MESH, number_nodes));

  fields.add(Ioss::Field("distribution_factors",
			 Ioss::Field::REAL, SCALAR(),
			 Ioss::Field::MESH, number_nodes));

}

int Ioss::NodeSet::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::NodeSet::internal_put_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property Ioss::NodeSet::get_implicit_property(const std::string& my_name) const
{
  return Ioss::GroupingEntity::get_implicit_property(my_name);
}

void Ioss::NodeSet::block_membership(std::vector<std::string> &block_members)
{
  block_members.push_back("nodeblock_1");
}
