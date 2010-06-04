/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_CommSet.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>

#include <string>
#include <assert.h>

Ioss::CommSet::CommSet(const Ioss::DatabaseIO *io_database,
		       const std::string& my_name,
		       const std::string& entity_type,
		       size_t entity_count)
  : Ioss::GroupingEntity(io_database, my_name)
{
  assert(entity_type == "node" || entity_type == "face");
  properties.add(Ioss::Property("entity_count", static_cast<int>(entity_count)));
  properties.add(Ioss::Property("entity_type",  entity_type));

  // Field contains a pair of type [entity_id, shared_cpu]
  fields.add(Ioss::Field("entity_processor", Ioss::Field::INTEGER, "pair",
			 Ioss::Field::COMMUNICATION, entity_count));
}

int Ioss::CommSet::internal_get_field_data(const Ioss::Field& field,
				 void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::CommSet::internal_put_field_data(const Ioss::Field& field,
				 void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::CommSet::get_implicit_property(const std::string& my_name) const
{
  return Ioss::GroupingEntity::get_implicit_property(my_name);
}
