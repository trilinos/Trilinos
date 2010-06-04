/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_ElementBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <Ioss_ElementTopology.h>

#include <string>

Ioss::ElementBlock::ElementBlock(const Ioss::DatabaseIO *io_database,
				 const std::string& my_name,
				 const std::string& element_type,
				 size_t number_elements,
				 size_t number_attributes)
  : Ioss::EntityBlock(io_database, my_name, element_type,
		      std::string("unknown"), number_elements),
    idOffset(0), elementCount(number_elements), attributeCount(number_attributes)
{
  properties.add(Ioss::Property(this, "attribute_count",
				Ioss::Property::INTEGER));
}

Ioss::ElementBlock::~ElementBlock() {}

void Ioss::ElementBlock::count_attributes() const
{
  if (attributeCount > 0)
    return;
  else {
    // If the block has a field named "attribute", then the number of
    // attributes is equal to the component count of that field...
    if (field_exists("attribute")) {
      Ioss::Field field = get_field("attribute");
      attributeCount = field.raw_storage()->component_count();
      return;
    } else {
      NameList results_fields;
      field_describe(Ioss::Field::ATTRIBUTE, &results_fields);

      Ioss::NameList::const_iterator IF;
      for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
	std::string field_name = *IF;
	Ioss::Field field = get_field(field_name);
	attributeCount += field.raw_storage()->component_count();
      }
    }
  }
}

Ioss::Property Ioss::ElementBlock::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "attribute_count") {
    count_attributes();
    return Ioss::Property(my_name, static_cast<int>(attributeCount));
  } else {
    return Ioss::EntityBlock::get_implicit_property(my_name);
  }
}

int Ioss::ElementBlock::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::ElementBlock::internal_put_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

void Ioss::ElementBlock::get_block_adjacencies(std::vector<std::string> &block_adjacency) const
{
  get_database()->get_block_adjacencies(this, block_adjacency);
}
