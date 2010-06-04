/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_NodeBlock.h>

#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_ElementTopology.h>
#include <assert.h>
#include <string>

namespace {
  const std::string VECTOR_2D() { return std::string("vector_2d");}
  const std::string VECTOR_3D() { return std::string("vector_3d");}
  const std::string UNKNOWN()   { return std::string("unknown");}
}

Ioss::NodeBlock::NodeBlock(const Ioss::DatabaseIO *io_database,
			   const std::string &my_name,
			   size_t node_count,
			   size_t degrees_of_freedom)
  : Ioss::EntityBlock(io_database, my_name, "node", UNKNOWN(), node_count)
{
  properties.add(Ioss::Property("component_degree",
				static_cast<int>(degrees_of_freedom)));

  std::string vector_name;
  assert(degrees_of_freedom == 2 || degrees_of_freedom == 3);
  if (degrees_of_freedom == 2)
    vector_name = VECTOR_2D();
  else if (degrees_of_freedom == 3)
    vector_name = VECTOR_3D();
  
  fields.add(Ioss::Field("mesh_model_coordinates",
			 Ioss::Field::REAL, vector_name,
			 Ioss::Field::MESH, node_count));
}

Ioss::NodeBlock::~NodeBlock() {}

Ioss::Property
Ioss::NodeBlock::get_implicit_property(const std::string& my_name) const
{
  return Ioss::EntityBlock::get_implicit_property(my_name);
}

int Ioss::NodeBlock::internal_get_field_data(const Ioss::Field& field,
				   void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::NodeBlock::internal_put_field_data(const Ioss::Field& field,
				   void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}
