/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_GroupingEntity.h>

#include <Ioss_Region.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Utils.h>
#include <Ioss_VariableType.h>
#include <string>

#include <iostream>
#include <assert.h>

Ioss::GroupingEntity::GroupingEntity(const Ioss::DatabaseIO *io_database,
					 const std::string& my_name)
  : entityName(my_name), database_(io_database), entityState(STATE_CLOSED)
{
  properties.add(Ioss::Property("name", my_name));
}

Ioss::GroupingEntity::~GroupingEntity()
{
  // Only deleted by owning entity (Ioss::Region)
  database_ = NULL;
}

// Default implementation is to do nothing. Redefined in Ioss::Region
// to actually delete the database.
void Ioss::GroupingEntity::delete_database()  {}

void Ioss::GroupingEntity::really_delete_database()
{
  Ioss::DatabaseIO *new_db = const_cast<Ioss::DatabaseIO*>(database_);
  delete new_db;
  new_db = NULL;
}

bool Ioss::GroupingEntity::is_alias(const std::string &my_name) const
{
  Region *region = database_->get_region();
  return region->get_alias(my_name) == entityName;
}

const Ioss::DatabaseIO* Ioss::GroupingEntity::get_database() const
{
  assert(database_ != NULL);
  return database_;
}

std::string Ioss::GroupingEntity::get_filename() const
{
  // Ok for database_ to be NULL at this point.
  if (database_ == NULL)
    return std::string();
  else
    return database_->get_filename();
}

void Ioss::GroupingEntity::set_database(const Ioss::DatabaseIO *io_database)
{
  assert(database_ == NULL);     // Must be unset if we are setting it.
  assert(io_database != NULL);  // Must be set to valid value
  database_ = io_database;
}

// Discuss Data Object functions:
// ---Affect the containing object:
//    open(in string object_name, out ?)
//    close()
//    destroy()
// ---Affect what the object contains:
//    set(in string propertyname, in any property_value)
//    get(in string propertyname, out any property_value)
//    add(in string propertyname);
//    delete(in string propertyname)
//    describe(out vector<Ioss::Properties>)
//
Ioss::State Ioss::GroupingEntity::get_state() const
{
  return entityState;
}

Ioss::Property
Ioss::GroupingEntity::get_implicit_property(const std::string& /* name */) const
{
  // Handle properties generic to all GroupingEntities.
  // These include:
  // --NONE--

  return Ioss::Property(); // Invalid property.
}

void Ioss::GroupingEntity::field_add(const Ioss::Field& new_field)
{
  size_t entity_size = get_property("entity_count").get_int();
  size_t field_size  = new_field.raw_count();
  if (entity_size != field_size && type() != REGION) {
    std::ostringstream errmsg;
    errmsg << "IO System error: The " << type_string() << " '"
	   << name() << "' has a size of "
	   << entity_size << ",\nbut the field '" << new_field.get_name()
	   << "' which is being output on that entity has a size of " << field_size
	   << ".\nThe sizes must match.  This is an internal error that should be reported.";
    IOSS_ERROR(errmsg);
  }
  fields.add(new_field);
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					void *data, size_t data_size) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  int retval = internal_get_field_data(field, data, data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(data);

  return retval;
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					void *data, size_t data_size) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.transform(data);
  return internal_put_field_data(field, data, data_size);
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<double>    &data) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::REAL);
  
  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(double);
  int retval = internal_get_field_data(field, &data[0], data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(&data[0]);

  return retval;
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<int>     &data) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::INTEGER);
  
  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(int);
  int retval = internal_get_field_data(field, &data[0], data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(&data[0]);

  return retval;
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<Complex> &data) const

{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::COMPLEX);

  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(Complex);
  int retval = internal_get_field_data(field, &data[0], data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(&data[0]);

  return retval;
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<double> &data) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::REAL);
  size_t data_size = data.size() * sizeof(double);
  field.transform(&data[0]);
  return internal_put_field_data(field, &data[0], data_size);
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<int> &data) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::INTEGER);
  size_t data_size = data.size() * sizeof(int);
  field.transform(&data[0]);
  return internal_put_field_data(field, &data[0], data_size);
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<Complex> &data) const
{
  assert(field_exists(field_name));
  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::COMPLEX);
  size_t data_size = data.size() * sizeof(Complex);
  field.transform(&data[0]);
  return internal_put_field_data(field, &data[0], data_size);
}

