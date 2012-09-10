// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Ioss_DatabaseIO.h>
#include <Ioss_GroupingEntity.h>
#include <Ioss_Property.h>
#include <Ioss_Region.h>
#include <Ioss_Utils.h>
#include <Ioss_VariableType.h>
#include <assert.h>
#include <stddef.h>
#include <iostream>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_EntityType.h"
#include "Ioss_Field.h"
#include "Ioss_FieldManager.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_State.h"

Ioss::GroupingEntity::GroupingEntity()
  : entityCount(0), entityName("invalid"), database_(NULL), entityState(STATE_CLOSED),
    attributeCount(0)
{}

Ioss::GroupingEntity::GroupingEntity(Ioss::DatabaseIO *io_database,
				     const std::string& my_name,
				     int64_t entity_count)
  : entityCount(entity_count), entityName(my_name), database_(io_database),
    entityState(STATE_CLOSED), attributeCount(0)    
{
  properties.add(Ioss::Property("name", my_name));

  properties.add(Ioss::Property("entity_count", entity_count));

  properties.add(Ioss::Property(this, "attribute_count",
				Ioss::Property::INTEGER));

  if (my_name != "null_entity") {
    Ioss::Field::BasicType int_type = Ioss::Field::INTEGER;
    if (io_database != NULL)
      int_type = field_int_type();
    fields.add(Ioss::Field("ids", int_type, "scalar",
			   Ioss::Field::MESH, entity_count));
  }
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

Ioss::DatabaseIO* Ioss::GroupingEntity::get_database() const
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

void Ioss::GroupingEntity::set_database(Ioss::DatabaseIO *io_database)
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
Ioss::GroupingEntity::get_implicit_property(const std::string& my_name) const
{
  // Handle properties generic to all GroupingEntities.
  // These include:
  if (my_name == "attribute_count") {
    count_attributes();
    return Ioss::Property(my_name, static_cast<int>(attributeCount));
  }

  // End of the line. No property of this name exists.
  std::ostringstream errmsg;
  errmsg << "\nERROR: Property '" << my_name << "' does not exist on "
         << type_string() << " " << name() << "\n\n";
  IOSS_ERROR(errmsg);
  return Ioss::Property();
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
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for input on "
	   << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

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
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for output on "
	   << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.transform(data);
  return internal_put_field_data(field, data, data_size);
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<double>    &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for input on "
	   << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::REAL);
  
  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(double);
  int retval = internal_get_field_data(field, TOPTR(data), data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(TOPTR(data));

  return retval;
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<int>     &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for input on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::INTEGER);
  
  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(int);
  int retval = internal_get_field_data(field, TOPTR(data), data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(TOPTR(data));

  return retval;
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<int64_t> &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for input on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::INT64);
  
  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(int64_t);
  int retval = internal_get_field_data(field, TOPTR(data), data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(TOPTR(data));

  return retval;
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<char>     &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for input on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::CHARACTER);
  
  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(char);
  int retval = internal_get_field_data(field, TOPTR(data), data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(TOPTR(data));

  return retval;
}

int Ioss::GroupingEntity::get_field_data(const std::string& field_name,
					 std::vector<Complex> &data) const

{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for input on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::COMPLEX);

  data.resize(field.raw_count() * field.raw_storage()->component_count());
  size_t data_size = data.size() * sizeof(Complex);
  int retval = internal_get_field_data(field, TOPTR(data), data_size);

  // At this point, transform the field if specified...
  if (retval >= 0)
    field.transform(TOPTR(data));

  return retval;
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<double> &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for output on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::REAL);
  size_t data_size = data.size() * sizeof(double);
  field.transform(TOPTR(data));
  return internal_put_field_data(field, TOPTR(data), data_size);
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<int> &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for output on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::INTEGER);
  size_t data_size = data.size() * sizeof(int);
  field.transform(TOPTR(data));
  return internal_put_field_data(field, TOPTR(data), data_size);
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<int64_t> &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for output on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::INT64);
  size_t data_size = data.size() * sizeof(int64_t);
  field.transform(TOPTR(data));
  return internal_put_field_data(field, TOPTR(data), data_size);
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<char> &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for output on "
	      << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::CHARACTER);
  size_t data_size = data.size() * sizeof(char);
  field.transform(TOPTR(data));
  return internal_put_field_data(field, TOPTR(data), data_size);
}

int Ioss::GroupingEntity::put_field_data(const std::string& field_name,
					 std::vector<Complex> &data) const
{
  if (!field_exists(field_name)) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Field '" << field_name << "' does not exist for output on "
	   << type_string() << " " << name() << "\n\n";
    IOSS_ERROR(errmsg);
  }

  Ioss::Field field = get_field(field_name);
  field.check_type(Ioss::Field::COMPLEX);
  size_t data_size = data.size() * sizeof(Complex);
  field.transform(TOPTR(data));
  return internal_put_field_data(field, TOPTR(data), data_size);
}

void Ioss::GroupingEntity::count_attributes() const
{
  if (attributeCount > 0)
    return;
  else {
    // If the set has a field named "attribute", then the number of
    // attributes is equal to the component count of that field...
    NameList results_fields;
    field_describe(Ioss::Field::ATTRIBUTE, &results_fields);

    Ioss::NameList::const_iterator IF;
    for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
      std::string field_name = *IF;
      if (field_name != "attribute" ||
	  (field_name == "attribute" && results_fields.size() == 1)) {
	Ioss::Field field = get_field(field_name);
	attributeCount += field.raw_storage()->component_count();
      }
    }
  }
}

