/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_FieldManager.h>
#include <Ioss_Field.h>

#include <assert.h>
#include <string>

Ioss::FieldManager::FieldManager() {}

Ioss::FieldManager::~FieldManager() {}

void Ioss::FieldManager::add(const Ioss::Field& new_field)
{
  if (!exists(new_field.get_name())) {
    fields.insert(FieldValuePair(new_field.get_name(), new_field));
  }
}

// Checks if a field with 'field_name' exists in the database.
bool Ioss::FieldManager::exists(const std::string& field_name) const
{
  return (fields.find(field_name) != fields.end());
}

Ioss::Field Ioss::FieldManager::get(const std::string& field_name) const
{
  FieldMapType::const_iterator iter = fields.find(field_name);
  assert(iter != fields.end());
  return (*iter).second;
}

// Assumes: Field 'name' must exist.
void Ioss::FieldManager::erase(const std::string& field_name)
{
  assert(exists(field_name));
  FieldMapType::iterator iter = fields.find(field_name);
  if (iter != fields.end()) {
    fields.erase(iter);
  }
}

// Returns the names of all fields
int Ioss::FieldManager::describe(NameList *names) const
{
  int the_count = 0;
  FieldMapType::const_iterator I;
  for (I = fields.begin(); I != fields.end(); ++I) {
    names->push_back((*I).first);
    the_count++;
  }
  return the_count;
}


// Returns the names of all fields
int Ioss::FieldManager::describe(Ioss::Field::RoleType role, NameList *names) const
{
  int the_count = 0;
  FieldMapType::const_iterator I;
  for (I = fields.begin(); I != fields.end(); ++I) {
    if ((*I).second.get_role() == role) {
      names->push_back((*I).first);
      the_count++;
    }
  }
  return the_count;
}


size_t Ioss::FieldManager::count() const
{
  return fields.size();
}

