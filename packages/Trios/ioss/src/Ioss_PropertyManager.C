/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_PropertyManager.h>

#include <assert.h>

#include <Ioss_Property.h>
#include <Ioss_Utils.h>
#include <string>

Ioss::PropertyManager::PropertyManager() {}

Ioss::PropertyManager::~PropertyManager()
{
  try {
    properties.clear();
  } catch (...) {
  }
}

void Ioss::PropertyManager::add(const Ioss::Property& new_prop)
{
  PropMapType::iterator iter = properties.find(new_prop.get_name());
  if (iter != properties.end()) {
    properties.erase(iter);
  }
  properties.insert(ValuePair(new_prop.get_name(), new_prop));
}

// Checks if a property with 'property_name' exists in the database.
bool Ioss::PropertyManager::exists(const std::string& property_name) const
{
  return (properties.find(property_name) != properties.end());
}

Ioss::Property Ioss::PropertyManager::get(const std::string& property_name) const
{
  PropMapType::const_iterator iter = properties.find(property_name);
  if (iter == properties.end()) {
    std::ostringstream errmsg;
    errmsg << "ERROR: Could not find property '" << property_name << "'\n";
    IOSS_ERROR(errmsg);
  }
  return (*iter).second;
}

void Ioss::PropertyManager::erase(const std::string& property_name)
{
  PropMapType::iterator iter = properties.find(property_name);
  if (iter != properties.end()) {
    properties.erase(iter);
  }
}

// Returns the names of all properties
int Ioss::PropertyManager::describe(NameList *names) const
{
  int the_count = 0;
  PropMapType::const_iterator I;
  for (I = properties.begin(); I != properties.end(); ++I) {
    names->push_back((*I).first);
    the_count++;
  }
  return the_count;
}


size_t Ioss::PropertyManager::count() const
{
  return properties.size();
}
