/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_PropertyManager_h
#define SIERRA_Ioss_PropertyManager_h

#include <Ioss_CodeTypes.h>
#include <map>
#include <vector>

#include <Ioss_Property.h>
#include <string>

namespace Ioss {
  typedef std::vector<std::string> NameList;
  typedef std::map<std::string, Property, std::less<std::string> > PropMapType;
  typedef PropMapType::value_type ValuePair;

  class PropertyManager {
  public:
    PropertyManager();
    ~PropertyManager();

    // Add the specified property to the list.
    void add(const Property& new_prop);

    // Assumes: Property 'name' must exist.
    void erase(const std::string& property_name);

    // Checks if a property with 'property_name' exists in the database.
    bool exists(const std::string& property_name) const;

    Property get(const std::string& property_name) const;

    // Returns the names of all properties
    int describe(NameList* names) const;

    size_t count() const;

  private:
    // Disallow copying; don't implement...
    PropertyManager(const PropertyManager&);
    PropertyManager& operator=(const PropertyManager& from); // do not implement

    PropMapType properties;
  };
}
#endif
