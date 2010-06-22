/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_FieldManager_h
#define IOSS_Ioss_FieldManager_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Utils.h>
#include <map>
#include <vector>
#include <utility>

#include <Ioss_Field.h>
#include <string>

namespace Ioss {
  typedef std::vector<std::string> NameList;

  // This performs a case-insensitive string comparison which will make
  // the FieldMapType do case-insensitive comparisons between field names.
  class StringCmp : public std::binary_function<std::string,std::string,bool>
    {
    public:
      StringCmp() {}
      bool operator() (const std::string &s1, const std::string &s2) const
      {
	return Utils::case_strcmp(s1, s2) < 0;
      }
    };
  typedef std::map<std::string, Field, StringCmp > FieldMapType;
  typedef FieldMapType::value_type FieldValuePair;

  class FieldManager {
  public:
    FieldManager();
    ~FieldManager();

    // Assumes: Field with the same 'name' does not exist.
    // Add the specified field to the list.
    void add(const Field& new_field);

    // Assumes: Field 'name' must exist.
    void erase(const std::string& field_name);

    // Checks if a field with 'field_name' exists in the database.
    bool exists(const std::string& field_name) const;

    Field get(const std::string& field_name) const;

    // Returns the names of all fields
    int describe(NameList* names) const;

    // Returns the names of all fields with the specified 'RoleType'
    int describe(Field::RoleType role, NameList* names) const;

    size_t count() const;

  private:
    // Disallow copying; don't implement...
    FieldManager(const FieldManager&);
    FieldManager& operator=(const FieldManager&); // Do not implement

    FieldMapType fields;
  };
}
#endif
