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
    const Field &getref(const std::string& field_name) const;

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
