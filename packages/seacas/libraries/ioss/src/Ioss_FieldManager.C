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

const Ioss::Field &Ioss::FieldManager::getref(const std::string& field_name) const
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

