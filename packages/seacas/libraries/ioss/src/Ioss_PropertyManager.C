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
