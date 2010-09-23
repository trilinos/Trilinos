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

#ifndef IOSS_Ioss_PropertyManager_h
#define IOSS_Ioss_PropertyManager_h

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
