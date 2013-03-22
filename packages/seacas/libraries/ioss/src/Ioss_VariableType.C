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

#include <Ioss_CompositeVariableType.h>
#include <Ioss_ConstructedVariableType.h>
#include <Ioss_NamedSuffixVariableType.h>
#include <Ioss_Utils.h>
#include <Ioss_VariableType.h>
#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
  class Deleter {
  public:
    void operator()(Ioss::VariableType* t) {delete t;}
  };
}

void Ioss::Registry::insert(const Ioss::VTM_ValuePair &value, bool delete_me)
{
  m_registry.insert(value);
  if (delete_me) {
    m_deleteThese.push_back(value.second);
  }
}

Ioss::Registry::~Registry()
{
  if (!m_deleteThese.empty()) 
    std::for_each(m_deleteThese.begin(), m_deleteThese.end(), Deleter());
}

Ioss::VariableType::~VariableType() {}

Ioss::VariableType::VariableType(const std::string& type, int comp_count, bool delete_me)
  : name_(type), componentCount(comp_count)
{
  std::string low_type = Ioss::Utils::lowercase(type);
  registry().insert(Ioss::VTM_ValuePair(low_type, this), delete_me);

  // Register uppercase version also
  std::string up_type = Ioss::Utils::uppercase(type);
  registry().insert(Ioss::VTM_ValuePair(up_type, this), false);

}

void Ioss::VariableType::alias(const std::string& base, const std::string& syn)
{
  registry().insert(Ioss::VTM_ValuePair(Ioss::Utils::lowercase(syn),
					(Ioss::VariableType*)factory(base)), false);
  // Register uppercase version also
  std::string up_type = Ioss::Utils::uppercase(syn);
  registry().insert(Ioss::VTM_ValuePair(up_type,
					(Ioss::VariableType*)factory(base)), false);

}

Ioss::Registry& Ioss::VariableType::registry()
{
  static Registry registry_;
  return registry_;
}

int Ioss::VariableType::describe(NameList *names)
{
  int count = 0;
  Ioss::VariableTypeMap::const_iterator I;
  for (I = registry().begin(); I != registry().end(); ++I) {
    names->push_back((*I).first);
    count++;
  }
  return count;
}

bool Ioss::VariableType::add_field_type_mapping(const std::string &raw_field, const std::string &raw_type)
{
  // See if storage type 'type' exists...
  std::string field = Ioss::Utils::lowercase(raw_field);
  std::string type  = Ioss::Utils::lowercase(raw_type);
  if (registry().find(type) == registry().end())
    return false;

  // Add mapping.
  if (registry().customFieldTypes.insert(std::make_pair(field, type)).second)
    return true;
  else
    return false;
}

bool Ioss::VariableType::create_named_suffix_field_type(const std::string &type_name, std::vector<std::string> &suffices)
{
  size_t count = suffices.size();
  if (count < 1)
    return false;

  std::string low_name = Ioss::Utils::lowercase(type_name);
  // See if the variable already exists...
  if (registry().find(low_name) != registry().end())
    return false;
  
  // Create the variable.  Note that the 'true' argument means Ioss will delete the pointer.
  Ioss::NamedSuffixVariableType *var_type = new Ioss::NamedSuffixVariableType(low_name, count, true);

  for (size_t i=0; i < count; i++) {
    var_type->add_suffix(i+1, suffices[i]);
  }
  return true;
}

bool Ioss::VariableType::get_field_type_mapping(const std::string &field, std::string *type)
{
  // Returns true if a mapping exists, 'type' contains the mapped type.
  // Returns false if no custom mapping exists for this field.
  std::string low_field = Ioss::Utils::lowercase(field);
  
  if (registry().customFieldTypes.find(low_field) == registry().customFieldTypes.end()) {
    return false;
  }
  else {
    *type = registry().customFieldTypes.find(low_field)->second;
    return true;
  }
}

const Ioss::VariableType* Ioss::VariableType::factory(const std::string& raw_name, int copies)
{
  Ioss::VariableType* inst = NULL;
  std::string name = Ioss::Utils::lowercase(raw_name);
  Ioss::VariableTypeMap::iterator iter = registry().find(name);
  if (iter == registry().end()) {
    bool can_construct = build_variable_type(name);
    if (can_construct) {
      iter = registry().find(name);
      assert(iter != registry().end());
      inst = (*iter).second;
    } else {
      std::ostringstream errmsg;
      errmsg << "ERROR: The variable type '" << raw_name << "' is not supported.\n";
      IOSS_ERROR(errmsg);
    }
  } else {
    inst = (*iter).second;
  }

  if (copies != 1) {
    inst = CompositeVariableType::composite_variable_type(inst, copies);
  }
  assert(inst != NULL);
  return inst;
}

const Ioss::VariableType* Ioss::VariableType::factory(const std::vector<Ioss::Suffix> &suffices)
{
  size_t size = suffices.size();
  const Ioss::VariableType* ivt = NULL;
  if (size <= 1)
    return NULL; // All storage types must have at least 2 components.

  Ioss::VariableTypeMap::const_iterator I  = registry().begin();
  Ioss::VariableTypeMap::const_iterator IE = registry().end();

  bool match = false;
  while (I != IE) {
    ivt = (*I++).second;
    if ( ivt->suffix_count() == (int)size) {
      if (ivt->match(suffices)) {
	match = true;
	break;
      }
    }
  }
  if (match == false) {
    match = true;
    // Check if the suffices form a sequence (1,2,3,...,N)
    // This indicates a "component" variable type that is
    // constructed "on-the-fly" for use in Sierra
    //
    // Maximum suffix size is 3.  If changed, then code in
    // src/framewk/Frio_RestartUtils::variable_name_kluge must also
    // be changed.

    char digits[6]; // Include trailing null
    assert(size < 100000);

    // Create a format for our use...
    char format[5];
    if (size < 10)
      std::strcpy(format, "%01d");
    else if (size < 100)
      std::strcpy(format, "%02d");
    else if (size < 1000)
      std::strcpy(format, "%03d");
    else if (size < 10000)
      std::strcpy(format, "%04d");
    else
      std::strcpy(format, "%05d");

    for (size_t i=0; i < size; i++) {
      std::sprintf(digits, format, i+1);
      if (std::strcmp(&suffices[i].data[0], &digits[0]) != 0) {
	match = false;
	break;
      }
    }
    if (match) {
      // Create a new type.  For now, assume that the base type is
      // "Real" (Sierra type).  The name of the new type is
      // "Real[component_count]"
      // Note that this type has not yet been constructed since
      // it would have been found above.
      ivt = new Ioss::ConstructedVariableType(size,true);
    } else {
      ivt = NULL;
    }
  }
  return ivt;
}

bool Ioss::VariableType::match(const std::vector<Ioss::Suffix> &suffices) const
{
  bool result = true;
  if ((int)suffices.size() == suffix_count()) {
    for (int i=0; i < suffix_count(); i++) {
      if (suffices[i] != label(i+1)) {
	result = false;
	break;
      }
    }
  } else {
    result = false;
  }
  return result;
}

std::string Ioss::VariableType::label_name(const std::string& base, int which,
				      const char suffix_sep) const
{
  static char tmp_sep[2];
  std::string my_name = base;
  std::string suffix = label(which, suffix_sep);
  if (!suffix.empty()) {
    if (suffix_sep != 0) {
      tmp_sep[0] = suffix_sep;
      my_name += tmp_sep;
    }
    my_name += suffix;
  }
  return my_name;
}

bool Ioss::VariableType::build_variable_type(const std::string& raw_type)
{
  // See if this is a multi-component instance of a base type.
  // An example would be REAL[2] which is a basic real type with
  // two components.  The suffices would be .0 and .1

  std::string type = Ioss::Utils::lowercase(raw_type);
  
  // Step 0:
  // See if the type contains '[' and ']'
  char const *typestr = type.c_str();
  char const *lbrace =  std::strchr(typestr, '[');
  char const *rbrace = std::strrchr(typestr, ']');

  if (lbrace == NULL || rbrace == NULL) return false;

  // Step 1:
  // First, we split off the basename (REAL/INTEGER) from the component count ([2])
  // and see if the basename is a valid variable type and the count is a
  // valid integer.
  size_t len = type.length() + 1;
  char *typecopy = new char[len];
  std::strcpy(typecopy, typestr);

  char *base = std::strtok(typecopy, "[]");
  assert (base != NULL);
  Ioss::VariableTypeMap::iterator iter = Ioss::VariableType::registry().find(base);
  if (iter == registry().end()) {
    delete [] typecopy;
    return false;
  }

  char *countstr = std::strtok(NULL, "[]");
  assert (countstr != NULL);
  int count = std::atoi(countstr);
  if (count <= 0) {
    delete [] typecopy;
    return false;
  }

  // We now know we have a valid base type and an integer
  // specifying the number of 'components' in our new type.
  // Create the new type and register it in the registry...
  new Ioss::ConstructedVariableType(type, count, true);
  delete [] typecopy;
  return true;
}

std::string Ioss::VariableType::numeric_label(int which, int ncomp, const std::string &name)
{
  char digits[8];
  // Create a format for our use...
  char format[5];
  if (ncomp <          10)
    std::strcpy(format, "%01d");
  else if (ncomp <    100)
    std::strcpy(format, "%02d");
  else if (ncomp <   1000)
    std::strcpy(format, "%03d");
  else if (ncomp <  10000)
    std::strcpy(format, "%04d");
  else if (ncomp < 100000)
    std::strcpy(format, "%05d");
  else {
    std::ostringstream errmsg;
    errmsg << "ERROR: Variable '" << name << "' has " << ncomp
	   << " components which is larger than the current maximum"
	   << " of 100,000. Please contact developer.\n";
    IOSS_ERROR(errmsg);
  }

  std::sprintf(digits, format, which);
  return std::string(digits);
}
