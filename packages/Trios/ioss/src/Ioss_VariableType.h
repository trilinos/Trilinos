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

#ifndef IOSS_Ioss_VariableType_h
#define IOSS_Ioss_VariableType_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <cstring>
#include <vector>
#include <map>

namespace Ioss {
  typedef std::vector<std::string> NameList;

  class VariableType;
  
  typedef std::map<std::string, VariableType*, std::less<std::string> > VariableTypeMap;
  typedef VariableTypeMap::value_type VTM_ValuePair;

  class Registry {
  public:
    void insert(const Ioss::VTM_ValuePair &value, bool delete_me);
    VariableTypeMap::iterator begin() {return m_registry.begin();}
    VariableTypeMap::iterator end()   {return m_registry.end();}
    VariableTypeMap::iterator find(const std::string &type) {return m_registry.find(type);}

    ~Registry();
    std::map<std::string, std::string> customFieldTypes;
  private:
    Ioss::VariableTypeMap m_registry;
    std::vector<Ioss::VariableType*> m_deleteThese;
  };

#define MAX_SUFFIX 8
  struct Suffix {
    Suffix(const char new_data[MAX_SUFFIX]) {std::strncpy(data, new_data,         MAX_SUFFIX); data[MAX_SUFFIX]='\0';}
    Suffix(const std::string &new_data)     {std::strncpy(data, new_data.c_str(), MAX_SUFFIX); data[MAX_SUFFIX]='\0';}
    bool operator==(const std::string &str) const
    {return std::strncmp(data, str.c_str(), MAX_SUFFIX) == 0;}
    bool operator!=(const std::string &str) const
    {return std::strncmp(data, str.c_str(), MAX_SUFFIX) != 0;}
    char data[MAX_SUFFIX+1];
  };

  class VariableType {
  public:

    static void alias(const std::string& base, const std::string& syn);
    static int describe(NameList *names);
    static bool create_named_suffix_field_type(const std::string &type_name, std::vector<std::string> &suffices);
    static bool get_field_type_mapping(const std::string &field, std::string *type);
    static bool add_field_type_mapping(const std::string &field, const std::string &type);

    virtual ~VariableType();
    int component_count() const;

    // Override this function if the derived class has no suffices
    // For example, a 'vector_2d' has suffices "x" and "y"
    // A 'quad4' has no suffices...
    virtual int suffix_count() const;
    std::string name() const;

    static std::string numeric_label(int which, int ncomp,
					      const std::string &name);
    virtual std::string label(int which, const char suffix_sep='_') const = 0;
    virtual std::string label_name(const std::string& base, int which,
				      const char suffix_sep='_') const;
    virtual bool match(const std::vector<Suffix> &suffices) const;

    static const VariableType* factory(const std::string& name, int copies = 1);
    static const VariableType* factory(const std::vector<Suffix> &suffices);

  protected:
    VariableType(const std::string& type, int comp_count, bool delete_me = false);
    static Registry& registry();

  private:

    const std::string name_;
    int componentCount;

    VariableType(const VariableType&); // Do not implement...
    VariableType& operator=(const VariableType&); // Do not implement...

    static bool build_variable_type(const std::string& type);
  };
}
inline std::string Ioss::VariableType::name() const
{return name_;}

inline int Ioss::VariableType::component_count() const
{return componentCount;}

inline int Ioss::VariableType::suffix_count() const
{return componentCount;}
#endif
