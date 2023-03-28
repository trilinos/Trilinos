// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_UTIL_COMMANDLINE_PARSEDOPTIONS_HPP
#define STK_UTIL_COMMANDLINE_PARSEDOPTIONS_HPP

#include "stk_util/stk_config.h"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/diag/StringUtil.hpp"
#include <sstream>
#include <string>
#include <map>
#include <vector>

namespace stk {

class VariableType
{
public:
  VariableType()
    : val(),
      isParsed(false)
  {}
  VariableType(const std::string& value)
    : val(value),
      isParsed(false)
  {}
  ~VariableType() {}

  operator const std::string&() const { return val; }
  operator std::string&() { return val; }

  bool empty() const { return val.empty(); }

  template<typename T>
  const T& as() const
  {
    STK_ThrowRequireMsg(!val.empty(), "Error in VariableType::as, internal value is empty.");
    STK_ThrowRequireMsg(sizeof(T) <= maxSize, std::string("Error: ") + sierra::demangle(typeid(T).name())
                    + " size too big in in VariableType::as()");
    std::istringstream iss(val);
    T* Tptr = new (asData) T(0);
    iss >> *Tptr;
    STK_ThrowRequireMsg(!iss.fail(), "Error in VariableType::as, failed to convert '" << val <<
                    "' to type " << sierra::demangle(typeid(T).name()) );
    return *Tptr;
  }

  const std::string& value() const { return val; }
  std::string& value() { return val; }

  void set_parsed() { isParsed = true; }
  bool is_parsed() const { return isParsed; }

private:
  static constexpr unsigned maxSize = 64;
  std::string val;
  mutable char asData[maxSize];
  bool isParsed;
};

template<>
inline
const std::string& VariableType::as<std::string>() const
{
  return val;
}

class ParsedOptions
{
public:
    ParsedOptions()
    : keyToIndex(), indexToValues() {}

    ~ParsedOptions() {}

    bool empty() const { return keyToIndex.empty(); }

    bool insert(const std::string& key, const std::string& val)
    {
      KeyMapType::iterator it = keyToIndex.find(key);
      
      if (it == keyToIndex.end()) {
        size_t idx = indexToValues.size();
        keyToIndex.insert(std::make_pair(key,idx));
        indexToValues.push_back(VariableType(val));
        return true;
      }
      else {
        indexToValues[it->second].value() = val;
      }

      return false;
    }

    bool insert(const std::string& name, const std::string& abbrev, const std::string& val)
    {
      KeyMapType::iterator it = keyToIndex.find(name);
      
      if (it == keyToIndex.end()) {
        size_t idx = indexToValues.size();
        keyToIndex.insert(std::make_pair(name,idx));
        if (!abbrev.empty()) {
          keyToIndex.insert(std::make_pair(abbrev,idx));
        }
        indexToValues.push_back(VariableType(val));
        return true;
      }
      else {
        indexToValues[it->second].value() = val;
      }

      return false;
    }

    void set_parsed(const std::string& key)
    {
      KeyMapType::iterator it = keyToIndex.find(key);

      if (it != keyToIndex.end()) {
        indexToValues[it->second].set_parsed();
      }
    }

    bool is_parsed(const std::string& key) const
    {
      KeyMapType::const_iterator it = keyToIndex.find(key);

      if (it != keyToIndex.end()) {
        return indexToValues[it->second].is_parsed();
      }
      else {
        return false;
      }
    }

    size_t count(const std::string& key) const { return keyToIndex.count(key); }

    const VariableType& operator[](const std::string& key) const
    {
      static VariableType s_var("");
      KeyMapType::const_iterator it = keyToIndex.find(key);
      if (it != keyToIndex.end()) {
        return indexToValues[it->second];
      }
      return s_var;
    }

private:
    using KeyMapType = std::map<std::string,unsigned>;

    std::map<std::string,unsigned> keyToIndex;
    std::vector<VariableType> indexToValues;
};

}

#endif //STK_UTIL_COMMANDLINE_PARSEDOPTIONS_HPP
