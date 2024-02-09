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
//

#ifndef stk_expreval_variable_hpp
#define stk_expreval_variable_hpp

#include "stk_util/util/string_case_compare.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <iostream>
#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>
#include <sstream>
#include <memory>


namespace stk {
namespace expreval {

class NgpNode;

class Variable
{
public:
  enum Type        {DOUBLE, INTEGER};
  enum Use         {DEPENDENT, INDEPENDENT};
  enum ArrayOffset {ZERO_BASED_INDEX, ONE_BASED_INDEX};

  explicit Variable(const std::string & name)
    : m_type(DOUBLE),
      m_use(INDEPENDENT),
      m_size(1),
      m_doublePtr(&m_doubleValue),
      m_doubleValue(0.0),
      m_name(name),
      m_index(0),
      m_isModifiable(true)
  {
  }

  Variable(Type type, const std::string & name)
    : m_type(type),
      m_use(INDEPENDENT),
      m_size(1),
      m_name(name),
      m_index(0),
      m_isModifiable(true)
  {
    switch (type) {
    case DOUBLE:
      m_doublePtr = &m_doubleValue;
      m_doubleValue = 0.0;
      break;
    case INTEGER:
      m_intPtr = &m_intValue;
      m_intValue = 0;
      break;
    }
  }

  Variable(const double & address, const std::string & name, unsigned definedLength=std::numeric_limits<int>::max())
    : m_type(DOUBLE),
      m_use(INDEPENDENT),
      m_size(definedLength),
      m_doublePtr(&address),
      m_doubleValue(0.0),
      m_name(name),
      m_index(0),
      m_isModifiable(false)
  {
  }

  Variable(double & address, const std::string & name, unsigned definedLength=std::numeric_limits<int>::max())
    : m_type(DOUBLE),
      m_use(INDEPENDENT),
      m_size(definedLength),
      m_doublePtr(&address),
      m_doubleValue(0.0),
      m_name(name),
      m_index(0),
      m_isModifiable(true)
  {
  }

  Variable(const int & address, const std::string & name, unsigned definedLength=std::numeric_limits<int>::max())
    : m_type(INTEGER),
      m_use(INDEPENDENT),
      m_size(definedLength),
      m_intPtr(&address),
      m_intValue(0),
      m_name(name),
      m_index(0),
      m_isModifiable(false)
  {
  }

  Variable(int & address, const std::string & name, unsigned definedLength=std::numeric_limits<int>::max())
    : m_type(INTEGER),
      m_use(INDEPENDENT),
      m_size(definedLength),
      m_intPtr(&address),
      m_intValue(0),
      m_name(name),
      m_index(0),
      m_isModifiable(true)
  {
  }

  Variable &operator=(double value) {
    STK_ThrowRequireMsg(m_size == 1 || m_size == std::numeric_limits<int>::max(),
                        "Invalid use of equal on the multi-component variable " + m_name);
    STK_ThrowRequireMsg(m_isModifiable, "Variable " + m_name + " is const and cannot be modified.");

    if (m_type == INTEGER) {
      *const_cast<int*>(m_intPtr) = static_cast<int>(value);
    }
    else if (m_type == DOUBLE) {
      *const_cast<double*>(m_doublePtr) = value;
    }

    return *this;
  }

private:
  Variable(const Variable &);
  Variable &operator=(const Variable &);

public:

  void setDependent() { m_use = DEPENDENT; }

  bool isDependent() const { return m_use == DEPENDENT; }

  inline double getArrayValue(int index, ArrayOffset offsetType) const
  {
    if (m_type == DOUBLE) {
      STK_ThrowRequireMsg(m_doublePtr != nullptr, "Variable " + m_name + " does not have a valid value.");

      if (offsetType == ZERO_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 0, "Provided index for variable array " + m_name + " is less than 0.");
        STK_ThrowRequireMsg(index < m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        return m_doublePtr[index];
      }
      else if (offsetType == ONE_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 1, "Provided index for variable array " + m_name + " is less than 1.");
        STK_ThrowRequireMsg(index <= m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        return m_doublePtr[index-1];
      }
      else {
        STK_ThrowErrorMsg("Invalid ArrayOffsetType for variable " + m_name);
        return m_doublePtr[0];
      }
    }
    else {
      STK_ThrowRequireMsg(m_intPtr != nullptr, "Variable " + m_name + " does not have a valid value.");

      if (offsetType == ZERO_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 0, "Provided index for variable array " + m_name + " is less than 0.");
        STK_ThrowRequireMsg(index < m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        return static_cast<double>(m_intPtr[index]);
      }
      else if (offsetType == ONE_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 1, "Provided index for variable array " + m_name + " is less than 1.");
        STK_ThrowRequireMsg(index <= m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        return static_cast<double>(m_intPtr[index-1]);
      }
      else {
        STK_ThrowErrorMsg("Invalid ArrayOffsetType for variable " + m_name);
        return static_cast<double>(m_intPtr[0]);
      }
    }
  }

  inline void assignArrayValue(int index, ArrayOffset offsetType, double value) const
  {
    STK_ThrowRequireMsg(m_isModifiable, "Variable " + m_name + " is const and cannot be modified.");

    if (m_type == DOUBLE) {
      STK_ThrowRequireMsg(m_doublePtr != nullptr, "Variable " + m_name + " does not have a valid value.");

      if (offsetType == ZERO_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 0, "Provided index for variable array " + m_name + " is less than 0.");
        STK_ThrowRequireMsg(index < m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        const_cast<double &>(m_doublePtr[index]) = value;
      }
      else if (offsetType == ONE_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 1, "Provided index for variable array " + m_name + " is less than 1.");
        STK_ThrowRequireMsg(index <= m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        const_cast<double &>(m_doublePtr[index-1]) = value;
      }
      else {
        STK_ThrowErrorMsg("Invalid ArrayOffsetType for variable " + m_name);
      }
    }
    else {
      STK_ThrowRequireMsg(m_intPtr != nullptr, "Variable " + m_name + " does not have a valid value.");

      if (offsetType == ZERO_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 0, "Provided index for variable array " + m_name + " is less than 0.");
        STK_ThrowRequireMsg(index < m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        const_cast<int &>(m_intPtr[index]) = static_cast<int>(value);
      }
      else if (offsetType == ONE_BASED_INDEX) {
        STK_ThrowRequireMsg(index >= 1, "Provided index for variable array " + m_name + " is less than 1.");
        STK_ThrowRequireMsg(index <= m_size, "Provided index for variable array " + m_name + " exceeds array upper bound.");

        const_cast<int &>(m_intPtr[index-1]) = static_cast<int>(value);
      }
      else {
        STK_ThrowErrorMsg("Invalid ArrayOffsetType for variable " + m_name);
      }
    }
  }

  inline Variable & bind(const double & value_ref, int definedLength=std::numeric_limits<int>::max()) {
    m_type = DOUBLE;
    m_doublePtr = &value_ref;
    m_size = definedLength;
    m_isModifiable = false;
    return *this;
  }

  inline Variable & bind(double & value_ref, int definedLength=std::numeric_limits<int>::max()) {
    m_type = DOUBLE;
    m_doublePtr = &value_ref;
    m_size = definedLength;
    m_isModifiable = true;
    return *this;
  }

  inline Variable & bind(const int & value_ref, int definedLength=std::numeric_limits<int>::max()) {
    m_type = INTEGER;
    m_intPtr = &value_ref;
    m_size = definedLength;
    m_isModifiable = false;
    return *this;
  }

  inline Variable & bind(int & value_ref, int definedLength=std::numeric_limits<int>::max()) {
    m_type = INTEGER;
    m_intPtr = &value_ref;
    m_size = definedLength;
    m_isModifiable = true;
    return *this;
  }

  inline Variable & unbind() {
    m_size = 1;
    switch (m_type) {
    case DOUBLE: {
      m_doublePtr = &m_doubleValue;
      m_doubleValue = 0.0;
      m_isModifiable = true;
      break;
    }
    case INTEGER: {
      m_intPtr = &m_intValue;
      m_intValue = 0;
      m_isModifiable = true;
      break;
    }
    }
    return *this;
  }

  inline Variable & deactivate() {
    m_size = 1;
    switch (m_type) {
    case DOUBLE: {
      m_doublePtr = nullptr;
      m_doubleValue = 0.0;
      m_isModifiable = false;
      break;
    }
    case INTEGER: {
      m_intPtr = nullptr;
      m_intValue = 0;
      m_isModifiable = false;
      break;
    }
    }
    return *this;
  }

  const double * getAddress() const { return m_doublePtr; }

  KOKKOS_FUNCTION
  int getLength() const { return m_size; }

  inline double getValue() const {
    STK_ThrowRequireMsg(m_size == 1 || m_size == std::numeric_limits<int>::max(),
     "Invalid direct access of array variable " + m_name + ", must access by index");

    STK_ThrowRequireMsg(((m_type == DOUBLE) && (m_doublePtr != nullptr)) || 
                        ((m_type == INTEGER) && (m_intPtr != nullptr)), 
                        "Variable " + m_name + " does not have a valid value."); 

    switch (m_type) {
    case DOUBLE:
      return *m_doublePtr;
    case INTEGER:
      return *m_intPtr;
    }

    STK_ThrowErrorMsg("Invalid variable type for variable " + m_name);
    return *m_doublePtr;
  }

  KOKKOS_FUNCTION
  int get_index() { return m_index; }

  void set_index(int index) { m_index = index; }

  KOKKOS_FUNCTION
  Type get_type() { return m_type; }

private:
  Type m_type;
  Use m_use;
  int m_size;

  union {
    const double * m_doublePtr;
    const int * m_intPtr;
  };

  union {
    double m_doubleValue;
    int m_intValue;
  };

  const std::string m_name;
  int m_index;
  bool m_isModifiable;
};


class VariableMap : public std::map<std::string, std::shared_ptr<Variable>, LessCase>
{
public:
  typedef std::map<std::string, std::shared_ptr<Variable>, LessCase >::value_type value_type;

  class Resolver
  {
  public:
    Resolver()
    {}

    virtual ~Resolver()
    {}

    virtual void resolve(VariableMap::iterator &)
    {
    }
  };

private:
  class DefaultResolver : public Resolver
  {
  public:
    DefaultResolver()
    {}

    virtual ~DefaultResolver()
    {}
  };

public:
  static Resolver &getDefaultResolver();

  VariableMap(Resolver & resolver = getDefaultResolver())
    : std::map<std::string, std::shared_ptr<Variable>, LessCase>(),
      m_resolver(resolver)
  {}

  virtual ~VariableMap() {}

  Variable* operator[](const std::string & s) {
    std::pair<iterator,bool> i = insert(std::pair<const std::string,
                                        std::shared_ptr<Variable>>(s, std::shared_ptr<Variable>(nullptr)));
    if (i.second) {
      (*i.first).second = std::make_shared<Variable>(s);
    }
    return (*i.first).second.get();
  }

  Resolver & getResolver() { return m_resolver; }

private:
  Resolver & m_resolver;
};

} // namespace expreval
} // namespace stk

#endif // stk_expreval_variable_hpp
