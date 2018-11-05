// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 

#ifndef stk_expreval_Variable_hpp
#define stk_expreval_Variable_hpp

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


#include <stk_util/util/string_case_compare.hpp>

namespace stk {
namespace expreval {

/**
 * @brief Class <b>Variable</b> defines a variable for the expression evaluation
 * virtual machine.  Typing is limited to int and double.  While all operations are
 * performed as doubles, the value and results may be retrieved and stored as integer
 * values.
 *
 * Array indexing is also supported, but no bounds checkin is available and only double
 * type array may be addressed.
 *
 * Variable are implemented as a pointer, "typed" by the <b>Type</b> enumeration
 * to the actual value.  At construction, this may point to a "local" double value or
 * integer value, or the address may be specified.  In either case, a type and address
 * for the variable may be provided via variable "binding".
 *
 */
class Variable
{
public:
  /**
   * @brief Enumeration <b>Type</b> lists the variable data types.  <b>double</b> and
   * <b>int</b> are currently supported.
   */
  enum Type        {DOUBLE, INTEGER};
  enum Use         {DEPENDENT, INDEPENDENT};
  enum ArrayOffset {ZERO_BASED_INDEX, ONE_BASED_INDEX};
  
  /**
   * Creates a new <b>Variable</b> instance.
   */
  Variable(const std::string& name)
    : m_type(DOUBLE),
      m_use(INDEPENDENT),
      m_size(1),
      m_doublePtr(&m_doubleValue),
      m_doubleValue(0.0),
      m_name(name)
  {
  }

  /**
   * Creates a new <b>Variable</b> instance.  The new variable will be local and
   * have the type specified by <b>type</b>.
   *
   * @param type		an <b>Type</b> value of the type of the new
   *				variable.
   */
  explicit Variable(Type type, const std::string& name)
    : m_type(type),
      m_use(INDEPENDENT),
      m_size(1),
      m_name(name)
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

  /**
   * Creates a new <b>Variable</b> instance.  The new variable will use the
   * address specified by <b>address</b> and be of type double.  This address must
   * remain in scope during the lifetime of the variable are until rebound to a new
   * address.
   *
   * @param address		an <b>double</b> reference to the value for this
   *				variable.
   *
   */
  explicit Variable(double &address, const std::string& name, unsigned definedLength=std::numeric_limits<int>::max())
    : m_type(DOUBLE),
      m_use(INDEPENDENT),
      m_size(definedLength),
      m_doublePtr(&address),
      m_doubleValue(0.0),
      m_name(name)
  {
  }

  /**
   * Creates a new <b>Variable</b> instance.  The new variable will use the
   * address specified by <b>address</b> and be of type int.  This address must
   * remain in scope during the lifetime of the variable are until rebound to a new
   * address.
   *
   * @param address		an <b>int</b> reference to the value for this
   *				variable.
   *
   */
  explicit Variable(int &address, const std::string& name, unsigned definedLength=std::numeric_limits<int>::max())
    : m_type(INTEGER),
      m_use(INDEPENDENT),
      m_size(definedLength),
      m_intPtr(&address),
      m_intValue(0),
      m_name(name)
  {
  }

  /**
   * @brief Member function <b>operator=</b> assigns a new value to the variable.
   * If the variable type is integer, the value is converted to integer before
   * assignment.
   *
   * @param value		a <b>double</b> value that is to be assigned to the
   *				variable.
   *
   * @return			a <b>Variable</b> reference to the variable.
   */
  Variable &operator=(const double &value) {
    m_type = this->m_type;
    m_use = this->m_use;
    if(m_size != 1 && m_size != std::numeric_limits<int>::max()) {
      std::stringstream error;
      error << "In analytic expression evaluator, invalid use of equal on multi-component array variable '"<<m_name<<"'.  ";
      throw std::runtime_error(error.str());
    }
    if (m_type == INTEGER)
      *m_intPtr = static_cast<int>(value);
    else if (m_type == DOUBLE)
      *m_doublePtr = value;
    return *this;
  }

  /**
   * @brief Member function <b>operator=</b> assigns a new value to the variable.
   * If the variable type is double, the value is converted to double before assignment.
   *
   * @param value		an <b>int</b> value that is to be assigned to the
   *				variable.
   *
   * @return			a <b>Variable</b> reference to the variable.
   */
  Variable &operator=(const int &value) {
    m_type = this->m_type;
    m_use = this->m_use;
    if(m_size != 1 && m_size != std::numeric_limits<int>::max()) {
      std::stringstream error;
      error << "In analytic expression evaluator, invalid use of equal on multi-component variable '"<<m_name<<"'.  ";
      throw std::runtime_error(error.str());
    }
    if (m_type == INTEGER)
      *m_intPtr = value;
    else if (m_type == DOUBLE)
      *m_doublePtr = static_cast<double>(value);
    return *this;
  }

private:
  Variable(const Variable &);
  Variable &operator=(const Variable &);

public:

  void setDependent() {
    m_use = DEPENDENT;
  }
  
  bool isDependent() const {
    return m_use == DEPENDENT;
  }

  /**
   * @brief Member function <b>operator[]</b> returns a value from an array of
   * double values.  
   *
   * @param index		a <b>int</b> value of the zero based index into the
   *				array to retrieve the value.
   *
   * @return			a <b>double</b> reference to the value.
   */
  inline double& getArrayValue(int index, ArrayOffset offsetType) const {
    if (m_type != DOUBLE) {
      std::stringstream error;
      error << "In analytic expression evaluator, only double arrays allowed for variable '"<<m_name<<"'.  ";
      throw std::runtime_error(error.str());
    }

    if (m_doublePtr == nullptr) {
      std::stringstream error;
      error << "In analytic expression evaluator, unbound variable '"<<m_name<<"'.  ";
      throw std::runtime_error(error.str());
    }

    if(offsetType == ZERO_BASED_INDEX) {
      if(index < 0 || (index+1) > m_size) {
        std::stringstream error;
        error << "In analytic expression evaluator, processing variable '"<<m_name<<"'.  ";
        error << "Attempting to access invalid component '"<<index<<"' in analytic function.  Valid components are 0 to '"<<m_size-1<<"'.  ";
        throw std::runtime_error(error.str());
      }
      return m_doublePtr[index];
    } else if (offsetType == ONE_BASED_INDEX) {
      if(index < 1 || (index) > m_size) {
        std::stringstream error;
        error << "In analytic expression evaluator, processing variable '"<<m_name<<"'.  ";
        error << "Attempting to access invalid component '"<<index<<"' in analytic function.  Valid components are 1 to '"<<m_size<<"'.  ";
        throw std::runtime_error(error.str());
      }
      return m_doublePtr[index-1];
    } else {
      std::stringstream error;
      error << "In analytic expression evaluator, processing variable '"<<m_name<<"'.  ";
      error << "Invalid internal state of expression evaluator";
      throw std::runtime_error(error.str());
    }
  } 

  /**
   * @brief Member function <b>bind</b> binds variable to the address of the
   * specified value.  The type is converted to double.  This address must remain in
   * scope during the lifetime of the variable are until rebound to a new address.
   *
   * @param value_ref	a <b>double</b> reference to be used for this variable.
   *
   * @return			a <b>Variable</b> reference to the variable.
   */
  inline Variable &bind(double &value_ref, int definedLength=std::numeric_limits<int>::max()) {
    m_type = DOUBLE;
    m_doublePtr = &value_ref;
    m_size = definedLength;
    return *this;
  }

  /**
   * @brief Member function <b>bind</b> binds variable to the address of the
   * specified value.  The type is converted to int.  This address must remain in scope
   * during the lifetime of the variable are until rebound to a new address.
   *
   * @param value_ref	a <b>int</b> reference to be used for this variable.
   *
   * @return			a <b>Variable</b> reference to the variable.
   */
  inline Variable &bind(int &value_ref, int definedLength=std::numeric_limits<int>::max()) {
    m_type = INTEGER;
    m_intPtr = &value_ref;
    m_size = definedLength;
    return *this;
  }

  /**
   * @brief Member function <b>unbind</b> binds the variable to the local value
   * and reinitializes that value to zero.  The type is left unchanged.
   *
   * @return			a <b>Variable</b> reference to the variable.
   */
  inline Variable &unbind() {
    switch (m_type) {
    case DOUBLE:
      m_doublePtr = &m_doubleValue;
      m_size = 1;
      m_doubleValue = 0.0;
      break;
    case INTEGER:
      m_intPtr = &m_intValue;
      m_size = 1;
      m_intValue = 0;
      break;
    }
    return *this;
  }

  //
  //  Get the variable pointer and its defined length
  //
  double* getAddress() const {
    return m_doublePtr;
  }
  int getLength() const {
    return m_size;
  }

  /**
   * @brief Member function <b>getValue</b> returns the variable value as a double.
   *
   * @return			a <b>double</b> value of the value of the variable.
   */
  inline double getValue() const {

    if(m_size != 1 && m_size != std::numeric_limits<int>::max()) {
      std::stringstream error;
      error << "In analytic expression evaluator, processing variable '"<<m_name<<"'.  ";
      error << "Invalid direct access of array variable, must access by index";
      throw std::runtime_error(error.str());
    }

    switch (m_type) {
    case DOUBLE:
      return *m_doublePtr;
    case INTEGER:
      return *m_intPtr;
    }

    std::stringstream error;
    error << "In analytic expression evaluator, processing variable '"<<m_name<<"'.  ";
    error << "Invalid variable type";
    throw std::runtime_error(error.str());
  }

private:
  Type	        m_type;                 ///< Variable data type
  Use           m_use;                  ///< Variable is dependent or independent
  int           m_size;                 ///< Size of defined values in the double or int pointer  
  union {
    double *	m_doublePtr;		///< Pointer to value as double
    int *	m_intPtr;		///< Pointer to value as integer
  };
  union {
    double	m_doubleValue;		///< Local variable value as double
    int	        m_intValue;             ///< Local variable value as integer
  };
  const std::string  m_name;            ///< Name given to the variable, used for error messaging
};


/**
 * @brief Class <b>VariableMap</b> implements a mapping from name to a pointer to
 * a <b>Variable</b> object.  The mapping is case insensitive.
 */
class VariableMap : public std::map<std::string, Variable *, LessCase>
{
public:
  /**
   * @brief Typedef <b>value_type</b> is the value_type of the
   * <b>std::map</b> subclass.  The mapping is case insensitive.
   */
  typedef std::map<std::string, Variable *, LessCase >::value_type value_type;

  /**
   * @brief Class <b>Resolver</b> is a base class for a variable name to value
   * resolver.  After parsing of an expression, a list of variable names is collected.
   * The <b>resolve</b> member function will be called for each of these during
   * the expression resolve stage.  This function may bind a new value to the variable,
   * ignore the variable by leaving the initialized local binding, or throw an exception
   * stating the variable cannot be resolved.
   *
   */
  class Resolver
  {
  public:
    /**
     * Creates a new <b>Resolver</b> instance.
     */
    Resolver()
    {}

    /**
     * Destroys a <b>Resolver</b> instance.
     */
    virtual ~Resolver()
    {}

    /**
     * @brief Member function <b>resolve</b> is the virtual member function for
     * the variable resolver.
     *
     * @param it		a <b>VariableMap::iterator</b> reference to the
     *			variable whose name is to be resolved.
     */
    virtual void resolve(VariableMap::iterator &)
    {
    }
  };

private:
  /**
   * @brief Member function <b>delete_variable</b> is a function which will delete
   * the variable pointed to by <b>t</b>.
   *
   * @param t			a <b>value_type</b> reference to the variable to be
   *				deleted.
   */
  static void delete_variable(value_type &t) {
    delete t.second;
  }

  /**
   * @brief Class <b>DefaultResolver</b> implement a default resolver.
   */
  class DefaultResolver : public Resolver
  {
  public:
    DefaultResolver()
    {}

    /**
     * Destroys a <b>DefaultResolver</b> instance.
     */
    virtual ~DefaultResolver()
    {}
  };

public:
  /**
   * @brief Member function <b>getDefaultResolver</b> returns a reference to the default resolver.
   */
  static Resolver &getDefaultResolver();

  /**
   * Creates a new <b>VariableMap</b> instance with the specified variable name
   * resolver.
   *
   * @param resolver		a <b>Resolver</b> reference to the variable name
   *				resolver for this variable map.
   */
  VariableMap(Resolver &resolver = getDefaultResolver())
    : std::map<std::string, Variable *, LessCase>(),
      m_resolver(resolver)
  {}

  /**
   * Destroys a <b>VariableMap</b> instance.  All variables are destroyed.
   */
  virtual ~VariableMap() {
    std::for_each(begin(), end(), &delete_variable);
  }

  /**
   * @brief Member function <b>operator[]</b> ...
   *
   * @param s			a <b>std::string</b> const reference to the variable's
   *				name.
   *
   * @return			a <b>Variable</b> pointer to the new variable.
   */
  Variable *operator[](const std::string &s) {
    std::pair<iterator,bool> i = insert(std::pair<const std::string, Variable *>(s, (Variable*)nullptr));
    if (i.second) {
      (*i.first).second = new Variable(s);
    }
    return (*i.first).second;
  }

  /**
   * @brief Member function <b>getResolver</b> returns a reference to the name
   * resolver.
   *
   * @return			a <b>Resolver</b> reference to the variable name
   *				resolver.
   */
  Resolver &getResolver() {
    return m_resolver;
  }

private:
  Resolver &		m_resolver;		///< Reference to this variable map name resolver
};

} // namespace expreval
} // namespace stk

#endif // stk_expreval_Variable_hpp
