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

#ifndef stk_expreval_evaluator_hpp
#define stk_expreval_evaluator_hpp

#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>
#include <memory>

#include <stk_expreval/Function.hpp>
#include <stk_expreval/Constants.hpp>
#include <stk_expreval/Variable.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <Kokkos_Core.hpp>
#include<stk_util/ngp/NgpSpaces.hpp>

namespace stk {
namespace expreval {

class Node;
class NgpNode;
class ParsedEval;
template <int NUMVARS=16>
class DeviceVariableMap;

enum class FunctionType {

  ABS,
  MAX,
  MIN,
  SIGN,
  IPART,
  FPART,
  CEIL,
  FLOOR,
  MOD,
  POW,
  SQRT,
  EXP,
  LN,
  LOG10,

  DEG,
  RAD,
  SIN,
  COS,
  TAN,
  ASIN,
  ACOS,
  ATAN,
  ATAN2,
  SINH,
  COSH,
  TANH,
  ASINH,
  ACOSH,
  ATANH,
  ERF,
  ERFC,
  POLTORECTX,
  POLTORECTY,
  RECTTOPOLR,
  RECTTOPOLA,

  UNIT_STEP,
  CYCLOIDAL_RAMP,
  COS_RAMP,
  HAVERSINE_PULSE,
  POINT2D,
  POINT3D,

  EXPONENTIAL_PDF,
  LOG_UNIFORM_PDF,
  NORMAL_PDF,
  WEIBULL_PDF,
  GAMMA_PDF,

  RAND,
  SRAND,
  RANDOM,
  TS_RANDOM,
  TS_NORMAL,
  TIME,
  
  UNDEFINED

};

class Eval
{
public:
  typedef std::set<std::string> UndefinedFunctionSet;

  Eval(VariableMap::Resolver &resolver = VariableMap::getDefaultResolver(), const std::string &expr = "", const Variable::ArrayOffset arrayOffsetType = Variable::ZERO_BASED_INDEX);
  explicit Eval(const std::string &expr, const Variable::ArrayOffset arrayOffsetType = Variable::ZERO_BASED_INDEX);

  Eval(const Eval &);

private:
  Eval &operator=(const Eval &);

public:
  ~Eval();

  const std::string &getExpression() const {
    return m_expression;
  }

  Eval &setExpression(const std::string &expression) {
    m_expression = expression;
    m_parseStatus = false;
    return *this;
  }

  VariableMap &getVariableMap() {
    return m_variableMap;
  }

  bool is_constant_expression() const;

  bool is_variable(const std::string& variableName) const; 

  bool is_scalar(const std::string& variableName) const;
 
  std::vector<std::string> get_variable_names() const;

  std::vector<std::string> get_dependent_variable_names() const;

  std::vector<std::string> get_independent_variable_names() const;

  int get_variable_index(const std::string& variableName) const;

  int get_variable_count() const { return m_variableMap.size(); }

  int get_node_count() const { return m_nodes.size(); }

  int get_head_node_index() const;

  Node* get_node(int i) const { return m_nodes[i].get(); }
 
  KOKKOS_FUNCTION 
  int get_num_variables() const
  {
    return m_variableMap.size();
  }
  
  UndefinedFunctionSet &getUndefinedFunctionSet() {
    return m_undefinedFunctionSet;
  }

  bool getSyntaxStatus() const {
    return m_syntaxStatus;
  }

  bool getParseStatus() const {
    return m_parseStatus;
  }

  Eval &setValue(const std::string &name, double* value, int definedLength) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it != m_variableMap.end()) {
      (*it).second->bind(*value, definedLength);
    }
    return *this;
  }

  double getValue(const std::string &name) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it == m_variableMap.end())
      throw std::runtime_error(std::string("Variable ") + name  + " not defined");
    return (*it).second->getValue();
  }

  Node *newNode(int op);
  
  Eval &bindVariable(const std::string &name, double &value_ref, int definedLength=std::numeric_limits<int>::max()) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it != m_variableMap.end()) {
      (*it).second->bind(value_ref, definedLength);
    }
    return *this;
  }

  Variable &getVariable(const std::string &name) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it == m_variableMap.end())
      throw std::runtime_error(std::string("Variable ") + name  + " not defined");
    return *(*it).second;
  }

  void syntaxCheck(const std::string &expr) {
    setExpression(expr);
    syntax();
  }

  void syntax();

  void parse(const std::string &expr) {
    setExpression(expr);
    parse();
  }

  void parse();

  void resolve();

  double evaluate() const;

  bool undefinedFunction() const;

  Variable::ArrayOffset getArrayOffsetType() {return m_arrayOffsetType;}

  ParsedEval& get_parsed_eval();

  FunctionType get_function_type(const std::string& functionName) const;

private:

  void initialize_function_map();

  VariableMap		m_variableMap;		
  UndefinedFunctionSet	m_undefinedFunctionSet;	
  std::map<std::string, FunctionType, LessCase> m_functionMap;

  std::string		m_expression;	
  bool			m_syntaxStatus;
  bool			m_parseStatus;

  Node*		m_headNode;
  std::vector<std::shared_ptr<Node>> m_nodes;                
  Variable::ArrayOffset m_arrayOffsetType;

  ParsedEval* m_parsedEval;
};

enum Opcode {
  OPCODE_UNDEFINED,
  OPCODE_CONSTANT,
  OPCODE_RVALUE,
  OPCODE_STATEMENT,
  OPCODE_ARGUMENT,

  OPCODE_TERNARY,
  
  OPCODE_MULTIPLY,
  OPCODE_DIVIDE,
  OPCODE_MODULUS,
  OPCODE_ADD,
  OPCODE_SUBTRACT,
  OPCODE_UNARY_MINUS,
  OPCODE_FUNCTION,

  OPCODE_EQUAL,
  OPCODE_NOT_EQUAL,
  OPCODE_LESS,
  OPCODE_GREATER,
  OPCODE_LESS_EQUAL,
  OPCODE_GREATER_EQUAL,

  OPCODE_UNARY_NOT,
  OPCODE_LOGICAL_AND,
  OPCODE_LOGICAL_OR,

  OPCODE_EXPONENIATION,

  OPCODE_ASSIGN
};

class NgpNode
{

public:
  enum { MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES = 5 };
  enum { MAXIMUM_FUNCTION_NAME_LENGTH = 32 };

  KOKKOS_FUNCTION
  NgpNode();

  KOKKOS_FUNCTION
  explicit NgpNode(const Node& node);

  KOKKOS_DEFAULTED_FUNCTION
  NgpNode(const NgpNode &) = default;

  KOKKOS_DEFAULTED_FUNCTION
  NgpNode &operator=(const NgpNode &) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~NgpNode() = default;

  template <int NUMVARS>
  KOKKOS_FUNCTION
  double
  eval(DeviceVariableMap<NUMVARS>& deviceVariableMap) const
  {
    switch (m_opcode) {
    case OPCODE_STATEMENT:
      {
        double value = 0.0;
        for (const NgpNode* statement = this; statement; statement = statement->get_right_node()) {
          value = statement->get_left_node()->eval(deviceVariableMap);
        }
        return value;
      }
  
    case OPCODE_CONSTANT:
      return m_data.constant.value;
  
    case OPCODE_RVALUE:
      if (get_left_node()) {
        return deviceVariableMap[m_data.variable.variableIndex].getArrayValue(get_left_node()->eval(deviceVariableMap), deviceVariableMap.get_array_offset_type());
      } 
      else {
        return deviceVariableMap[m_data.variable.variableIndex].getValue();
      }
  
    case OPCODE_MULTIPLY:
      return get_left_node()->eval(deviceVariableMap)*get_right_node()->eval(deviceVariableMap);
  
    case OPCODE_EXPONENIATION:
      return std::pow(get_left_node()->eval(deviceVariableMap), get_right_node()->eval(deviceVariableMap));
  
    case OPCODE_DIVIDE:
      return get_left_node()->eval(deviceVariableMap)/get_right_node()->eval(deviceVariableMap);
  
    case OPCODE_MODULUS:
      return std::fmod(get_left_node()->eval(deviceVariableMap), get_right_node()->eval(deviceVariableMap));
  
    case OPCODE_ADD:
      return get_left_node()->eval(deviceVariableMap) + get_right_node()->eval(deviceVariableMap);
      
    case OPCODE_SUBTRACT:
      return get_left_node()->eval(deviceVariableMap) - get_right_node()->eval(deviceVariableMap);
  
    case OPCODE_EQUAL:
      return get_left_node()->eval(deviceVariableMap) == get_right_node()->eval(deviceVariableMap) ? 1.0 : 0.0;
  
    case OPCODE_NOT_EQUAL:
      return get_left_node()->eval(deviceVariableMap) != get_right_node()->eval(deviceVariableMap) ? 1.0 : 0.0;
  
    case OPCODE_LESS:
      return get_left_node()->eval(deviceVariableMap) < get_right_node()->eval(deviceVariableMap) ? 1.0 : 0.0;
  
    case OPCODE_GREATER:
      return get_left_node()->eval(deviceVariableMap) > get_right_node()->eval(deviceVariableMap) ? 1.0 : 0.0;
  
    case OPCODE_LESS_EQUAL:
      return get_left_node()->eval(deviceVariableMap) <= get_right_node()->eval(deviceVariableMap) ? 1.0 : 0.0;
  
    case OPCODE_GREATER_EQUAL:
      return get_left_node()->eval(deviceVariableMap) >= get_right_node()->eval(deviceVariableMap) ? 1.0 : 0.0;
  
    case OPCODE_LOGICAL_AND: {
      double left = get_left_node()->eval(deviceVariableMap);
      double right = get_right_node()->eval(deviceVariableMap);
      return  (left != 0.0) && (right != 0.0) ? 1.0 : 0.0;
    }
    case OPCODE_LOGICAL_OR: {
      double left = get_left_node()->eval(deviceVariableMap);
      double right = get_right_node()->eval(deviceVariableMap);
      return (left != 0.0) || (right != 0.0) ? 1.0 : 0.0;
    }
    case OPCODE_TERNARY:
      return get_left_node()->eval(deviceVariableMap) != 0.0 ? get_right_node()->eval(deviceVariableMap) : get_other_node()->eval(deviceVariableMap);
  
    case OPCODE_UNARY_MINUS:
      return -get_right_node()->eval(deviceVariableMap);
  
    case OPCODE_UNARY_NOT:
      return get_right_node()->eval(deviceVariableMap) == 0.0 ? 1.0 : 0.0;
  
    case OPCODE_ASSIGN:
      if (get_left_node())
        return deviceVariableMap[m_data.variable.variableIndex].getArrayValue(get_left_node()->eval(deviceVariableMap), deviceVariableMap.get_array_offset_type()) = get_right_node()->eval(deviceVariableMap);
      else {
        deviceVariableMap[m_data.variable.variableIndex] = get_right_node()->eval(deviceVariableMap);
        return deviceVariableMap[m_data.variable.variableIndex].getValue();
      }
  
    case OPCODE_FUNCTION:
      {
        double arguments[20];
  
        int argumentCount = 0;
        for (const NgpNode* arg = get_right_node(); arg; arg = arg->get_right_node())
        {
          arguments[argumentCount++] = arg->get_left_node()->eval(deviceVariableMap);
        }

        return evaluate_function(argumentCount, arguments); 
      }
  
    default: // Unknown opcode
    return 0.0;
    }
    
    return 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  const NgpNode* get_left_node() const {
    return m_leftNodeIndex != -1 ? (this + (m_leftNodeIndex - m_currentNodeIndex)) : nullptr; 
  }

  KOKKOS_INLINE_FUNCTION
  const NgpNode* get_right_node() const {
    return m_rightNodeIndex != -1 ? (this + (m_rightNodeIndex - m_currentNodeIndex)) : nullptr; 
  }

  KOKKOS_INLINE_FUNCTION
  const NgpNode* get_other_node() const {
    return m_otherNodeIndex != -1 ? (this + (m_otherNodeIndex - m_currentNodeIndex)) : nullptr; 
  }

  Opcode m_opcode;

  union _data
  {
    struct _constant
    {
      double value;
    } constant;

    struct _variable
    {
      int variableIndex;
      Variable::Type variableType;
      int variableSize;
    } variable;

    struct _function
    {
      FunctionType functionType;
    } function;
  } m_data;

  int m_currentNodeIndex;
  int m_leftNodeIndex;
  int m_rightNodeIndex;
  int m_otherNodeIndex;
 
  KOKKOS_FUNCTION
  double evaluate_function(int argumentCount, double* arguments) const;
};

class ParsedEval 
{

  using NodeView = Kokkos::View<NgpNode*, stk::ngp::MemSpace>;

public:
  ParsedEval(Eval& eval); 

  KOKKOS_DEFAULTED_FUNCTION ParsedEval(const ParsedEval&) = default;

  KOKKOS_DEFAULTED_FUNCTION ~ParsedEval() = default;

  KOKKOS_INLINE_FUNCTION
  int get_num_variables() const { return m_numVariables; }

  KOKKOS_INLINE_FUNCTION
  Variable::ArrayOffset get_array_offset_type() const { return m_arrayOffsetType; }

  template <int NUMVARS>
  KOKKOS_INLINE_FUNCTION
  double
  evaluate(DeviceVariableMap<NUMVARS>& deviceVariableMap) const
  {
    return m_deviceNodes[m_headNodeIndex].eval(deviceVariableMap);
  }

private:
  template <int NUMVARS>
  friend class DeviceVariableMap;

  int m_numVariables;
  Variable::ArrayOffset m_arrayOffsetType;
  int m_headNodeIndex;
  NodeView m_deviceNodes;
  NodeView::HostMirror m_hostNodes;
  
};

template <int NUMVARS>
class DeviceVariableMap 
{
  using VariableMapView = Kokkos::View<DeviceVariable*, stk::ngp::MemSpace>;

public:

  KOKKOS_DEFAULTED_FUNCTION ~DeviceVariableMap() = default;

  KOKKOS_INLINE_FUNCTION
  explicit DeviceVariableMap(const ParsedEval& parsedEval) 
    : m_arrayOffsetType(parsedEval.m_arrayOffsetType) 
  {
    NGP_ThrowRequireMsg(parsedEval.get_num_variables() <= NUMVARS, "Size of DeviceVariableMap is too small");

    const ParsedEval::NodeView& deviceNodes = parsedEval.m_deviceNodes;
    for (unsigned nodeIndex = 0u; nodeIndex < deviceNodes.extent(0); ++nodeIndex)
    {
      if (deviceNodes(nodeIndex).m_opcode == OPCODE_ASSIGN) {
        m_deviceVariableMap[deviceNodes(nodeIndex).m_data.variable.variableIndex] = DeviceVariable(deviceNodes(nodeIndex).m_data.variable.variableType, deviceNodes(nodeIndex).m_data.variable.variableSize); 
      }
    }   
  }

  KOKKOS_INLINE_FUNCTION
  void bind(int variableIndex, double& value_ref, int definedLength=1, int strideLength=1) {
    m_deviceVariableMap[variableIndex].bind(value_ref, definedLength, strideLength);
  }

  KOKKOS_INLINE_FUNCTION
  void bind(int variableIndex, int& value_ref, int definedLength=1, int strideLength=1) {
    m_deviceVariableMap[variableIndex].bind(value_ref, definedLength, strideLength);
  }

  KOKKOS_INLINE_FUNCTION
  DeviceVariable& operator[](int index) { return m_deviceVariableMap[index]; }

  KOKKOS_INLINE_FUNCTION
  Variable::ArrayOffset get_array_offset_type() { return m_arrayOffsetType; } 

private:
  DeviceVariable m_deviceVariableMap[NUMVARS]; 
  Variable::ArrayOffset m_arrayOffsetType;

};
} // namespace expreval
} // namespace stk

#endif //stk_expreval_evaluator_hpp
