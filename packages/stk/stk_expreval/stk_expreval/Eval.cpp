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

#include "stk_expreval/Eval.hpp"
#include "stk_expreval/Lexer.hpp"
#include "stk_expreval/Parser.hpp"
#include "stk_expreval/ParsedEval.hpp"
#include <queue>

namespace stk {
namespace expreval {

Eval::Eval(VariableMap::Resolver & resolver, const std::string & expression, Variable::ArrayOffset arrayOffsetType)
  : m_variableMap(resolver),
    m_expression(expression),
    m_syntaxStatus(false),
    m_parseStatus(false),
    m_headNode(nullptr),
    m_arrayOffsetType(arrayOffsetType),
    m_parsedEval(nullptr)
{
  initialize_function_map();
}

Eval::Eval(const std::string & expression, Variable::ArrayOffset arrayOffsetType)
  : m_variableMap(VariableMap::getDefaultResolver()),
    m_expression(expression),
    m_syntaxStatus(false),
    m_parseStatus(false),
    m_headNode(nullptr),
    m_arrayOffsetType(arrayOffsetType),
    m_parsedEval(nullptr)
{
  initialize_function_map();
}

Eval::Eval(const Eval& otherEval)
  : m_variableMap(otherEval.m_variableMap),
    m_undefinedFunctionSet(otherEval.m_undefinedFunctionSet),
    m_functionMap(otherEval.m_functionMap),
    m_expression(otherEval.m_expression),
    m_syntaxStatus(otherEval.m_syntaxStatus),
    m_parseStatus(otherEval.m_parseStatus),
    m_headNode(otherEval.m_headNode),
    m_nodes(otherEval.m_nodes),
    m_evaluationNodes(otherEval.m_evaluationNodes),
    m_arrayOffsetType(otherEval.m_arrayOffsetType),
    m_resultBuffer(otherEval.m_resultBuffer),
    m_parsedEval(nullptr)
{}

Eval::~Eval()
{
  delete m_parsedEval;
}

Node*
Eval::newNode(int opcode)
{
  m_nodes.push_back(std::make_shared<Node>(static_cast<Opcode>(opcode), this));
  m_nodes.back().get()->m_currentNodeIndex = m_nodes.size() - 1;
  return m_nodes.back().get();
}

void
Eval::syntax()
{
  m_syntaxStatus = false;
  m_parseStatus = false;

  try {
    // Validate the characters
    LexemVector lex_vector = tokenize(m_expression);

    // Call the multiparse routine to parse subexpressions
    m_headNode = Parser::parseStatements(*this, lex_vector.begin(), lex_vector.end());

    m_syntaxStatus = true;
  }
  catch (std::runtime_error &) {
  }
}

void
Eval::ternaryFixup()
{
  for (Node * node : m_evaluationNodes) {
    if (node->m_left && node->m_left->m_opcode == OPCODE_TERNARY_PREDICATE) {
      Node * joinNode = node->m_left->m_ternaryOther;
      node->m_left = joinNode;
    }

    if (node->m_right && node->m_right->m_opcode == OPCODE_TERNARY_PREDICATE) {
      Node * joinNode = node->m_right->m_ternaryOther;
      node->m_right = joinNode;
    }
  }
}

void
Eval::parse()
{
  try {
    syntax();

    if (m_syntaxStatus) {
      if (!m_undefinedFunctionSet.empty()) {
        std::ostringstream strout;
        strout << "In expression '" << m_expression << "', the following functions are not defined:" << std::endl;
        for (const auto& it : m_undefinedFunctionSet)
          strout << it << std::endl;
        throw std::runtime_error(strout.str());
      }

      resolve();

      if (m_headNode) {
        NodeWeightMap nodeWeights;
        m_headNode->computeNodeWeight(nodeWeights);
        m_headNode->evalTrace(nodeWeights, m_evaluationNodes);
        ternaryFixup();

        assign_result_buffer_indices();
      }

      m_parseStatus = true;
    } else {
      throw std::runtime_error("The following expression has a syntax error in it.\n" + m_expression);
    }
  }
  catch (std::runtime_error & ) {
    throw;
  }

  int index = 0;
  for (VariableMap::iterator it = m_variableMap.begin(); it != m_variableMap.end(); ++it) {
    it->second->set_index(index);
    ++index;
  }

}

void
Eval::resolve()
{
  for (VariableMap::iterator it = m_variableMap.begin(); it != m_variableMap.end(); ++it) {
    m_variableMap.getResolver().resolve(it);
  }
}

double
Eval::evaluate() const
{
  if (!m_parseStatus) {
    throw std::runtime_error(std::string("Expression '") + m_expression + "' did not parse successfully");
  }

  double returnValue = 0.0;
  try
  {

    if (m_headNode) {
      int nodeIndex = m_evaluationNodes.front()->m_currentNodeIndex;
      while (nodeIndex >= 0) {
        Node * node = m_nodes[nodeIndex].get();
        node->eval();
        nodeIndex = node->getNextNodeIndex();
      }
      returnValue = m_evaluationNodes.back()->getResult();
    }

  }
  catch(expression_evaluation_exception &)
  {
    throw std::runtime_error(std::string("Expression '") + m_expression + "' did not evaluate successfully");
  }
  return returnValue;
}

bool
Eval::undefinedFunction() const
{
  /* Check for an undefined function in any allocated node */
  return !m_undefinedFunctionSet.empty();
}

bool
Eval::is_constant_expression() const
{
  return m_variableMap.empty();
}

bool
Eval::is_variable(const std::string& variableName) const
{
  return (m_variableMap.count(variableName) > 0);
}

bool
Eval::is_scalar(const std::string& variableName) const
{
  auto variableIterator = m_variableMap.find(variableName);

  if (variableIterator == m_variableMap.end()) {
    return false;
  }

  int variableLength = variableIterator->second->getLength();
  return variableLength == 1 || variableLength == std::numeric_limits<int>::max();
}

std::vector<std::string>
Eval::get_variable_names() const
{
  std::vector<std::string> variableList;
  for(auto& currentVariable : m_variableMap) {
    std::string variableName = currentVariable.first;
    variableList.push_back(variableName);
  }

  return variableList;
}

std::vector<std::string>
Eval::get_dependent_variable_names() const
{
  std::vector<std::string> dependentVariableList;
  for(auto& currentVariable : m_variableMap) {
    std::string variableName = currentVariable.first;
    stk::expreval::Variable* variable = currentVariable.second.get();
    if (variable->isDependent()) {
      dependentVariableList.push_back(variableName);
    }
  }

  return dependentVariableList;
}

std::vector<std::string>
Eval::get_independent_variable_names() const
{
  std::vector<std::string> independentVariableList;
  for(auto& currentVariable : m_variableMap) {
    std::string variableName = currentVariable.first;
    stk::expreval::Variable* variable = currentVariable.second.get();
    if (!(variable->isDependent())) {
      independentVariableList.push_back(variableName);
    }
  }

  return independentVariableList;
}

int
Eval::get_variable_index(const std::string & variable) const
{
  const auto variableIter = m_variableMap.find(variable);
  STK_ThrowRequireMsg(variableIter != m_variableMap.end(), "Variable " + variable + " Not Found in VariableMap");
  return variableIter->second->get_index();
}

int
Eval::get_head_node_index() const
{
  return (m_headNode) ? m_headNode->m_currentNodeIndex : -1;
}

int
Eval::get_first_node_index() const
{
  return (!m_evaluationNodes.empty()) ? m_evaluationNodes.front()->m_currentNodeIndex : -1;
}

int
Eval::get_last_node_index() const
{
  return (!m_evaluationNodes.empty()) ? m_evaluationNodes.back()->m_currentNodeIndex : -1;
}

FunctionType
Eval::get_function_type(const std::string& functionName) const
{
  const auto functionIt = m_functionMap.find(functionName);
  if (functionIt != m_functionMap.end()) {
    return functionIt->second;
  }
  else {
    return FunctionType::UNDEFINED;
  }
}

void
Eval::initialize_function_map()
{
  m_functionMap["abs"] = FunctionType::ABS;
  m_functionMap["fabs"] = FunctionType::ABS;
  m_functionMap["max"] = FunctionType::MAX;
  m_functionMap["min"] = FunctionType::MIN;
  m_functionMap["sign"] = FunctionType::SIGN;
  m_functionMap["ipart"] = FunctionType::IPART;
  m_functionMap["fpart"] = FunctionType::FPART;
  m_functionMap["ceil"] = FunctionType::CEIL;
  m_functionMap["floor"] = FunctionType::FLOOR;
  m_functionMap["mod"] = FunctionType::MOD;
  m_functionMap["fmod"] = FunctionType::MOD;
  m_functionMap["pow"] = FunctionType::POW;
  m_functionMap["sqrt"] = FunctionType::SQRT;
  m_functionMap["exp"] = FunctionType::EXP;
  m_functionMap["ln"] = FunctionType::LN;
  m_functionMap["log"] = FunctionType::LN;
  m_functionMap["log10"] = FunctionType::LOG10;
  m_functionMap["deg"] = FunctionType::DEG;
  m_functionMap["rad"] = FunctionType::RAD;
  m_functionMap["sin"] = FunctionType::SIN;
  m_functionMap["cos"] = FunctionType::COS;
  m_functionMap["tan"] = FunctionType::TAN;
  m_functionMap["asin"] = FunctionType::ASIN;
  m_functionMap["acos"] = FunctionType::ACOS;
  m_functionMap["atan"] = FunctionType::ATAN;
  m_functionMap["atan2"] = FunctionType::ATAN2;
  m_functionMap["sinh"] = FunctionType::SINH;
  m_functionMap["cosh"] = FunctionType::COSH;
  m_functionMap["tanh"] = FunctionType::TANH;
  m_functionMap["asinh"] = FunctionType::ASINH;
  m_functionMap["acosh"] = FunctionType::ACOSH;
  m_functionMap["atanh"] = FunctionType::ATANH;
  m_functionMap["erf"] = FunctionType::ERF;
  m_functionMap["erfc"] = FunctionType::ERFC;
  m_functionMap["poltorectx"] = FunctionType::POLTORECTX;
  m_functionMap["poltorecty"] = FunctionType::POLTORECTY;
  m_functionMap["recttopolr"] = FunctionType::RECTTOPOLR;
  m_functionMap["recttopola"] = FunctionType::RECTTOPOLA;

  m_functionMap["unit_step"] = FunctionType::UNIT_STEP;
  m_functionMap["cycloidal_ramp"] = FunctionType::CYCLOIDAL_RAMP;
  m_functionMap["cos_ramp"] = FunctionType::COS_RAMP;
  m_functionMap["cosine_ramp"] = FunctionType::COS_RAMP;
  m_functionMap["linear_ramp"] = FunctionType::LINEAR_RAMP;
  m_functionMap["haversine_pulse"] = FunctionType::HAVERSINE_PULSE;
  m_functionMap["point2d"] = FunctionType::POINT2D;
  m_functionMap["point3d"] = FunctionType::POINT3D;

  m_functionMap["exponential_pdf"] = FunctionType::EXPONENTIAL_PDF;
  m_functionMap["log_uniform_pdf"] = FunctionType::LOG_UNIFORM_PDF;
  m_functionMap["normal_pdf"] = FunctionType::NORMAL_PDF;
  m_functionMap["weibull_pdf"] = FunctionType::WEIBULL_PDF;
  m_functionMap["gamma_pdf"] = FunctionType::GAMMA_PDF;


  m_functionMap["ts_random"] = FunctionType::TS_RANDOM;
  m_functionMap["ts_normal"] = FunctionType::TS_NORMAL;
}

Eval &
Eval::setExpression(const std::string & expression)
{
  m_expression = expression;
  m_parseStatus = false;
  return *this;
}

class ResultBufferIndices
{
public:
  int get_free_index()
  {
    if  (m_freeList.size() == 0) {
      release_index(m_resultBufferSize++);
    }

    auto idx = m_freeList.front();
    m_freeList.pop();
    return idx;
  }

  void release_index(const int idx)
  {
    STK_ThrowRequireMsg(idx >= 0, "Attempting to free negative index");
    m_freeList.push(idx);
  }

  int get_result_buffer_size() const { return m_resultBufferSize;}

private:
  std::queue<int> m_freeList;
  int m_resultBufferSize = 0;
};

int
Eval::assign_result_buffer_indices()
{
  ResultBufferIndices indexAssigner;

  for (auto node : m_evaluationNodes) {
    if (node->m_left && node->m_left->m_resultIdx >= 0) {
      indexAssigner.release_index(node->m_left->m_resultIdx);
    }

    if (node->m_right && node->m_right->m_resultIdx >= 0) {
      indexAssigner.release_index(node->m_right->m_resultIdx);
    }

    if (node->m_ternaryOther && node->m_ternaryOther->m_resultIdx >= 0) {
      indexAssigner.release_index(node->m_ternaryOther->m_resultIdx);
    }

    node->m_resultIdx = indexAssigner.get_free_index();
  }

  m_resultBuffer.resize(indexAssigner.get_result_buffer_size());
  return indexAssigner.get_result_buffer_size();
}

Eval &
Eval::setValue(const std::string &name, double* value, int definedLength)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->bind(*value, definedLength);
  }
  return *this;
}

double
Eval::getValue(const std::string &name)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it == m_variableMap.end()) {
    throw std::runtime_error(std::string("Variable ") + name  + " not defined");
  }
  return (*it).second->getValue();
}

Eval &
Eval::bindVariable(const std::string &name, const double &value_ref, int definedLength)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->bind(value_ref, definedLength);
  }
  return *this;
}

Eval &
Eval::bindVariable(const std::string &name, double &value_ref, int definedLength)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->bind(value_ref, definedLength);
  }
  return *this;
}

Eval &
Eval::bindVariable(const std::string &name, const int &value_ref, int definedLength)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->bind(value_ref, definedLength);
  }
  return *this;
}

Eval &
Eval::bindVariable(const std::string &name, int &value_ref, int definedLength)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->bind(value_ref, definedLength);
  }
  return *this;
}

Eval &
Eval::unbindVariable(const std::string &name)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->unbind();
  }
  return *this;
}

Eval &
Eval::deactivateVariable(const std::string &name)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it != m_variableMap.end()) {
    (*it).second->deactivate();
  }
  return *this;
}

Variable &
Eval::getVariable(const std::string &name)
{
  VariableMap::iterator it = m_variableMap.find(name);
  if (it == m_variableMap.end()) {
    throw std::runtime_error(std::string("Variable ") + name  + " not defined");
  }
  return *(*it).second;
}

void
Eval::syntaxCheck(const std::string &expr)
{
  setExpression(expr);
  syntax();
}

void
Eval::parse(const std::string &expr)
{
  setExpression(expr);
  parse();
}

}
}
