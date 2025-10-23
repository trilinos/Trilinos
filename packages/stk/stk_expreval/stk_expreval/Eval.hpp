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

#ifndef EVAL_HPP
#define EVAL_HPP

#include "stk_expreval/Variable.hpp"
#include "stk_expreval/Function.hpp"
#include "stk_expreval/Node.hpp"
#include "stk_expreval/ParsedEvalBase.hpp"
#include <set>
#include <string>
#include <vector>
#include <memory>

namespace stk {
namespace expreval {

class Node;

template <int RESULT_BUFFER_SIZE=DEFAULT_RESULT_BUFFER_SIZE>
class ParsedEval;

class Eval
{
public:
  typedef std::set<std::string> UndefinedFunctionSet;

  enum class FPErrorBehavior {
    Ignore,
    Warn,
    WarnOnce,
    Error
  };

  Eval(VariableMap::Resolver &resolver = VariableMap::getDefaultResolver(),
       const std::string &expr = "",
       const Variable::ArrayOffset arrayOffsetType = Variable::ZERO_BASED_INDEX);

  explicit Eval(const std::string &expr,
                const Variable::ArrayOffset arrayOffsetType = Variable::ZERO_BASED_INDEX);

  Eval(const Eval &);

private:
  Eval &operator=(const Eval &);

public:
  ~Eval();

  const std::string &getExpression() const { return m_expression; }

  Eval &setExpression(const std::string &expression);

  VariableMap &getVariableMap() { return m_variableMap; }

  int assign_result_buffer_indices();

  bool is_constant_expression() const;

  bool is_variable(const std::string& variableName) const;

  bool is_scalar(const std::string& variableName) const;

  std::vector<std::string> get_variable_names() const;

  std::vector<std::string> get_dependent_variable_names() const;

  std::vector<std::string> get_independent_variable_names() const;

  int get_variable_index(const std::string& variableName) const;

  int get_variable_count() const { return m_variableMap.size(); }

  int get_node_count() const { return m_nodes.size(); }

  int get_result_buffer_size() const { return m_resultBuffer.size(); }

  int get_head_node_index() const;

  int get_first_node_index() const;

  int get_last_node_index() const;

  Node* get_node(int i) const { return m_nodes[i].get(); }

  KOKKOS_FUNCTION
  int get_num_variables() const { return m_variableMap.size(); }

  UndefinedFunctionSet &getUndefinedFunctionSet() { return m_undefinedFunctionSet; }

  void set_fp_error_behavior(FPErrorBehavior flag);

  FPErrorBehavior get_fp_error_behavior() const { return m_fpErrorBehavior; }
  
  void set_fp_warning_issued() { m_fpWarningIssued = true; }
  
  bool get_fp_warning_issued() const { return m_fpWarningIssued; }

  bool getSyntaxStatus() const { return m_syntaxStatus; }

  bool getParseStatus() const { return m_parseStatus; }

  Eval &setValue(const std::string &name, double* value, int definedLength);

  double getValue(const std::string &name);

  Node *newNode(int op);

  Eval &bindVariable(const std::string &name, const double &value_ref,
                     int definedLength=std::numeric_limits<int>::max());

  Eval &bindVariable(const std::string &name, double &value_ref,
                     int definedLength=std::numeric_limits<int>::max());

  Eval &bindVariable(const std::string &name, const int &value_ref,
                     int definedLength=std::numeric_limits<int>::max());

  Eval &bindVariable(const std::string &name, int &value_ref,
                     int definedLength=std::numeric_limits<int>::max());

  Eval &unbindVariable(const std::string &name);

  Eval &deactivateVariable(const std::string &name);

  Variable &getVariable(const std::string &name);

  void syntaxCheck(const std::string &expr);

  void syntax();

  void ternaryFixup();

  void parse(const std::string &expr);

  void parse();

  void resolve();

  double evaluate() const;

  bool undefinedFunction() const;

  Variable::ArrayOffset getArrayOffsetType() const {return m_arrayOffsetType;}

  template <int RESULT_BUFFER_SIZE=DEFAULT_RESULT_BUFFER_SIZE>
  ParsedEval<RESULT_BUFFER_SIZE> &
  get_parsed_eval()
  {
    if (m_parsedEval == nullptr) {
      m_parsedEval = new ParsedEval<RESULT_BUFFER_SIZE>(*this);
    }

    if (m_parsedEval->get_result_buffer_size() != RESULT_BUFFER_SIZE) {
      delete m_parsedEval;
      m_parsedEval = new ParsedEval<RESULT_BUFFER_SIZE>(*this);
    }
    return *static_cast<ParsedEval<RESULT_BUFFER_SIZE>*>(m_parsedEval);
  }

  template <int RESULT_BUFFER_SIZE=DEFAULT_RESULT_BUFFER_SIZE>
 const ParsedEval<RESULT_BUFFER_SIZE>
  get_standalone_parsed_eval() const
  {
    return ParsedEval<RESULT_BUFFER_SIZE>(*this);
  }

  FunctionType get_function_type(const std::string& functionName) const;

  double& get_result_buffer_value(const int idx) { return m_resultBuffer[idx];}

private:
  void print_expression_if_fp_warning(bool fpWarningPreviouslyIssued) const;
  
  friend void check_node_order(const std::string & expression);
  friend void check_evaluation_node_order(const std::string & expression);

  void initialize_function_map();

  VariableMap m_variableMap;
  UndefinedFunctionSet m_undefinedFunctionSet;
  std::map<std::string, FunctionType, LessCase> m_functionMap;

  std::string m_expression;
  bool m_syntaxStatus;
  bool m_parseStatus;
  FPErrorBehavior m_fpErrorBehavior;
  mutable bool m_fpWarningIssued;

  Node* m_headNode;
  std::vector<std::shared_ptr<Node>> m_nodes;
  EvalNodesType m_evaluationNodes;
  Variable::ArrayOffset m_arrayOffsetType;
  std::vector<double> m_resultBuffer;

  ParsedEvalBase * m_parsedEval;
};

Eval::FPErrorBehavior fp_error_behavior_string_to_enum(const std::string& str);

}
}

#endif // EVAL_HPP
