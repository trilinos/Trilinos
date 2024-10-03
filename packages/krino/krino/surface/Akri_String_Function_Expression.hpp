// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_STRING_FUNCTION_EXPRESSION_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_STRING_FUNCTION_EXPRESSION_HPP_

#include <stk_expreval/Eval.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

class String_Function_Expression : public stk::expreval::VariableMap::Resolver
{
public:
  String_Function_Expression(const std::string & expression);
  String_Function_Expression (const String_Function_Expression&) { throw std::runtime_error("copying String_Function_Expression not allowed because it will not be resolved correctly."); }
  String_Function_Expression& operator= (const String_Function_Expression&) { throw std::runtime_error("copying String_Function_Expression not allowed because it will not be resolved correctly."); }
  void resolve(stk::expreval::VariableMap::iterator & varIt) override;
  double evaluate(const stk::math::Vector3d &coords) const;
  double evaluate(const double time, const stk::math::Vector3d &coord) const;
private:
  void parse(const std::string & expression);
  stk::expreval::Eval myEvaluator;
  mutable bool myDoesUseTime{false};
  mutable double myTime{0.0};
  mutable stk::math::Vector3d myQueryCoords{stk::math::Vector3d::ZERO};
};

void initialize_expression_vector(const std::vector<std::string> & stringVec, std::vector<String_Function_Expression> & exprVec);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_STRING_FUNCTION_EXPRESSION_HPP_ */
