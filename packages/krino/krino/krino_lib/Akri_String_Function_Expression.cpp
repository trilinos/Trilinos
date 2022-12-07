// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_String_Function_Expression.hpp>
#include <stk_expreval/Eval.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

String_Function_Expression::String_Function_Expression(const std::string & expression)
: myEvaluator(*this)
{
  parse(expression);
}

void String_Function_Expression::parse(const std::string & expression)
{
  try {
    myEvaluator.parse(expression);
  }
  catch (std::runtime_error &x) {
    stk::RuntimeDoomedSymmetric() << "In expression '" << expression << "':" << std::endl << x.what() << std::endl;
  }
}

void String_Function_Expression::resolve(stk::expreval::VariableMap::iterator & varIt)
{
  std::string name = (*varIt).first;

  if (!(name).compare("x"))
    (*varIt).second->bind(myQueryCoords[0]);
  else if (!(name).compare("y"))
    (*varIt).second->bind(myQueryCoords[1]);
  else if (!(name).compare("z"))
    (*varIt).second->bind(myQueryCoords[2]);
  else
  {
    std::ostringstream msg;
    msg << "  Unable to resolve symbol: " << name;
    throw std::runtime_error(msg.str());
  }
}

double
String_Function_Expression::evaluate(const Vector3d &coord) const
{
  myQueryCoords = coord;
  return myEvaluator.evaluate();
}

}
