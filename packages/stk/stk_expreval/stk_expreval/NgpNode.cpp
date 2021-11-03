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

#include "stk_expreval/NgpNode.hpp"
#include "stk_expreval/Function.hpp"

namespace stk {
namespace expreval {

KOKKOS_FUNCTION
NgpNode::NgpNode()
  : m_opcode(OPCODE_UNDEFINED),
    m_data{{0.0}},
    m_currentNodeIndex(-1),
    m_nextNodeIndex(-1),
    m_ternaryTrueNextNodeIndex(-1),
    m_ternaryFalseNextNodeIndex(-1),
    m_leftNodeIndex(-1),
    m_rightNodeIndex(-1),
    m_ternaryOtherNodeIndex(-1)
{}

NgpNode::NgpNode(const Node& node)
  : m_opcode(node.m_opcode),
    m_resultIdx(node.m_resultIdx),
    m_currentNodeIndex(node.m_currentNodeIndex),
    m_nextNodeIndex(node.m_nextNodeIndex),
    m_ternaryTrueNextNodeIndex(node.m_ternaryTrueNextNodeIndex),
    m_ternaryFalseNextNodeIndex(node.m_ternaryFalseNextNodeIndex),
    m_leftNodeIndex((node.m_left != nullptr) ? node.m_left->m_currentNodeIndex : -1),
    m_rightNodeIndex((node.m_right != nullptr) ? node.m_right->m_currentNodeIndex : -1),
    m_ternaryOtherNodeIndex((node.m_ternaryOther != nullptr) ? node.m_ternaryOther->m_currentNodeIndex : -1)
{
  if (m_opcode == OPCODE_CONSTANT) {
    m_data.constant.value = node.m_data.constant.value;
  }
  else if (m_opcode == OPCODE_RVALUE) {
    m_data.variable.variableIndex = node.m_data.variable.variable->get_index();
  }
  else if (m_opcode == OPCODE_ASSIGN) {
    m_data.variable.variableIndex = node.m_data.variable.variable->get_index();
    m_data.variable.variableType = node.m_data.variable.variable->get_type();
    m_data.variable.variableSize = node.m_data.variable.variable->getLength();
  }
  else if (m_opcode == OPCODE_FUNCTION) {
    m_data.function.functionType = node.m_data.function.functionType;
  }
}

KOKKOS_FUNCTION
double
NgpNode::evaluate_function(int argumentCount, double * arguments) const
{
  switch (m_data.function.functionType) {
  case FunctionType::ABS : {
    if (argumentCount == 1) {
      return std::fabs(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for abs function");
    break;
  }
  case FunctionType::MAX : {
    if (argumentCount == 2) {
      return max_2(arguments[0], arguments[1]);
    }
    else if (argumentCount == 3) {
      return max_3(arguments[0], arguments[1], arguments[2]);
    }
    else if (argumentCount == 4) {
      return max_4(arguments[0], arguments[1], arguments[2], arguments[3]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for max function");
    break;
  }
  case FunctionType::MIN : {
    if (argumentCount == 2) {
      return min_2(arguments[0], arguments[1]);
    }
    else if (argumentCount == 3) {
      return min_3(arguments[0], arguments[1], arguments[2]);
    }
    else if (argumentCount == 4) {
      return min_4(arguments[0], arguments[1], arguments[2], arguments[3]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for min function");
    break;
  }
  case FunctionType::SIGN : {
    if (argumentCount == 1) {
      return sign(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for sign function");
    break;
  }
  case FunctionType::IPART : {
    if (argumentCount == 1) {
      return ipart(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for ipart function");
    break;
  }
  case FunctionType::FPART : {
    if (argumentCount == 1) {
      return fpart(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for fpart function");
    break;
  }
  case FunctionType::CEIL : {
    if (argumentCount == 1) {
      return std::ceil(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for ceil function");
    break;
  }
  case FunctionType::FLOOR : {
    if (argumentCount == 1) {
      return std::floor(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for floor function");
    break;
  }
  case FunctionType::MOD : {
    if (argumentCount == 2) {
      return std::fmod(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for fmod or mod function");
    break;
  }
  case FunctionType::POW : {
    if (argumentCount == 2) {
      return std::pow(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for pow function");
    break;
  }
  case FunctionType::SQRT : {
    if (argumentCount == 1) {
      return std::sqrt(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for sqrt function");
    break;
  }
  case FunctionType::EXP : {
    if (argumentCount == 1) {
      return std::exp(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for exp function");
    break;
  }
  case FunctionType::LN : {
    if (argumentCount == 1) {
      return std::log(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for ln or log function");
    break;
  }
  case FunctionType::LOG10 : {
    if (argumentCount == 1) {
      return std::log10(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for log10 function");
    break;
  }
  case FunctionType::DEG : {
    if (argumentCount == 1) {
      return radian_to_degree()*arguments[0];
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for deg function");
    break;
  }
  case FunctionType::RAD : {
    if (argumentCount == 1) {
      return degree_to_radian()*arguments[0];
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for rad function");
    break;
  }
  case FunctionType::SIN : {
    if (argumentCount == 1) {
      return std::sin(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for sin function");
    break;
  }
  case FunctionType::COS : {
    if (argumentCount == 1) {
      return std::cos(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for cos function");
    break;
  }
  case FunctionType::TAN : {
    if (argumentCount == 1) {
      return std::tan(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for tan function");
    break;
  }
  case FunctionType::ASIN : {
    if (argumentCount == 1) {
      return std::asin(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for asin function");
    break;
  }
  case FunctionType::ACOS : {
    if (argumentCount == 1) {
      return std::acos(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for acos function");
    break;
  }
  case FunctionType::ATAN : {
    if (argumentCount == 1) {
      return std::atan(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for atan function");
    break;
  }
  case FunctionType::ATAN2 : {
    if (argumentCount == 2) {
      return std::atan2(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for atan2 function");
    break;
  }
  case FunctionType::SINH : {
    if (argumentCount == 1) {
      return std::sinh(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for sinh function");
    break;
  }
  case FunctionType::COSH : {
    if (argumentCount == 1) {
      return std::cosh(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arugments for cosh function");
    break;
  }
  case FunctionType::TANH : {
    if (argumentCount == 1) {
      return std::tanh(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arugments for tanh function");
    break;
  }
  case FunctionType::ASINH : {
    if (argumentCount == 1) {
      return std::asinh(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for asinh function");
    break;
  }
  case FunctionType::ACOSH : {
    if (argumentCount == 1) {
      return std::acosh(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for acosh function");
    break;
  }
  case FunctionType::ATANH : {
    if (argumentCount == 1) {
      return std::atanh(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for atanh function");
    break;
  }
  case FunctionType::ERF : {
    if (argumentCount == 1) {
      return std::erf(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for erf function");
    break;
  }
  case FunctionType::ERFC : {
    if (argumentCount == 1) {
      return std::erfc(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for erfc function");
    break;
  }
  case FunctionType::POLTORECTX : {
    if (argumentCount == 2) {
      return poltorectx(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for poltorectx function");
    break;
  }
  case FunctionType::POLTORECTY : {
    if (argumentCount == 2) {
      return poltorecty(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for poltorecty function");
    break;
  }
  case FunctionType::RECTTOPOLR : {
    if (argumentCount == 2) {
      return recttopolr(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for recttopolr function");
    break;
  }
  case FunctionType::RECTTOPOLA : {
    if (argumentCount == 2) {
      return recttopola(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for recttopola function");
    break;
  }
  case FunctionType::UNIT_STEP : {
    if (argumentCount == 3) {
      return unit_step3(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for unit_step function");
    break;
  }
  case FunctionType::CYCLOIDAL_RAMP : {
    if (argumentCount == 3) {
      return cycloidal_ramp(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for cycloidal_ramp function");
    break;
  }
  case FunctionType::COS_RAMP : {
    if (argumentCount == 1) {
      return cosine_ramp1(arguments[0]);
    }
    else if (argumentCount == 2) {
      return cosine_ramp2(arguments[0], arguments[1]);
    }
    else if (argumentCount == 3) {
      return cosine_ramp3(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for cos_ramp or cosine_ramp function");
    break;
  }
  case FunctionType::HAVERSINE_PULSE : {
    if (argumentCount == 3) {
      return haversine_pulse(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for haversine_pulse function");
    break;
  }
  case FunctionType::POINT2D : {
    if (argumentCount == 4) {
      return point_2(arguments[0], arguments[1], arguments[2], arguments[3]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for pulse_2 function");
    break;
  }
  case FunctionType::POINT3D : {
    if (argumentCount == 5) {
      return point_3(arguments[0], arguments[1], arguments[2], arguments[3], arguments[4]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for pulse_3 function");
    break;
  }
  case FunctionType::EXPONENTIAL_PDF : {
    if (argumentCount == 2) {
      return exponential_pdf(arguments[0], arguments[1]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for exponential_pdf function");
    break;
  }
  case FunctionType::LOG_UNIFORM_PDF : {
    if (argumentCount == 3) {
      return log_uniform_pdf(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for log_uniform_pdf function");
    break;
  }
  case FunctionType::NORMAL_PDF : {
    if (argumentCount == 3) {
      return normal_pdf(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for normal_pdf function");
    break;
  }
  case FunctionType::WEIBULL_PDF : {
    if (argumentCount == 3) {
      return weibull_pdf(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for weibull_pdf function");
    break;
  }
  case FunctionType::GAMMA_PDF : {
    if (argumentCount == 3) {
      return gamma_pdf(arguments[0], arguments[1], arguments[2]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for gamma_pdf function");
    break;
  }
  case FunctionType::RAND : {
    if (argumentCount == 0) {
      return real_rand();
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for rand function");
    break;
  }
  case FunctionType::SRAND : {
    if (argumentCount == 1) {
      return real_srand(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for srand function");
    break;
  }
  case FunctionType::RANDOM : {
    if (argumentCount == 0) {
      return random0();
    }
    else if (argumentCount == 1) {
      return random1(arguments[0]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for random function");
    break;
  }
  case FunctionType::TS_RANDOM : {
    if (argumentCount == 4) {
      return time_space_random(arguments[0], arguments[1], arguments[2], arguments[3]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for ts_random function");
    break;
  }
  case FunctionType::TS_NORMAL : {
    if (argumentCount == 8) {
      return time_space_normal(arguments[0], arguments[1], arguments[2], arguments[3], arguments[4],
          arguments[5], arguments[6], arguments[7]);
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for ts_normal function");
    break;
  }
  case FunctionType::TIME : {
    if (argumentCount == 0) {
      return current_time();
    }
    NGP_ThrowErrorMsg("Incorrect number of arguments for time function");
    break;
  }
  case FunctionType::UNDEFINED : {
    NGP_ThrowErrorMsg("Undefined function type");
    break;
  }
  default : {
    break;
  }
  }

  return 0.0;
}

}
}
