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

#ifndef NGPNODE_HPP
#define NGPNODE_HPP

#include "Kokkos_Core.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_expreval/Variable.hpp"
#include "stk_expreval/Node.hpp"
#include "stk_expreval/Function.hpp"

namespace stk {
namespace expreval {

static constexpr int DEFAULT_MAX_BOUND_VARIABLES = 16;

class Node;

template <int MAX_BOUND_VARIABLES=DEFAULT_MAX_BOUND_VARIABLES>
class DeviceVariableMap;

class NgpNode
{
public:
  KOKKOS_FUNCTION
  NgpNode()
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

  explicit NgpNode(const Node& node);

  KOKKOS_DEFAULTED_FUNCTION
  NgpNode(const NgpNode &) = default;

  KOKKOS_DEFAULTED_FUNCTION
  NgpNode &operator=(const NgpNode &) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~NgpNode() = default;

  template <int MAX_BOUND_VARIABLES>
  KOKKOS_FUNCTION
  double
  eval(DeviceVariableMap<MAX_BOUND_VARIABLES> & deviceVariableMap, double resultBuffer[])
  {
    switch (m_opcode) {
    case OPCODE_STATEMENT: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_CONSTANT: {
      setResult(resultBuffer) = m_data.constant.value;
      break;
    }
    case OPCODE_RVALUE: {
      if (get_left_node()) {
        setResult(resultBuffer) = deviceVariableMap[m_data.variable.variableIndex].getArrayValue(get_left_node()->getResult(resultBuffer), deviceVariableMap.get_array_offset_type());
      }
      else {
        setResult(resultBuffer) = deviceVariableMap[m_data.variable.variableIndex].getValue();
      }
      break;
    }
    case OPCODE_MULTIPLY: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer)*get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_EXPONENIATION: {
      setResult(resultBuffer) = std::pow(get_left_node()->getResult(resultBuffer),get_right_node()->getResult(resultBuffer));
      break;
    }
    case OPCODE_DIVIDE: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer)/get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_MODULUS: {
      setResult(resultBuffer) = std::fmod(get_left_node()->getResult(resultBuffer), get_right_node()->getResult(resultBuffer));
      break;
    }
    case OPCODE_ADD: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) + get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_SUBTRACT: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) - get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) == get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_NOT_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) != get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LESS: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) < get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_GREATER: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) > get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LESS_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) <= get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_GREATER_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) >= get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LOGICAL_AND: {
      double left = get_left_node()->getResult(resultBuffer);
      double right = get_right_node()->getResult(resultBuffer);
      setResult(resultBuffer) = (left != 0.0) && (right != 0.0) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LOGICAL_OR: {
      double left = get_left_node()->getResult(resultBuffer);
      double right = get_right_node()->getResult(resultBuffer);
      setResult(resultBuffer) = (left != 0.0) || (right != 0.0) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_TERNARY_PREDICATE: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_TERNARY_JOIN: {
      if (get_ternary_other_node()->getResult(resultBuffer) == 0.0) {
        setResult(resultBuffer) = get_right_node()->getResult(resultBuffer);
      }
      else {
        setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      }
      break;
    }
    case OPCODE_UNARY_MINUS: {
      setResult(resultBuffer) = -get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_UNARY_NOT: {
      setResult(resultBuffer) = get_right_node()->getResult(resultBuffer) == 0.0 ? 1.0 : 0.0;
      break;
    }
    case OPCODE_ASSIGN: {
      if (get_left_node()) {
        deviceVariableMap[m_data.variable.variableIndex].assignArrayValue(get_left_node()->getResult(resultBuffer),
                                                                          deviceVariableMap.get_array_offset_type(),
                                                                          get_right_node()->getResult(resultBuffer));
      }
      else {
        deviceVariableMap[m_data.variable.variableIndex] = get_right_node()->getResult(resultBuffer);
      }
      setResult(resultBuffer) = get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_ARGUMENT: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_FUNCTION: {
      double arguments[20];

      int argumentCount = 0;
      for (const NgpNode* arg = get_right_node(); arg; arg = arg->get_right_node()) {
        arguments[argumentCount++] = arg->getResult(resultBuffer);
      }

      setResult(resultBuffer) = evaluate_function(argumentCount, arguments);
      break;
    }
    default: { //  Unknown opcode
      setResult(resultBuffer) = 0.0;
    }
    }

    return 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  int getNextNodeIndex(double resultBuffer[]) {
    if (m_opcode == OPCODE_TERNARY_PREDICATE) {
      if (setResult(resultBuffer) == 0.0) {
        return m_ternaryFalseNextNodeIndex;
      }
      else {
        return m_ternaryTrueNextNodeIndex;
      }
    }
    return m_nextNodeIndex;
  }

  KOKKOS_INLINE_FUNCTION
  double getResult(double resultBuffer[]) const { return resultBuffer[m_resultIdx]; }

  KOKKOS_INLINE_FUNCTION
  double& setResult(double resultBuffer[]) const { return resultBuffer[m_resultIdx]; }

  KOKKOS_INLINE_FUNCTION
  const NgpNode* get_node(int nodeIndex) const {
    return nodeIndex != -1 ? (this + (nodeIndex - m_currentNodeIndex)) : nullptr;
  }

  KOKKOS_INLINE_FUNCTION const NgpNode* get_left_node() const { return get_node(m_leftNodeIndex); }
  KOKKOS_INLINE_FUNCTION const NgpNode* get_right_node() const { return get_node(m_rightNodeIndex); }
  KOKKOS_INLINE_FUNCTION const NgpNode* get_ternary_other_node() const { return get_node(m_ternaryOtherNodeIndex); }

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

  int m_resultIdx;
  int m_currentNodeIndex;
  int m_nextNodeIndex;
  int m_ternaryTrueNextNodeIndex;
  int m_ternaryFalseNextNodeIndex;
  int m_leftNodeIndex;
  int m_rightNodeIndex;
  int m_ternaryOtherNodeIndex;

  KOKKOS_FUNCTION
  double evaluate_function(int argumentCount, double* arguments) const
  {
    switch (m_data.function.functionType) {
    case FunctionType::ABS : {
      if (argumentCount == 1) {
        return std::fabs(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for abs function");
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
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for max function");
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
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for min function");
      break;
    }
    case FunctionType::SIGN : {
      if (argumentCount == 1) {
        return sign(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for sign function");
      break;
    }
    case FunctionType::IPART : {
      if (argumentCount == 1) {
        return ipart(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for ipart function");
      break;
    }
    case FunctionType::FPART : {
      if (argumentCount == 1) {
        return fpart(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for fpart function");
      break;
    }
    case FunctionType::CEIL : {
      if (argumentCount == 1) {
        return std::ceil(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for ceil function");
      break;
    }
    case FunctionType::FLOOR : {
      if (argumentCount == 1) {
        return std::floor(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for floor function");
      break;
    }
    case FunctionType::MOD : {
      if (argumentCount == 2) {
        return std::fmod(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for fmod or mod function");
      break;
    }
    case FunctionType::POW : {
      if (argumentCount == 2) {
        return std::pow(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for pow function");
      break;
    }
    case FunctionType::SQRT : {
      if (argumentCount == 1) {
        return std::sqrt(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for sqrt function");
      break;
    }
    case FunctionType::EXP : {
      if (argumentCount == 1) {
        return std::exp(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for exp function");
      break;
    }
    case FunctionType::LN : {
      if (argumentCount == 1) {
        return std::log(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for ln or log function");
      break;
    }
    case FunctionType::LOG10 : {
      if (argumentCount == 1) {
        return std::log10(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for log10 function");
      break;
    }
    case FunctionType::DEG : {
      if (argumentCount == 1) {
        return radian_to_degree()*arguments[0];
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for deg function");
      break;
    }
    case FunctionType::RAD : {
      if (argumentCount == 1) {
        return degree_to_radian()*arguments[0];
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for rad function");
      break;
    }
    case FunctionType::SIN : {
      if (argumentCount == 1) {
        return std::sin(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for sin function");
      break;
    }
    case FunctionType::COS : {
      if (argumentCount == 1) {
        return std::cos(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for cos function");
      break;
    }
    case FunctionType::TAN : {
      if (argumentCount == 1) {
        return std::tan(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for tan function");
      break;
    }
    case FunctionType::ASIN : {
      if (argumentCount == 1) {
        return std::asin(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for asin function");
      break;
    }
    case FunctionType::ACOS : {
      if (argumentCount == 1) {
        return std::acos(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for acos function");
      break;
    }
    case FunctionType::ATAN : {
      if (argumentCount == 1) {
        return std::atan(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for atan function");
      break;
    }
    case FunctionType::ATAN2 : {
      if (argumentCount == 2) {
        return std::atan2(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for atan2 function");
      break;
    }
    case FunctionType::SINH : {
      if (argumentCount == 1) {
        return std::sinh(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for sinh function");
      break;
    }
    case FunctionType::COSH : {
      if (argumentCount == 1) {
        return std::cosh(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arugments for cosh function");
      break;
    }
    case FunctionType::TANH : {
      if (argumentCount == 1) {
        return std::tanh(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arugments for tanh function");
      break;
    }
    case FunctionType::ASINH : {
      if (argumentCount == 1) {
        return std::asinh(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for asinh function");
      break;
    }
    case FunctionType::ACOSH : {
      if (argumentCount == 1) {
        return std::acosh(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for acosh function");
      break;
    }
    case FunctionType::ATANH : {
      if (argumentCount == 1) {
        return std::atanh(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for atanh function");
      break;
    }
    case FunctionType::ERF : {
      if (argumentCount == 1) {
        return std::erf(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for erf function");
      break;
    }
    case FunctionType::ERFC : {
      if (argumentCount == 1) {
        return std::erfc(arguments[0]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for erfc function");
      break;
    }
    case FunctionType::POLTORECTX : {
      if (argumentCount == 2) {
        return poltorectx(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for poltorectx function");
      break;
    }
    case FunctionType::POLTORECTY : {
      if (argumentCount == 2) {
        return poltorecty(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for poltorecty function");
      break;
    }
    case FunctionType::RECTTOPOLR : {
      if (argumentCount == 2) {
        return recttopolr(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for recttopolr function");
      break;
    }
    case FunctionType::RECTTOPOLA : {
      if (argumentCount == 2) {
        return recttopola(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for recttopola function");
      break;
    }
    case FunctionType::UNIT_STEP : {
      if (argumentCount == 3) {
        return unit_step3(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for unit_step function");
      break;
    }
    case FunctionType::CYCLOIDAL_RAMP : {
      if (argumentCount == 3) {
        return cycloidal_ramp(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for cycloidal_ramp function");
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
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for cos_ramp or cosine_ramp function");
      break;
    }
    case FunctionType::LINEAR_RAMP : {
      if (argumentCount == 3) {
        return linear_ramp3(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for linear_ramp function");
      break;
    }
    case FunctionType::HAVERSINE_PULSE : {
      if (argumentCount == 3) {
        return haversine_pulse(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for haversine_pulse function");
      break;
    }
    case FunctionType::POINT2D : {
      if (argumentCount == 4) {
        return point_2(arguments[0], arguments[1], arguments[2], arguments[3]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for pulse_2 function");
      break;
    }
    case FunctionType::POINT3D : {
      if (argumentCount == 5) {
        return point_3(arguments[0], arguments[1], arguments[2], arguments[3], arguments[4]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for pulse_3 function");
      break;
    }
    case FunctionType::EXPONENTIAL_PDF : {
      if (argumentCount == 2) {
        return exponential_pdf(arguments[0], arguments[1]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for exponential_pdf function");
      break;
    }
    case FunctionType::LOG_UNIFORM_PDF : {
      if (argumentCount == 3) {
        return log_uniform_pdf(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for log_uniform_pdf function");
      break;
    }
    case FunctionType::NORMAL_PDF : {
      if (argumentCount == 3) {
        return normal_pdf(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for normal_pdf function");
      break;
    }
    case FunctionType::WEIBULL_PDF : {
      if (argumentCount == 3) {
        return weibull_pdf(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for weibull_pdf function");
      break;
    }
    case FunctionType::GAMMA_PDF : {
      if (argumentCount == 3) {
        return gamma_pdf(arguments[0], arguments[1], arguments[2]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for gamma_pdf function");
      break;
    }
    case FunctionType::TS_RANDOM : {
      if (argumentCount == 4) {
        return time_space_random(arguments[0], arguments[1], arguments[2], arguments[3]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for ts_random function");
      break;
    }
    case FunctionType::TS_NORMAL : {
      if (argumentCount == 8) {
        return time_space_normal(arguments[0], arguments[1], arguments[2], arguments[3], arguments[4],
            arguments[5], arguments[6], arguments[7]);
      }
      STK_NGP_ThrowErrorMsg("Incorrect number of arguments for ts_normal function");
      break;
    }
    case FunctionType::UNDEFINED : {
      STK_NGP_ThrowErrorMsg("Undefined function type");
      break;
    }
    default : {
      break;
    }
    }

    return 0.0;
  }
};

}
}

#endif // NGPNODE_HPP
