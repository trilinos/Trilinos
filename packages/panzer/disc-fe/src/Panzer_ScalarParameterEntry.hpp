// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCALAR_PARAMETER_ENTRY_HPP
#define PANZER_SCALAR_PARAMETER_ENTRY_HPP

#include "Panzer_Traits.hpp"
#include "Sacado_ScalarParameterEntry.hpp"
#include "Panzer_EvaluationTraits.hpp"

namespace panzer {

  /** Holds a single scalar parameter's value for a given evaluation
    * type, as tracked by the Sacado parameter library. Used to expose
    * a named, possibly AD-differentiable, scalar (e.g. a physical
    * constant) that can be looked up and updated through Panzer's
    * scalar parameter registry.
    *
    * \tparam EvalType the evaluation type (e.g. Residual, Jacobian) whose scalar type this entry stores.
    */
  template <typename EvalType>
  class ScalarParameterEntry : public  Sacado::ScalarParameterEntry<EvalType,panzer::EvaluationTraits> {

  public:

    typedef typename Sacado::ScalarParameterEntry<EvalType,panzer::EvaluationTraits>::ScalarT ScalarT;

    /// \brief Sets the parameter's value from a plain double, converting to the evaluation type's scalar type.
    void setRealValue(double value)
    {
      m_parameter = ScalarT(value);
    }

    /// \brief Sets the parameter's value directly, in the evaluation type's scalar type.
    void setValue(const ScalarT& value)
    {
      m_parameter = value;
    }

    /// \brief Returns the parameter's current value.
    const ScalarT& getValue() const
    {
      return m_parameter;
    }

  private:
    
    ScalarT m_parameter;

  };

}

#endif
