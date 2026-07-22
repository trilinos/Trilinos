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

  template <typename EvalType>
  class ScalarParameterEntry : public  Sacado::ScalarParameterEntry<EvalType,panzer::EvaluationTraits> {
    
  public:

    typedef typename Sacado::ScalarParameterEntry<EvalType,panzer::EvaluationTraits>::ScalarT ScalarT;
    
    void setRealValue(double value)
    {
      m_parameter = ScalarT(value);
    }
    
    void setValue(const ScalarT& value)
    {
      m_parameter = value;
    }

    const ScalarT& getValue() const
    {
      return m_parameter;
    }

  private:
    
    ScalarT m_parameter;

  };

}

#endif
