#ifndef PANZER_SCALAR_PARAMETER_ENTRY_HPP
#define PANZER_SCALAR_PARAMETER_ENTRY_HPP

#include "Panzer_Traits.hpp"
#include "Sacado_ScalarParameterEntry.hpp"

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
