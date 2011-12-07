#ifndef PANZER_PARAMETER_LIBRARY_HPP
#define PANZER_PARAMETER_LIBRARY_HPP

#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterVector.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_EvaluationTraits.hpp"

namespace panzer {
  
  typedef Sacado::ScalarParameterLibrary<panzer::EvaluationTraits> ParamLib;
  typedef Sacado::ScalarParameterVector<panzer::EvaluationTraits> ParamVec;

}

#endif
