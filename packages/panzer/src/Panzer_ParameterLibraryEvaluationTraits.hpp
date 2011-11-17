#ifndef PANZER_SACADO_PARAM_LIBRARY_TRAITS_HPP
#define PANZER_SACADO_PARAM_LIBRARY_TRAITS_HPP

#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterVector.hpp"
#include "Panzer_Traits.hpp"

namespace panzer {
  
  struct EvaluationTraits {
    template <class T> struct apply {
      typedef typename T::ScalarT type;
    };
  };

  typedef Sacado::ScalarParameterLibrary<panzer::EvaluationTraits> ParamLib;
  typedef Sacado::ScalarParameterVector<panzer::EvaluationTraits> ParamVec;

}

#endif
