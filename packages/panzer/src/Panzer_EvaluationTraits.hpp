#ifndef PANZER_EVALUATION_TRAITS_HPP
#define PANZER_EVALUATION_TRAITS_HPP

namespace panzer {
  
  struct EvaluationTraits {
    template <class T> struct apply {
      typedef typename T::ScalarT type;
    };
  };

}

#endif
