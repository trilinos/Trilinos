
#ifndef PHX_EVALUATOR_UTILITIES_H
#define PHX_EVALUATOR_UTILITIES_H

#include <vector>

#include "Phalanx_Field.hpp"
#include "Phalanx_FieldManager.hpp"

namespace PHX {

  /*! @brief Utilities to hide templating in concrete Evaluators.
   
  */
  template<typename EvalT, typename Traits> 
  struct EvaluatorUtilities {
    
    template <typename DataT>
    void setFieldData(PHX::Field<DataT>& f, PHX::FieldManager<Traits>& fm) 
    {
      fm.template getFieldData<DataT,EvalT>(f);
    }

  };
}

#endif
