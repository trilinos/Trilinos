// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_UTILITIES_H
#define PHX_EVALUATOR_UTILITIES_H

#include <vector>

#include "Phalanx_FieldManager.hpp"

namespace PHX {

  // Forward declarations
  template<typename DataT, int Rank, typename Layout> class Field;
  template <typename DataT,typename...Props> class MDField;

  /*! @brief Utilities to hide templating in concrete Evaluators. */
  template<typename EvalT, typename Traits>
  struct EvaluatorUtilities {

    template <typename DataT,typename...Props>
    void setFieldData(PHX::MDField<DataT,Props...>& f,
                      PHX::FieldManager<Traits>& fm)
    {
      fm.template getFieldData<EvalT>(f);
    }

    template <typename DataT,int Rank>
    void setFieldData(PHX::Field<DataT,Rank>& f,
                      PHX::FieldManager<Traits>& fm)
    {
      fm.template getFieldData<EvalT,DataT>(f);
    }

    template<typename DataT,typename... Props>
    void setFieldData(const PHX::FieldTag& /* ft */,
                      Kokkos::View<DataT,Props...>& /* f */,
                      PHX::FieldManager<Traits>& /* fm */)
    {}

  };
}

#endif
