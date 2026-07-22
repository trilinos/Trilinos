// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_MULTIVARIATE_PARAMETER_DECL_HPP
#define PANZER_EVALUATOR_MULTIVARIATE_PARAMETER_DECL_HPP

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_ParameterLibrary.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

  template <typename EvalT> class ScalarParameterEntry;

  //! Constant parameter from sacado parameter library
  /*!
   * Represents a sum of n independent parameters.  Useful for SA/UQ scaling
   * studies where the number of parameters must be varied.
   */
  template<typename EvalT, typename TRAITS>
  class MultiVariateParameter :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS> {

  public:

    MultiVariateParameter(const std::string parameter_name,
                          const int num_param,
                          const std::string field_name,
                          const Teuchos::RCP<PHX::DataLayout>& data_layout,
                          panzer::ParamLib& param_lib);

    void evaluateFields(typename TRAITS::EvalData ud);

  private:

    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT, Cell, Point> target_field;

    std::size_t cell_data_size;

    Teuchos::Array< Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > > param;
  };

}

#endif
