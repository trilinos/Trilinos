// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __MiniEM_ClosureModelFactory_hpp__
#define __MiniEM_ClosureModelFactory_hpp__

#include "Panzer_ClosureModel_Factory.hpp"

namespace panzer {
  class InputEquationSet;
}

namespace mini_em {

/** \brief mini-em's application-specific panzer::ClosureModelFactory,
  * dispatching on each requested model's "Type" parameter to build the
  * source terms, analytic solutions, and conductivity models used by
  * the mini-em Darcy/Maxwell equation sets (e.g. GaussianPulse,
  * RandomForcing, DarcyAnalyticForcing/Solution,
  * MaxwellAnalyticForcing/Solution, TensorConductivity,
  * VariableTensorConductivity, PiecewiseConstant, RTC), in addition to
  * generic models like Constant, DotProduct, Product, and Sum.
  */
template<typename EvalT>
class ModelFactory : public panzer::ClosureModelFactory<EvalT> {

public:

   /** \brief Builds the evaluators for the named model, dispatching on its "Type" parameter as documented in the class description.
     *
     * \param[in] model_id the closure model ID to build; selects which sublist of models to use.
     * \param[in] models ParameterList of named closure model sublists, keyed by model ID.
     * \param[in] fl the field layouts available for the current physics block.
     * \param[in] ir the integration rule the built evaluators should use.
     * \param[in] default_params equation-set-level parameters used as defaults when building closures.
     * \param[in] user_data user-supplied parameters passed through to closures that need additional application data.
     * \param[in] global_data global data (parameter library, output stream) shared across the problem.
     * \param[in] fm the field manager the built evaluators will be registered into.
     */
   Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
   buildClosureModels(const std::string& model_id,
                      const Teuchos::ParameterList& models,
		      const panzer::FieldLayoutLibrary& fl,
		      const Teuchos::RCP<panzer::IntegrationRule>& ir,
                      const Teuchos::ParameterList& default_params,
                      const Teuchos::ParameterList& user_data,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      PHX::FieldManager<panzer::Traits>& fm) const;

};

}

#include "MiniEM_ClosureModel_Factory_impl.hpp"

#endif
