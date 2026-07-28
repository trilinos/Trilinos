// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_DECL_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_DECL_HPP

#include "Panzer_ClosureModel_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"

namespace panzer {

  /** \brief A panzer::ClosureModelFactory that delegates to a list of
    * other user defined, internally owned closure model factories.
    *
    * buildClosureModels() asks each of its constituent factories, in
    * order, to build the requested model, and collects the evaluators
    * from every factory that produces any; this allows closure models
    * from independently-registered factories to be combined for a
    * single model ID.
    */
  template<typename EvalT>
  class ClosureModelFactoryComposite : public panzer::ClosureModelFactory<EvalT> {

  public:

    /// \brief Construct from the list of factories to delegate to, in the order they will be tried.
    ClosureModelFactoryComposite(const std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >& factories);

    /** \brief Builds the evaluators for the named model by delegating to each constituent factory and combining their results.
      *
      * \param[in] model_id the closure model ID to build; selects which sublist of models to use.
      * \param[in] models ParameterList of named closure model sublists, keyed by model ID.
      * \param[in] fl the field layouts available for the current physics block.
      * \param[in] ir the integration rule the built evaluators should evaluate on.
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

    /// \brief Returns the list of constituent factories this composite delegates to.
    std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits>>>&
    getFactories(){return m_factories;}

  private:

    std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits>>> m_factories;

  };

}

#endif
