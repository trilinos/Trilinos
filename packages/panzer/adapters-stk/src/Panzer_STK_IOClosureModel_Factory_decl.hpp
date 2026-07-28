// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_DECL_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory.hpp"

#include "Panzer_STK_Interface.hpp"

#include <vector>
#include <string>

namespace panzer {
  class InputEquationSet;
}

namespace panzer_stk {

  /** \brief A panzer::ClosureModelFactory decorator that provides
    * evaluators for writing fields to exodus files. Wraps a user
    * closure model factory. buildClosureModels() delegates to the
    * wrapped factory for every evaluation type; the
    * panzer::Traits::Residual specialization additionally builds the
    * scatter evaluators needed to write selected nodal/cell fields
    * out to the STK mesh for visualization/output, since output only
    * needs to happen once per workset regardless of evaluation type.
    */
  template<typename EvalT>
  class IOClosureModelFactory : public panzer::ClosureModelFactory<EvalT> {
  public:

    /** \brief Construct with the fields to write specified via an output ParameterList (parsed with parseOutputList()).
      *
      * \param[in] userCMF_ the wrapped closure model factory to delegate buildClosureModels() calls to.
      * \param[in] mesh the STK mesh to write output fields to.
      * \param[in] outputList ParameterList naming the fields to write out.
      */
    IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_,
                          const Teuchos::RCP<STK_Interface> & mesh,
                          const Teuchos::ParameterList & outputList);

    /** \brief Construct with the fields to write specified explicitly, split into nodal and cell fields per element block.
      *
      * \param[in] userCMF_ the wrapped closure model factory to delegate buildClosureModels() calls to.
      * \param[in] mesh the STK mesh to write output fields to.
      * \param[in] nodalFields nodal field names to write out, keyed by element block ID.
      * \param[in] cellFields cell (element-averaged) field names to write out, keyed by element block ID.
      */
    IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_,
                          const Teuchos::RCP<STK_Interface> & mesh,
                          const std::map<std::string,std::vector<std::string> > & nodalFields,
                          const std::map<std::string,std::vector<std::string> > & cellFields);

    /** \brief Delegates to the wrapped factory's buildClosureModels() for the named model.
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

  private:
    void parseOutputList(const Teuchos::ParameterList & pl,
                         std::map<std::string,std::vector<std::string> > & blockIdToFields) const;

    //! Mesh pointer, will be passed around
    Teuchos::RCP<STK_Interface> mesh_;
 
    //! Map showing which cell averaged fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToCellAvgFields_;

    //! Map showing which cell averaged vector fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToCellAvgVectors_;

    //! Map showing which cell fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToCellFields_;

    //! Map showing which nodal fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToNodalFields_;

    /** Map stating if an evaluator for a particular block ID has been included.
      *
      * This is a bit of hack that is done to gurantee only one evaluator is
      * added to each field manager. However, if an instantiation of this closure model factory is
      * used in multiple places then the appropriate evaluator will be added and
      * required only once. So its likely that not every field manager will have (and require)
      * the scatter cell evaluators.
      */    
    mutable std::map<std::string,bool> blockIdEvaluated_;

    //! we will reuse the drekar closure model factory
    Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_;
  };

  /** \brief Delegates to the wrapped factory's buildClosureModels() for the named model, then adds the scatter-to-mesh evaluators needed for the fields requested at construction.
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
  template < >
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
  panzer_stk::IOClosureModelFactory<panzer::Traits::Residual>::buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

}

#endif
