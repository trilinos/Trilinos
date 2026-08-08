// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseFactory_BCStrategyAdapter_hpp__
#define __Panzer_ResponseFactory_BCStrategyAdapter_hpp__

#include <vector>
#include <string>
#include <tuple>

#include "Teuchos_RCP.hpp"

#include "Panzer_BCStrategy.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Normals.hpp"

#include "Panzer_Traits.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_MDField.hpp"

/**
 * \file Panzer_ResponseFactory_BCStrategyAdapter.hpp
 * \brief Adapts ResponseEvaluatorFactory objects so they can be driven through the BCStrategy interface.
 *
 * Used only by the ResponseLibrary, to unify the interface between
 * volumetric response objects and surface (boundary condition side
 * set) response objects: surface responses are registered through the
 * BCStrategy machinery, so this file wraps a set of
 * ResponseEvaluatorFactory objects to look like a BCStrategy.
 */

namespace panzer {
namespace response_bc_adapters {

  /**
   * \brief Adapts a set of named ResponseEvaluatorFactory objects into a single BCStrategy, so responses can be evaluated on a boundary condition side set.
   * \tparam EvalT the evaluation type this adapter builds evaluators for.
   */
  template <typename EvalT>
  class ResponseFactory_BCStrategyAdapter : public panzer::BCStrategy<EvalT>
  {
    /// (response name, factory) pairs to register evaluators for.
    std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > refVec_;

  public:

    /// \brief Construct from the boundary condition being adapted and the responses to build on it.
    ResponseFactory_BCStrategyAdapter(const panzer::BC & bc,const std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > & refVec)
      : panzer::BCStrategy<EvalT>(bc), refVec_(refVec) {}

    virtual ~ResponseFactory_BCStrategyAdapter() {}

    //! \name Derived from BCStrategy
    //@{

    /// \brief No-op: this adapter requires no additional setup.
    virtual void setup(const panzer::PhysicsBlock& /* side_pb */, const Teuchos::ParameterList& /* user_data */) {}

    /** \brief Registers the side physics block's equation-set/closure-model evaluators, then delegates to each wrapped response factory (for evaluation types it supports) to register its response evaluators.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] side_pb the physics block for the boundary condition's side set.
      * \param[in] factory closure model factory used to build the side physics block's closure model evaluators.
      * \param[in] models the closure model ParameterList to build from.
      * \param[in] user_data user-supplied parameters passed through to the equation-set/closure-model evaluators.
      */
    virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					    const panzer::PhysicsBlock& side_pb,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					    const Teuchos::ParameterList& models,
					    const Teuchos::ParameterList& user_data) const
    {
      side_pb.buildAndRegisterEquationSetEvaluators(fm, user_data);
      side_pb.buildAndRegisterClosureModelEvaluatorsForType<EvalT>(fm,factory,models,user_data);

      for(std::size_t i=0;i<refVec_.size();i++) {
        Teuchos::RCP<const ResponseEvaluatorFactoryBase> respEvalFact = refVec_[i].second->template getAsBase<EvalT>();

        // only register evaluators if the type is supported
        if(respEvalFact!=Teuchos::null && respEvalFact->typeSupported())
          respEvalFact->buildAndRegisterEvaluators(refVec_[i].first,fm,side_pb,user_data);
      }
    }

    /// \brief No-op: response adapters do not scatter residual/Jacobian contributions.
    virtual void
    buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
				      const panzer::PhysicsBlock& /* side_pb */,
				      const LinearObjFactory<panzer::Traits> & /* lof */,
				      const Teuchos::ParameterList& /* user_data */) const {}

    /** \brief Registers gather/orientation/DOF-projection evaluators for the side, plus a side-normal evaluator for each integration rule on the side.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] side_pb the physics block for the boundary condition's side set.
      * \param[in] lof linear object factory used to build the gather evaluators.
      * \param[in] user_data user-supplied parameters passed through to the gather/orientation evaluators.
      */
    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					           const panzer::PhysicsBlock& side_pb,
						   const LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;

      side_pb.buildAndRegisterGatherAndOrientationEvaluators(fm,lof,user_data); 
      side_pb.buildAndRegisterDOFProjectionsToIPEvaluators(fm,Teuchos::ptrFromRef(lof),user_data); 

      // add in side normals
      const std::map<int,Teuchos::RCP<panzer::IntegrationRule> > & int_rules = side_pb.getIntegrationRules();
      for(std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator itr=int_rules.begin();
          itr!=int_rules.end();++itr) {
         
        std::stringstream s;
        s << "Side Normal:" << side_pb.cellData().side();
        Teuchos::ParameterList p(s.str());
        p.set<std::string>("Name","Side Normal");
        p.set<int>("Side ID",side_pb.cellData().side());
        p.set< Teuchos::RCP<panzer::IntegrationRule> >("IR", Teuchos::rcp_const_cast<panzer::IntegrationRule>(itr->second));
        p.set<bool>("Normalize",true);

        RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));

        fm.template registerEvaluator<EvalT>(op);
      }

    }

    //@}

  private:

  };

  /**
   * \brief Builder functor used with BCStrategy_TemplateManager::buildObjects() to construct a ResponseFactory_BCStrategyAdapter for each evaluation type.
   * \sa ResponseFactory_BCStrategyAdapter
   */
  class BCStrategy_TM_ResponseAdapterBuilder {
  public:
    typedef std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > RespFact_TM_Vector;
    const BC & bc_;
    const RespFact_TM_Vector & vec_;

    /// \brief Construct from the boundary condition and the (response name, factory) pairs to adapt.
    BCStrategy_TM_ResponseAdapterBuilder(const BC & bc,const RespFact_TM_Vector & vec)
      : bc_(bc), vec_(vec) {}

    /// \brief Builds a ResponseFactory_BCStrategyAdapter<EvalT> for the given evaluation type.
    template <typename EvalT>
    Teuchos::RCP<BCStrategyBase> build() const
    { return Teuchos::rcp(new ResponseFactory_BCStrategyAdapter<EvalT>(bc_,vec_)); }
  };

  /**
   * \brief BCStrategyFactory that produces a BCStrategy_TemplateManager wrapping the responses registered for each boundary condition.
   *
   * Constructed once with the full map from boundary condition to its
   * associated response factories; buildBCStrategy() then looks up the
   * entry for a given BC and builds the adapter template manager for it.
   */
  class BCFactoryResponse : public panzer::BCStrategyFactory {
  public:
    typedef std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > RespFact_TM_Vector;
    typedef std::unordered_map<BC,Teuchos::RCP<RespFact_TM_Vector>,BC::BCHash,BC::BCEquality> BCHashMap;

    /// \brief Construct from the full map of boundary condition to its (response name, factory) pairs.
    BCFactoryResponse(const BCHashMap & hashMap)
      : hashMap_(hashMap) {}

    /// \brief Looks up \p bc in the hash map and builds a BCStrategy_TemplateManager adapting its response factories.
    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& /* global_data */) const
    {
      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcstrategy_tm
          = Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);

      BCHashMap::const_iterator itr = hashMap_.find(bc);
      TEUCHOS_ASSERT(itr!=hashMap_.end());

      BCStrategy_TM_ResponseAdapterBuilder builder(bc,*itr->second);
      bcstrategy_tm->buildObjects(builder);

      return bcstrategy_tm;
    }

  private:
    /// Map from boundary condition to the responses registered on it.
    BCHashMap hashMap_;
  };
  

} // response_bc_adapters
} // end panzer

#endif
