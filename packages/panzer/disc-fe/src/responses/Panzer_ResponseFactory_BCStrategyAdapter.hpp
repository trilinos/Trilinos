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

// This file is used to make a set of ResponseEvaluatorFactory objects look 
// like BCStrategy objects. It is only being used by the ResponseLibrary and
// is primarily there to unify the interface between volumetric response objects
// and surface response objects. It is a little bit painful!

namespace panzer {
namespace response_bc_adapters {

  // An adapter that turns a ResponseEvaluatorFactory object into a
  // BCStrategy object.
  template <typename EvalT>
  class ResponseFactory_BCStrategyAdapter : public panzer::BCStrategy<EvalT>
  {
    std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > refVec_;

  public:    
    
    ResponseFactory_BCStrategyAdapter(const panzer::BC & bc,const std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > & refVec)
      : panzer::BCStrategy<EvalT>(bc), refVec_(refVec) {}
    
    virtual ~ResponseFactory_BCStrategyAdapter() {}
    
    //! \name Derived from BCStrategy
    //@{ 

    virtual void setup(const panzer::PhysicsBlock& /* side_pb */, const Teuchos::ParameterList& /* user_data */) {}
      
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

    virtual void 
    buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
				      const panzer::PhysicsBlock& /* side_pb */,
				      const LinearObjFactory<panzer::Traits> & /* lof */,
				      const Teuchos::ParameterList& /* user_data */) const {}

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

  class BCStrategy_TM_ResponseAdapterBuilder {
  public:
    typedef std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > RespFact_TM_Vector;
    const BC & bc_;
    const RespFact_TM_Vector & vec_;
 
    BCStrategy_TM_ResponseAdapterBuilder(const BC & bc,const RespFact_TM_Vector & vec)
      : bc_(bc), vec_(vec) {}

    template <typename EvalT>
    Teuchos::RCP<BCStrategyBase> build() const
    { return Teuchos::rcp(new ResponseFactory_BCStrategyAdapter<EvalT>(bc_,vec_)); }
  };

  class BCFactoryResponse : public panzer::BCStrategyFactory {
  public:
    typedef std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<panzer::Traits> > > > RespFact_TM_Vector;
    typedef std::unordered_map<BC,Teuchos::RCP<RespFact_TM_Vector>,BC::BCHash,BC::BCEquality> BCHashMap;

    BCFactoryResponse(const BCHashMap & hashMap)
      : hashMap_(hashMap) {}

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
    BCHashMap hashMap_;
  };
  

} // response_bc_adapters
} // end panzer

#endif
