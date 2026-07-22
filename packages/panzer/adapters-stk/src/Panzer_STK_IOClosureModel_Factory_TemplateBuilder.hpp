// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_STK_IOClosureModel_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"

namespace panzer_stk {

  template <typename TraitsT>
  class IOClosureModelFactory_TemplateBuilder {

  public:
    IOClosureModelFactory_TemplateBuilder(const panzer::ClosureModelFactory_TemplateManager<TraitsT> & cmf_tm,
                                          const Teuchos::RCP<STK_Interface> & mesh,
                                          const Teuchos::ParameterList & outputList)
       : cmf_tm_(cmf_tm), mesh_(mesh), outputList_(outputList), plConstr_(true) {}

    IOClosureModelFactory_TemplateBuilder(const panzer::ClosureModelFactory_TemplateManager<TraitsT> & cmf_tm,
                                          const Teuchos::RCP<STK_Interface> & mesh,
                                          const std::map<std::string,std::vector<std::string> > & nodalFields,
                                          const std::map<std::string,std::vector<std::string> > & cellFields)

       : cmf_tm_(cmf_tm), mesh_(mesh), nodalFields_(nodalFields), cellFields_(cellFields), plConstr_(false) {}
    
    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      if(plConstr_)
        return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
                            (new panzer_stk::IOClosureModelFactory<EvalT>(cmf_tm_.template getAsObject<EvalT>(),mesh_,outputList_)) );
      else
        return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
                            (new panzer_stk::IOClosureModelFactory<EvalT>(cmf_tm_.template getAsObject<EvalT>(),mesh_,nodalFields_,cellFields_)) );
    }
    
  private:
     const panzer::ClosureModelFactory_TemplateManager<TraitsT> & cmf_tm_;
     Teuchos::RCP<STK_Interface> mesh_;
     Teuchos::ParameterList outputList_;
     std::map<std::string,std::vector<std::string> > nodalFields_;
     std::map<std::string,std::vector<std::string> > cellFields_;
     bool plConstr_;
  };
  
}

#endif 
