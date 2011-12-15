#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP

#include <string>
#include "boost/mpl/apply.hpp"
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
       : cmf_tm_(cmf_tm), mesh_(mesh), outputList_(outputList) {}
    
    template <typename EvalT>
    Teuchos::RCP<panzer::Base> build() const {
      return Teuchos::rcp( static_cast<panzer::Base*>
			   (new panzer_stk::IOClosureModelFactory<EvalT>(cmf_tm_.template getAsObject<EvalT>(),mesh_,outputList_)) );
    }
    
  private:
     const panzer::ClosureModelFactory_TemplateManager<TraitsT> & cmf_tm_;
     Teuchos::RCP<STK_Interface> mesh_;
     Teuchos::ParameterList outputList_;
   
  };
  
}

#endif 
