#ifndef __Example_ClosureModel_Factory_TemplateBuilder_hpp__
#define __Example_ClosureModel_Factory_TemplateBuilder_hpp__

#include <string>
#include "boost/mpl/apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_Base.hpp"
#include "Example_ClosureModel_Factory.hpp"

namespace Example {

class ClosureModelFactory_TemplateBuilder {
public:
    
   template <typename EvalT>
   Teuchos::RCP<panzer::Base> build() const 
   {
      return Teuchos::rcp( static_cast<panzer::Base*>(new Example::ModelFactory<EvalT>) );
   }
    
};
  
}

#endif 
