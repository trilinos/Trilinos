#ifndef PANZER_ASSEMBLY_ENGINE_TEMPLATE_BUILDER_HPP
#define PANZER_ASSEMBLY_ENGINE_TEMPLATE_BUILDER_HPP

#include <string>
#include "boost/mpl/apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine.hpp"

namespace panzer {

  class AssemblyEngine_TemplateBuilder {
    
    Teuchos::RCP<panzer::FieldManagerBuilder> m_region_data;
    
  public:
    
    AssemblyEngine_TemplateBuilder(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb) :
      m_fmb(fmb) 
      {}
      
    template <typename EvalT>
      Teuchos::RCP<panzer::Base> build() const {
      return Teuchos::rcp( static_cast<panzer::Base*>
			   (new panzer::AssemblyEngine<EvalT>(m_fmb)) );
    }
    
  };
  
}

#endif 
