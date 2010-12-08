#ifndef PANZER_ASSEMBLY_ENGINE_TEMPLATE_MANAGER_HPP
#define PANZER_ASSEMBLY_ENGINE_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_AssemblyEngine.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

  template<typename Traits,typename LO,typename GO>
  class AssemblyEngine_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::Base,
                                panzer::AssemblyEngine<_,LO,GO> > {

  public:

    AssemblyEngine_TemplateManager() {}

    ~AssemblyEngine_TemplateManager() {}

  };

} 

#endif 
