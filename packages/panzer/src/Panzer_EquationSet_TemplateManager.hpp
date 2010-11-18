#ifndef PANZER_EQUATION_SET_TEMPLATE_MANAGER_H
#define PANZER_EQUATION_SET_TEMPLATE_MANAGER_H

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_EquationSet.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

  template<typename Traits>
  class EquationSet_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::EquationSetBase,
                                panzer::EquationSet<_> > {

  public:

    EquationSet_TemplateManager() {}

    ~EquationSet_TemplateManager() {}

  };

} 

#endif 
