#ifndef __Panzer_AuxiliaryEvaluator_Factory_TemplateManager_hpp__
#define __Panzer_AuxiliaryEvaluator_Factory_TemplateManager_hpp__

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Panzer_VectorTemplateManager.hpp"
#include "Panzer_AuxiliaryEvaluator_Factory.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

template <typename Traits>
class AuxiliaryEvaluator_TemplateManager : 
   public VectorTemplateManager<typename Traits::EvalTypes,
                                panzer::Base,
                                panzer::AuxiliaryEvaluator_Factory<_> > {
public:
   AuxiliaryEvaluator_TemplateManager() {}

   ~AuxiliaryEvaluator_TemplateManager() {}
};

}

#endif
