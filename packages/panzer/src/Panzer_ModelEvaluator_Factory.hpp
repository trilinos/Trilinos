#ifndef PANZER_MODEL_EVALUATOR_FACTORY_DECL_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_DECL_HPP

#include "Teuchos_RCP.hpp"

namespace Thyra {
  template<typename Scalar> class ModelEvaluator;
}

namespace panzer {
  
  template <typename ScalarT, typename LO, typename GO>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > 
  buildModelEvaluator(const RCP<panzer::FieldManagerBuilder<LO,GO>& fmb,
		      const RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof);
    
}

#include "Panzer_ModelEvaluator_Factory_impl.hpp"

#endif
