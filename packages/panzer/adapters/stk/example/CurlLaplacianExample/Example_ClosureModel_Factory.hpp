#ifndef __Example_ClosureModelFactory_hpp__
#define __Example_ClosureModelFactory_hpp__

#include "Panzer_ClosureModel_Factory.hpp"

namespace panzer {
  class InputEquationSet;
}

namespace Example {

template<typename EvalT>
class ModelFactory : public panzer::ClosureModelFactory<EvalT> {

public:

   Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
   buildClosureModels(const std::string& model_id,
                      const panzer::InputEquationSet& set,
                      const Teuchos::ParameterList& models,
                      const Teuchos::ParameterList& default_params,
                      const Teuchos::ParameterList& user_data,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      PHX::FieldManager<panzer::Traits>& fm) const;

};

}

#include "Example_ClosureModel_Factory_impl.hpp"

#endif
