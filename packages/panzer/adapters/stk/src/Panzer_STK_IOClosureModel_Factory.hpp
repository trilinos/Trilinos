#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_HPP

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory.hpp"

#include "Panzer_STK_Interface.hpp"

#include <vector>
#include <string>

namespace panzer {
  class InputEquationSet;
}

namespace panzer_stk {

  template<typename EvalT>
  class IOClosureModelFactory : public panzer::ClosureModelFactory<EvalT> {
  public:

    IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_,
                          const Teuchos::RCP<STK_Interface> & mesh,
                          const Teuchos::ParameterList & outputList);

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const panzer::InputEquationSet& set,
		       const Teuchos::ParameterList& models,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  private:
    void parseOutputList(const Teuchos::ParameterList & pl);

    //! Mesh pointer, will be passed around
    Teuchos::RCP<STK_Interface> mesh_;
 
    //! Map showing which fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToFields_;

    //! we will reuse the drekar closure model factory
    Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_;
  };

  template < >
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
  panzer_stk::IOClosureModelFactory<panzer::Traits::Residual>::buildClosureModels(const std::string& model_id,
		       const panzer::InputEquationSet& set,
		       const Teuchos::ParameterList& models,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

}

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION
#include "Panzer_STK_IOClosureModel_FactoryT.hpp"
#endif

#endif
