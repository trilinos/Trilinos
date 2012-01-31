#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_DECL_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_DECL_HPP

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
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  private:
    void parseOutputList(const Teuchos::ParameterList & pl);

    //! Mesh pointer, will be passed around
    Teuchos::RCP<STK_Interface> mesh_;
 
    //! Map showing which fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToFields_;

    /** Map stating if an evaluator for a particular block ID has been included.
      *
      * This is a bit of hack that is done to gurantee only one evaluator is
      * added to each field manager. However, if an instantiation of this closure model factory is
      * used in multiple places then the appropriate evaluator will be added and
      * required only once. So its likely that not every field manager will have (and require)
      * the scatter cell evaluators.
      */    
    mutable std::map<std::string,bool> blockIdEvaluated_;

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
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

}

#endif
