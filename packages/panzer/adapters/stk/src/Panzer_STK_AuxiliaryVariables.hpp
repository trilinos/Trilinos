#ifndef __Panzer_STK_AuxiliaryVariables_hpp__
#define __Panzer_STK_AuxiliaryVariables_hpp__

#include "Panzer_AuxiliaryEvaluator_Factory.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer_stk {

/** This class provides a user with the opportunity
  * to inject there own (likely) mesh dependent evaluators
  * into the field manager.
  */
template <typename EvalT>
class AuxiliaryVariables : public panzer::AuxiliaryEvaluator_Factory<EvalT> {
public:
   AuxiliaryVariables(const Teuchos::RCP<const STK_Interface> & mesh,
                      const Teuchos::ParameterList & pl);
   
   virtual ~AuxiliaryVariables() {}

   virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits> & fm) const;

private:
   Teuchos::RCP<panzer::Basis> basis_;
   Teuchos::RCP<std::vector<std::string> > fieldNames_;
   Teuchos::RCP<const panzer_stk::STK_Interface> mesh_;

   AuxiliaryVariables();
   AuxiliaryVariables(const AuxiliaryVariables &);
};

}

#include "Panzer_STK_AuxiliaryVariablesT.hpp"

#endif
