#ifndef __Panzer_AuxiliaryEvaluator_Factory_hpp__
#define __Panzer_AuxiliaryEvaluator_Factory_hpp__

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Base.hpp"
#include "Panzer_Traits.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

/** This class provides a user with the opportunity
  * to inject there own (likely) mesh dependent evaluators
  * into the field manager.
  */
template <typename EvalT>
class AuxiliaryEvaluator_Factory : public panzer::Base {
public:
   virtual ~AuxiliaryEvaluator_Factory() {}

   virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits> & fm) = 0;
};

}

#endif
