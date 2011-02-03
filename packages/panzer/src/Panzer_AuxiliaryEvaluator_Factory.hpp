#ifndef __Panzer_AuxiliaryEvaluator_Factory_hpp__
#define __Panzer_AuxiliaryEvaluator_Factory_hpp__

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Base.hpp"
#include "Panzer_Traits.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

/** Base class for auxiliary factory.
  */
class AuxiliaryEvaluator_FactoryBase {
public:
   virtual ~AuxiliaryEvaluator_FactoryBase() {}

   virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits> & fm) const = 0;
};

/** This class provides a user with the opportunity
  * to inject there own (likely) mesh dependent evaluators
  * into the field manager.
  */
template <typename EvalT>
class AuxiliaryEvaluator_Factory : public panzer::AuxiliaryEvaluator_FactoryBase {
public:
   virtual ~AuxiliaryEvaluator_Factory() {}

   virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits> & fm) const = 0;
};

}

#endif
