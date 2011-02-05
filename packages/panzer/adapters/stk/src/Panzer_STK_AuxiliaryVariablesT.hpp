#include "Panzer_STK_GatherFields.hpp"
#include "Panzer_Basis.hpp"

namespace panzer_stk {

/** This class provides a user with the opportunity
  * to inject there own (likely) mesh dependent evaluators
  * into the field manager.
  */
template <typename EvalT>
AuxiliaryVariables<EvalT>::AuxiliaryVariables(const Teuchos::RCP<const STK_Interface> & mesh,
                                       const Teuchos::ParameterList & pl)
{
   // read in parameter lists
   // determine values and parameters to read

   basis_ = pl.get<Teuchos::RCP<panzer::Basis> > ("Basis");
   fieldNames_ = pl.get<Teuchos::RCP<std::vector<std::string> > >("Field Names");
   mesh_ = mesh;
}

template <typename EvalT>
void AuxiliaryVariables<EvalT>::buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits> & fm) const
{
   Teuchos::ParameterList pl;
   pl.set("Basis",basis_);
   pl.set("Field Names",fieldNames_);
  
   Teuchos::RCP< PHX::Evaluator<panzer::Traits> > op =
      Teuchos::rcp(new panzer_stk::GatherFields<EvalT,panzer::Traits>(mesh_,pl));

   fm.template registerEvaluator<EvalT>(op);

//   for(std::vector<Teuchos::RCP<PHX::FieldTag> >::const_iterator itr=op->evaluatedFields().begin();
//       itr!=op->evaluatedFields().end();++itr) {
//      fm.requireField<EvalT>(**itr);
//   }
}

}
