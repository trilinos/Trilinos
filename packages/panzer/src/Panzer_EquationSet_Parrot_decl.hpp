#ifndef __Panzer_EquationSet_Parrot_decl_hpp__
#define __Panzer_EquationSet_Parrot_decl_hpp__

#include <vector>
#include <string>

#include "Panzer_config.hpp"

#include "Panzer_EquationSet_Base.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"

namespace panzer {

  /** This class allows the creation of an equation set that contains
    * no extra evalutors nor supports operations but yet provides access
    * to all the required gathered fields of a previously constructed
    * equation set. Thus it is a "Parrot" equation set. This class implements
    * the default equation set and provides gather operations, but ignores
    * any scatter operations. It is useful for unifying surface responses
    * and volume response interfaces.
    */
  template <typename EvalT>
  class EquationSet_Parrot : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:    

    EquationSet_Parrot(const EquationSetBase & eqSet,
                       const panzer::InputEquationSet& ies,
                       const panzer::CellData& cell_data,
                       const Teuchos::RCP<panzer::GlobalData>& global_data,
                       const bool build_transient_support);
    
    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm, 
                                               const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs, 
                                               const Teuchos::ParameterList& user_data) const {}

  };

}

#endif
