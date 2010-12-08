
#ifndef USER_APP_EQUATIONSET_ENERGY_H
#define USER_APP_EQUATIONSET_ENERGY_H

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace user_app {

  template <typename EvalT>
    class EquationSet_Energy : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:    

    EquationSet_Energy(const panzer::InputEquationSet& ies,
		       const panzer::CellData& cell_data);
    
    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm, const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const;

  };

}

#include "user_app_EquationSet_EnergyT.hpp"

#endif
