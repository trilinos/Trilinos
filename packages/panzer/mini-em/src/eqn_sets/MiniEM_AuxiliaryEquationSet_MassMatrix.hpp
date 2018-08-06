#ifndef _MiniEM_AuxiliaryEquationSet_MassMatrix_hpp_
#define _MiniEM_AuxiliaryEquationSet_MassMatrix_hpp_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace mini_em {

  template <typename EvalT>
  class AuxiliaryEquationSet_MassMatrix : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:    

    AuxiliaryEquationSet_MassMatrix(const Teuchos::RCP<panzer::GlobalEvaluationDataContainer> & gedc,
                             const Teuchos::RCP<Teuchos::ParameterList>& params,
			     const int& default_integration_order,
			     const panzer::CellData& cell_data,
		             const Teuchos::RCP<panzer::GlobalData>& global_data,
		             const bool build_transient_support);
    
    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					       const panzer::FieldLibrary& /* field_library */,
                                               const Teuchos::ParameterList& /* user_data */) const;

    void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					   const panzer::FieldLibrary& field_library,
                                           const panzer::LinearObjFactory<panzer::Traits> & lof,
                                           const Teuchos::ParameterList& user_data) const;

  protected:
    std::string dof_name;
    double multiplier;
    Teuchos::RCP<const std::vector<std::string> > fieldMultipliers;
    Teuchos::RCP<panzer::GlobalEvaluationDataContainer> m_gedc;
    Teuchos::RCP<std::vector<std::string> > m_dof_names;
  };

template < >
void AuxiliaryEquationSet_MassMatrix<panzer::Traits::Jacobian>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				  const panzer::FieldLibrary& field_library,
                                  const panzer::LinearObjFactory<panzer::Traits> & lof,
                                  const Teuchos::ParameterList& user_data) const;

}

#endif
