#ifndef PANZER_EQUATION_SET_DEFAULTIMPL_H
#define PANZER_EQUATION_SET_DEFAULTIMPL_H

#include "Panzer_EquationSet.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"

namespace PHX {
  template<typename Traits> class FieldManager;
}

namespace panzer {

  template <typename EvalT>
  class EquationSet_DefaultImpl : public panzer::EquationSet<EvalT> {
    
  public:    
    
    EquationSet_DefaultImpl(const panzer::InputEquationSet& ies, const panzer::CellData& cell_data);
    
    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const = 0;

    virtual void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							 const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const = 0;
    
    // virtual void buildAndRegisterMaterialModelEvaluators(int physics_id, PHX::FieldManager<panzer::Traits>& fm, const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const;
    
    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const = 0;
    
    virtual const std::vector<std::string> & getDOFNames() const = 0;
    
  protected:
    
    const panzer::InputEquationSet m_input_eq_set;
    const panzer::CellData m_cell_data;
  };
  
}

#include "Panzer_EquationSet_DefaultImplT.hpp"

#endif
