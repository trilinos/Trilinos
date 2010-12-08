#ifndef PANZER_EQUATION_SET_DEFAULTIMPL_HPP
#define PANZER_EQUATION_SET_DEFAULTIMPL_HPP

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
    
    virtual ~EquationSet_DefaultImpl() {}
    
    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const = 0;

    virtual void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							 const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const;
    
    virtual void buildAndRegisterModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						 const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
						 const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
						 const std::vector<Teuchos::ParameterList>& models) const;
    
    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const;
    
    virtual const std::vector<std::string> & getDOFNames() const;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & getProvidedDOFs() const;

  protected:
    
    const panzer::InputEquationSet m_input_eq_set;
    const panzer::CellData m_cell_data;
    
    std::string m_eqset_prefix;
    
    Teuchos::RCP<panzer::IntegrationRule> m_int_rule;
    Teuchos::RCP<panzer::Basis> m_basis;
    
    std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > >  m_provided_dofs;
    Teuchos::RCP< std::vector<std::string> > m_dof_names;
    Teuchos::RCP< std::vector<std::string> > m_dof_gradient_names;
    Teuchos::RCP< std::vector<std::string> > m_residual_names;
    std::string m_scatter_name;
    Teuchos::RCP<Teuchos::ParameterList> m_eval_plist;

  };
  
}

#include "Panzer_EquationSet_DefaultImplT.hpp"

#endif
