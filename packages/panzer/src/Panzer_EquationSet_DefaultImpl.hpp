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
    
    EquationSet_DefaultImpl(const panzer::InputEquationSet& ies, const panzer::CellData& cell_data, const bool build_transient_support);
    
    virtual ~EquationSet_DefaultImpl() {}
    
    //! Builds the integration rule, basis, DOFs, and default parameter list
    virtual void setupDOFs(int equation_dimension);

    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const = 0;

    virtual void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							 const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
                                                         const LinearObjFactory<panzer::Traits> & lof) const;
    
    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
							const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							const Teuchos::ParameterList& models) const;

    virtual void buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							    const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
							    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							    const Teuchos::ParameterList& models,
							    const LinearObjFactory<panzer::Traits> & lof,
							    const Teuchos::ParameterList& user_data) const;
    
    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const;
    
    virtual const std::vector<std::string> & getDOFNames() const;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & getProvidedDOFs() const;

  protected:
    
    const panzer::InputEquationSet m_input_eq_set;
    const panzer::CellData m_cell_data;
    const bool m_build_transient_support;
    
    std::string m_eqset_prefix;
    
    Teuchos::RCP<panzer::IntegrationRule> m_int_rule;
    Teuchos::RCP<panzer::Basis> m_basis;
    
    std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > >  m_provided_dofs;
    Teuchos::RCP< std::vector<std::string> > m_dof_names;
    Teuchos::RCP< std::vector<std::string> > m_dof_gradient_names;
    Teuchos::RCP< std::vector<std::string> > m_dof_time_derivative_names;
    Teuchos::RCP< std::vector<std::string> > m_residual_names;
    std::string m_scatter_name;
    Teuchos::RCP<Teuchos::ParameterList> m_eval_plist;

  };
  
}

#include "Panzer_EquationSet_DefaultImplT.hpp"

#endif
