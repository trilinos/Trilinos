#ifndef PANZER_EQUATION_SET_HPP
#define PANZER_EQUATION_SET_HPP

#include "Panzer_EquationSet_Base.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_CellData.hpp"

namespace PHX {
  template<typename Traits> class FieldManager;
}

namespace panzer {

  template <typename EvalT>
  class EquationSet : public panzer::EquationSetBase {
    
  public:    
    
    EquationSet() {}

    virtual ~EquationSet() {}

    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const = 0;

    virtual void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							 const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const = 0;
    
    virtual void buildAndRegisterModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						 const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
						 const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
						 const std::vector<Teuchos::ParameterList>& models) const = 0;
    
    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const = 0;
    
    virtual const std::vector<std::string> & getDOFNames() const = 0;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & getProvidedDOFs() const = 0;

  };
  
}

#endif
