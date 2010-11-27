#ifndef PANZER_EQUATION_SET_BASE_HPP
#define PANZER_EQUATION_SET_BASE_HPP

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ModelFactory_TemplateManager.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  class Basis; // forward declaration for getProvidedDOFs();

  //! Non-templated empty base class for EquationSet objects
  class EquationSetBase {
    
  public:
    
    EquationSetBase() {}
    
    virtual ~EquationSetBase() {}
    
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
