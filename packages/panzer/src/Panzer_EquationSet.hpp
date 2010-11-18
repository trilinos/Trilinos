#ifndef PANZER_EQUATION_SET_H
#define PANZER_EQUATION_SET_H

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
    
   //  virtual void buildAndRegisterMaterialModelEvaluators(int physics_id, PHX::FieldManager<panzer::Traits>& fm, const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const;
    
    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const = 0;
    
    virtual const std::vector<std::string> & getDOFNames() const = 0;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & getProvidedDOFs() const = 0;

  };
  
}

#endif
