#ifndef PANZER_PHYSICS_BLOCK_HPP
#define PANZER_PHYSICS_BLOCK_HPP

#include <string>
#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_EquationSet.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"

namespace shards {
  class CellTopology;
}

namespace panzer {
  class RegionFillData;
  class MaterialModel;
  class Basis;
  class EquationSetFactory;
}

namespace panzer {

  class InputPhysicsBlock;

  class PhysicsBlock {

  public:    
    typedef std::pair<std::string,Teuchos::RCP<panzer::Basis> > StrBasisPair;
    
    PhysicsBlock(const panzer::InputPhysicsBlock& ipb,
                 const std::string & element_block_id,
		 const panzer::CellData cell_data,
		 const panzer::EquationSetFactory& factory);

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm) const;

    void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm) const;

    void buildAndRegisterModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					 const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
					 const std::vector<Teuchos::ParameterList>& models) const;

    template<typename EvalT>
    void buildAndRegisterEquationSetEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm) const;

    template<typename EvalT>
    void buildAndRegisterGatherScatterEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm) const;

    template<typename EvalT>
    void buildAndRegisterModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
						const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
						const std::vector<Teuchos::ParameterList>& models) const;

    const std::vector<std::string>& getDOFNames() const;
    const std::vector<StrBasisPair>& getProvidedDOFs() const;
    const std::map<std::string,Teuchos::RCP<panzer::Basis> >& getBases() const;
    
    const shards::CellTopology getBaseCellTopology() const;

    std::string physicsBlockID() const;
    std::string elementBlockID() const;

    const panzer::CellData cellData() const;

  protected:

    std::string m_physics_id;
    std::string m_element_block_id;
    panzer::CellData m_cell_data;

    std::vector<std::string> m_dof_names;
    std::vector<StrBasisPair> m_provided_dofs;
    std::map<std::string,Teuchos::RCP<panzer::Basis> > m_bases;
    
    std::vector< Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > > m_equation_sets;

  };
  
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterEquationSetEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    eqstm.getAsObject<EvalT>()->buildAndRegisterEquationSetEvaluators(fm, m_provided_dofs);

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterGatherScatterEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    eqstm.getAsObject<EvalT>()->buildAndRegisterGatherScatterEvaluators(fm, m_provided_dofs);

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
								  const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
								  const std::vector<Teuchos::ParameterList>& models) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    eqstm.getAsObject<EvalT>()->buildAndRegisterModelEvaluators(fm, m_provided_dofs, factories, models);

  }
}

#endif
