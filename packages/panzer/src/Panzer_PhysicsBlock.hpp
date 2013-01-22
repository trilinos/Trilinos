// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

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
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_FieldLibrary.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace shards {
  class CellTopology;
}

namespace panzer {
  class RegionFillData;
  class MaterialModel;
  class PureBasis;
  class IntegrationRule;
  class EquationSetFactory;
  class GlobalData;
  class PhysicsBlock;
}

namespace panzer {


  /** Non-member function for building the physics blocks
      \relates panzer::PhysicsBlock 
  */
  void buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                          const std::map<std::string,Teuchos::RCP<const shards::CellTopology> >& block_ids_to_cell_topo,
                          const Teuchos::RCP<Teuchos::ParameterList>& physics_blocks_plist,
			  const int default_integration_order,
                          const std::size_t workset_size,
                          const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
			  const Teuchos::RCP<panzer::GlobalData>& global_data,
                          const bool build_transient_support,
                          std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks);
  
  //! Object that contains information on the physics and discretization of a block of elements with the SAME topology.
  class PhysicsBlock {

  public:    
    /** for testing purposes only */
    explicit PhysicsBlock() 
       : m_build_transient_support(false), m_global_data(Teuchos::null)
    { std::cout << "WARNING: Default constructor for panzer::PhysicsBlock is for testing purposes only!" << std::endl; } 

    PhysicsBlock(const Teuchos::RCP<Teuchos::ParameterList>& physics_block_plist,
                 const std::string & element_block_id,
		 const int default_integration_order,
		 const panzer::CellData & cell_data,
		 const Teuchos::RCP<const panzer::EquationSetFactory>& factory,
		 const Teuchos::RCP<panzer::GlobalData>& global_data,
		 const bool build_transient_support);

    PhysicsBlock(const panzer::PhysicsBlock & pb,
                 const panzer::CellData & cell_data);

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					       const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const panzer::LinearObjFactory<panzer::Traits> & lof,
							const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						      const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					   const panzer::LinearObjFactory<panzer::Traits> & lof,
					   const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
						const Teuchos::ParameterList& models,
						const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
						    const std::string& model_name,
						    const Teuchos::ParameterList& models,
						    const panzer::LinearObjFactory<panzer::Traits> & lof,
						    const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				            const std::string& model_name,
					    const Teuchos::ParameterList& models,
					    const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterEquationSetEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
						      const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterGatherAndOrientationEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
							       const LinearObjFactory<panzer::Traits> & lof,
							       const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterDOFProjectionsToIPEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
							     const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterScatterEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
						  const LinearObjFactory<panzer::Traits> & lof,
						  const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterClosureModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
						       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
						       const Teuchos::ParameterList& models,
						       const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterInitialConditionEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
							   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							   const std::string& model_name,
							   const Teuchos::ParameterList& models,
							   const panzer::LinearObjFactory<panzer::Traits> & lof,
							   const Teuchos::ParameterList& user_data) const;

    const std::vector<std::string>& getDOFNames() const;
    const std::vector<StrPureBasisPair>& getProvidedDOFs() const;

    //! Returns the unique set of bases, key is the unique panzer::PureBasis::name() of the basis
    const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& getBases() const;

    //! Returns the unique set of point rules, key is the unique panzer::PointRule::name()
    const std::map<int,Teuchos::RCP<panzer::IntegrationRule> >& getIntegrationRules() const;

    const shards::CellTopology getBaseCellTopology() const;

    std::string physicsBlockID() const;
    std::string elementBlockID() const;

    const panzer::CellData & cellData() const;

    /** Build a copy of this physics block object but use the cell_data
      * passed in by the user. This is useful for creating physics block objects
      * that reside on the boundaries from physics blocks created for the volume.
      */
    Teuchos::RCP<PhysicsBlock> copyWithCellData(const panzer::CellData & cell_data) const;

    Teuchos::RCP<panzer::GlobalData> globalData() const;

    Teuchos::RCP<const FieldLibrary> getFieldLibrary() const 
    { return m_field_lib.getConst(); }
    
    Teuchos::RCP<const FieldLibraryBase> getFieldLibraryBase() const 
    { return m_field_lib.getConst(); }

  protected:
    void initialize(const Teuchos::RCP<Teuchos::ParameterList>& input_parameters,
		    const int& default_integration_order,
                    const std::string & element_block_id,
   		    const panzer::CellData & cell_data,
		    const bool build_transient_support);

    std::string m_physics_id;
    std::string m_element_block_id;
    int m_default_integration_order;
    panzer::CellData m_cell_data;
    //! store the input parameter list for copy ctors
    Teuchos::RCP<Teuchos::ParameterList> m_input_parameters;
    const bool m_build_transient_support;
    const Teuchos::RCP<panzer::GlobalData> m_global_data;

    std::vector<std::string> m_dof_names;
    std::vector<StrPureBasisPair> m_provided_dofs;
    //! map of unique bases, key is the panzer::PureBasis::name() corresponding to its value
    std::map<std::string,Teuchos::RCP<panzer::PureBasis> > m_bases;
    //! map of unique integration rules, key is panzer::IntegrationRule::order() corresponding to its value
    std::map<int,Teuchos::RCP<panzer::IntegrationRule> > m_integration_rules;
    
    std::vector< Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > > m_equation_sets;
    Teuchos::RCP<FieldLibrary> m_field_lib;
    Teuchos::RCP<const panzer::EquationSetFactory> m_eqset_factory;
  };
  
}

// ************************************************************
// template implementations
// ************************************************************

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterEquationSetEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
									const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    // Loop over integration rules
    for (std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator ir_iter = m_integration_rules.begin();
	 ir_iter != m_integration_rules.end(); ++ ir_iter) {
      
      Teuchos::RCP<panzer::IntegrationRule> ir = ir_iter->second;
      
      eqstm.getAsObject<EvalT>()->buildAndRegisterEquationSetEvaluators(fm, *m_field_lib->buildFieldLayoutLibrary(*ir), ir, user_data);
    }

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterGatherAndOrientationEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
										 const LinearObjFactory<panzer::Traits> & lof,
										 const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    eqstm.getAsObject<EvalT>()->buildAndRegisterGatherAndOrientationEvaluators(fm,*m_field_lib,lof,user_data);

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterDOFProjectionsToIPEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
									       const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    // Loop over integration rules
    for (std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator ir_iter = m_integration_rules.begin();
	 ir_iter != m_integration_rules.end(); ++ ir_iter) {
      
      Teuchos::RCP<panzer::IntegrationRule> ir = ir_iter->second;
      
      eqstm.getAsObject<EvalT>()->buildAndRegisterDOFProjectionsToIPEvaluators(fm,*m_field_lib->buildFieldLayoutLibrary(*ir),ir,user_data);

    }

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterScatterEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
								    const LinearObjFactory<panzer::Traits> & lof,
								    const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    eqstm.getAsObject<EvalT>()->buildAndRegisterScatterEvaluators(fm,*m_field_lib,lof,user_data);

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterClosureModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
									 const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
									 const Teuchos::ParameterList& models,
									 const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    // Loop over integration rules
    for (std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator ir_iter = m_integration_rules.begin();
	 ir_iter != m_integration_rules.end(); ++ ir_iter) {
      
      Teuchos::RCP<panzer::IntegrationRule> ir = ir_iter->second;
      
      eqstm.getAsObject<EvalT>()->buildAndRegisterClosureModelEvaluators(fm,*m_field_lib->buildFieldLayoutLibrary(*ir),ir,factory,models,user_data);
    }

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterInitialConditionEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
									     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
									     const std::string& model_name,
									     const Teuchos::ParameterList& models,
									     const panzer::LinearObjFactory<panzer::Traits> & lof,
									     const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    eqstm.getAsObject<EvalT>()->buildAndRegisterInitialConditionEvaluators(fm, *m_field_lib, factory, model_name, models, lof, user_data);

  }
}

#endif
