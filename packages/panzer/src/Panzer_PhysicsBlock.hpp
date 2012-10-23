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

namespace shards {
  class CellTopology;
}

namespace panzer {
  class RegionFillData;
  class MaterialModel;
  class PureBasis;
  class EquationSetFactory;
  class GlobalData;
}

namespace panzer {

  class InputPhysicsBlock;

  class PhysicsBlock {

  public:    
    /** for testing purposes only */
    explicit PhysicsBlock() 
       : m_build_transient_support(false), m_global_data(Teuchos::null)
    { std::cout << "WARNING: Default constructor for panzer::PhysicsBlock is for testing purposes only!" << std::endl; } 

    PhysicsBlock(const panzer::InputPhysicsBlock& ipb,
                 const std::string & element_block_id,
		 const panzer::CellData & cell_data,
		 const panzer::EquationSetFactory& factory,
		 const Teuchos::RCP<panzer::GlobalData>& global_data,
		 const bool build_transient_support);

    PhysicsBlock(const panzer::PhysicsBlock & pb,
                 const panzer::CellData & cell_data,
                 const panzer::EquationSetFactory& factory);

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
    const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& getBases() const;

    const panzer::InputPhysicsBlock & getInputPhysicsBlock() const;
    
    const shards::CellTopology getBaseCellTopology() const;

    std::string physicsBlockID() const;
    std::string elementBlockID() const;

    const panzer::CellData & cellData() const;

    /** Build a copy of this physics block object but use the cell_data
      * passed in by the user. This is useful for creating physics block objects
      * that reside on the boundaries.
      */
    Teuchos::RCP<PhysicsBlock> copyWithCellData(const panzer::CellData & cell_data,
                                                const panzer::EquationSetFactory & factory) const;

    /** Copy using the "Parrot" equation sets.
      */
    Teuchos::RCP<PhysicsBlock> copyWithCellData(const panzer::CellData & cell_data) const;

    Teuchos::RCP<panzer::GlobalData> globalData() const;

    Teuchos::RCP<const FieldLibrary> getFieldLibrary() const 
    { return m_field_lib.getConst(); }

    Teuchos::RCP<const FieldLibraryBase> getFieldLibraryBase() const 
    { return m_field_lib.getConst(); }

  protected:
    void initialize(const panzer::InputPhysicsBlock & ipb,
                    const std::string & element_block_id,
   		    const panzer::CellData & cell_data,
   		    const panzer::EquationSetFactory& factory,
		    const bool build_transient_support);

    std::string m_physics_id;
    std::string m_element_block_id;
    panzer::CellData m_cell_data;
    panzer::InputPhysicsBlock m_initializer;
    const bool m_build_transient_support;
    const Teuchos::RCP<panzer::GlobalData> m_global_data;

    std::vector<std::string> m_dof_names;
    std::vector<StrPureBasisPair> m_provided_dofs;
    std::map<std::string,Teuchos::RCP<panzer::PureBasis> > m_bases;
    
    std::vector< Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > > m_equation_sets;
    Teuchos::RCP<FieldLibrary> m_field_lib;

  };
  
}

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
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    if(providedDOFs.size()==0) {
       Teuchos::RCP<IntegrationRule> intRule = eqstm.getAsObject<EvalT>()->getIntegrationRule();
       for(std::size_t i=0;i<m_provided_dofs.size();i++) {
          Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
          providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
       }
    }

    eqstm.getAsObject<EvalT>()->buildAndRegisterEquationSetEvaluators(fm, providedDOFs, user_data);

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
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    if(providedDOFs.size()==0) {
       Teuchos::RCP<IntegrationRule> intRule = eqstm.getAsObject<EvalT>()->getIntegrationRule();
       for(std::size_t i=0;i<m_provided_dofs.size();i++) {
          Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
          providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
       }
    }

    eqstm.getAsObject<EvalT>()->buildAndRegisterGatherAndOrientationEvaluators(fm, providedDOFs,lof,user_data);

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
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    if(providedDOFs.size()==0) {
       Teuchos::RCP<IntegrationRule> intRule = eqstm.getAsObject<EvalT>()->getIntegrationRule();
       for(std::size_t i=0;i<m_provided_dofs.size();i++) {
          Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
          providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
       }
    }

    eqstm.getAsObject<EvalT>()->buildAndRegisterDOFProjectionsToIPEvaluators(fm, providedDOFs,user_data);

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
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    if(providedDOFs.size()==0) {
       Teuchos::RCP<IntegrationRule> intRule = eqstm.getAsObject<EvalT>()->getIntegrationRule();
       for(std::size_t i=0;i<m_provided_dofs.size();i++) {
          Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
          providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
       }
    }

    eqstm.getAsObject<EvalT>()->buildAndRegisterScatterEvaluators(fm, providedDOFs,lof,user_data);

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
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    if(providedDOFs.size()==0) {
       Teuchos::RCP<IntegrationRule> intRule = eqstm.getAsObject<EvalT>()->getIntegrationRule();
       for(std::size_t i=0;i<m_provided_dofs.size();i++) {
          Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
          providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
       }
    }

    eqstm.getAsObject<EvalT>()->buildAndRegisterClosureModelEvaluators(fm, providedDOFs, factory, models, user_data);

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

    if(providedDOFs.size()==0) {
       Teuchos::RCP<IntegrationRule> intRule = eqstm.getAsObject<EvalT>()->getIntegrationRule();
       for(std::size_t i=0;i<m_provided_dofs.size();i++) {
          Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout(m_provided_dofs[i].second,*intRule));
          providedDOFs.push_back(std::make_pair(m_provided_dofs[i].first,basis));
       }
    }

    eqstm.getAsObject<EvalT>()->buildAndRegisterInitialConditionEvaluators(fm, providedDOFs, factory, model_name, models, lof, user_data);

  }
}

#endif
