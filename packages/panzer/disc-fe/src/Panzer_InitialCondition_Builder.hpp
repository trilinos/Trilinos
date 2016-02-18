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

#ifndef PANZER_INITIAL_CONDITION_BUILDER_HPP
#define PANZER_INITIAL_CONDITION_BUILDER_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <vector>
#include <map>

namespace panzer {

  /** \brief Builds PHX::FieldManager objects for inital conditions and registers evaluators. 

  \param[in] wkstContainer Worksets for the corresponding evaluations.
  \param[in] physicsBlock PhysicsBlocks created by FieldManagerBuilder.
  \param[in] cm_factory Factory that provides all closure models required by the initial condition evaluation.
  \param[in] closure_models List of closure models for each element block required for initial conditions.
  \param[in] lo_factory LinearObjFactory corresponding to the problem.
  \param[in] user_data ParameterList with optional user data.
  \param[out] phx_ic_field_managers Allocated PHX::FieldManagers for each element block.

  */
  void setupInitialConditionFieldManagers(WorksetContainer & wkstContainer,
                                          const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                          const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                                          const Teuchos::ParameterList& ic_block_closure_models,
                                          const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
                                          const Teuchos::ParameterList& user_data,
                                          const bool write_graphviz_file,
                                          const std::string& graphviz_file_prefix,
                                          std::map<std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers);

  /** Structure that defines a field, including name, basis type and order.
    */
  typedef struct {
    std::string fieldName;
    std::string basisName;
    int basisOrder;
  } ICFieldDescriptor;

  /** Fill a vector with the initial conditions specified using the closure models. This is intended to work
    * with out physics blocks and so is most applicable for controls.
    *
    * \param[in] block_ids_to_cell_topo Map from element blocks to cell topologies.
    * \param[in] block_ids_to_fields Map from element blocks to required fields.
    * \param[in] wkstContainer Worksets for the corresponding evaluations.
    * \param[in] lof LinearObjFactory corresponding to the problem.
    * \param[in] cm_factory Factory that provides all closure models required by the initial condition evaluation.
    * \param[in] closure_models List of closure models for each element block required for initial conditions.
    * \param[in] user_data ParameterList with optional user data.
    * \param[in] workset_size Size of the worksets.
    * \param[in] t0 Initial time.
    * \param[out] vec Vector defined by the <code>lof</code> containing the intial conditions.
    */
  void setupControlInitialCondition(const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_ids_to_cell_topo,
                                    const std::map<std::string,std::vector<ICFieldDescriptor> > & block_ids_to_fields,
                                    panzer::WorksetContainer & wkstContainer,
                                    const panzer::LinearObjFactory<panzer::Traits> & lof,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                    const Teuchos::ParameterList & ic_closure_models,
                                    const Teuchos::ParameterList & user_data,
                                    int workset_size,
                                    double t0,
                                    const Teuchos::RCP<Thyra::VectorBase<double> > & vec);

  /** A helper function that builds the physics blocks specified by the field descriptors. This is 
    * used inside the setup function for controls to build mock physics blocks.
    *
    * \param[in] block_ids_to_cell_topo Map from element blocks to cell topologies.
    * \param[in] block_ids_to_fields Map from element blocks to required fields.
    * \param[in] workset_size Size of the worksets.
    * \param[out] physics_blocks Physics blocks used to contrust field managers.
    */
  void buildICPhysicsBlocks(const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_ids_to_cell_topo,
                            const std::map<std::string,std::vector<ICFieldDescriptor> > & block_ids_to_fields,
                            int workset_size,
                            std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physics_blocks);

  /** Execute the construction of an initial condition.
    */
  void evaluateInitialCondition(WorksetContainer & wkstContainer,
                                const std::map<std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers,
                                Teuchos::RCP<panzer::LinearObjContainer> loc,
                                const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
                                const double time_stamp);
}

#endif
