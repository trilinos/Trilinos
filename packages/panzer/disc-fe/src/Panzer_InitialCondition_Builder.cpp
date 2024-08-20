// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_InitialCondition_Builder.hpp"

#include "Teuchos_Assert.hpp"

#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

namespace panzer {

void 
setupInitialConditionFieldManagers(WorksetContainer & wkstContainer,
                                   const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                                   const Teuchos::ParameterList& ic_block_closure_models,
                                   const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
                                   const Teuchos::ParameterList& user_data,
                                   const bool write_graphviz_file,
                                   const std::string& graphviz_file_prefix,
                                   std::map< std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers)
{
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
    Teuchos::RCP<panzer::PhysicsBlock> pb = *blkItr;
    std::string blockId = pb->elementBlockID();

    // build a field manager object
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
          = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    
    // Choose model sublist for this element block
    std::string closure_model_name = "";
    if (ic_block_closure_models.isSublist(blockId))
      closure_model_name = blockId;
    else if (ic_block_closure_models.isSublist("Default"))
      closure_model_name = "Default";
    else 
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to find initial condition for element block \"" << blockId 
                                                      << "\".  You must provide an initial condition for each element block or set a default!" 
                                                      << ic_block_closure_models);

    // Only the residual is used by initial conditions
    std::vector<bool> active_evaluation_types(Sacado::mpl::size<panzer::Traits::EvalTypes>::value,false);
    int residual_index = Sacado::mpl::find<panzer::Traits::EvalTypes,panzer::Traits::Residual>::value;
    active_evaluation_types[residual_index] = true;
    pb->setActiveEvaluationTypes(active_evaluation_types);

    // use the physics block to register evaluators
    pb->buildAndRegisterInitialConditionEvaluators(*fm, cm_factory, closure_model_name, ic_block_closure_models, lo_factory, user_data);

    pb->activateAllEvaluationTypes();

    // build the setup data using passed in information
    Traits::SD setupData;
    const WorksetDescriptor wd = blockDescriptor(blockId);
    setupData.worksets_ = wkstContainer.getWorksets(wd);
    setupData.orientations_ = wkstContainer.getOrientations();

    int i=0;
    for (auto eval_type=fm->begin(); eval_type != fm->end(); ++eval_type,++i) {
      if (active_evaluation_types[i])
        eval_type->postRegistrationSetup(setupData,*fm,false,false,nullptr);
    }

    phx_ic_field_managers[blockId] = fm;

    if (write_graphviz_file)
      fm->writeGraphvizFile(graphviz_file_prefix+"_IC_"+blockId);
  }
}

void 
setupInitialConditionFieldManagers(WorksetContainer & wkstContainer,
                                   const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                                   const Teuchos::ParameterList& closure_models,
                                   const Teuchos::ParameterList& ic_block_closure_models,
                                   const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
                                   const Teuchos::ParameterList& user_data,
                                   const bool write_graphviz_file,
                                   const std::string& graphviz_file_prefix,
                                   std::map< std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers)
{
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
    Teuchos::RCP<panzer::PhysicsBlock> pb = *blkItr;
    std::string blockId = pb->elementBlockID();

    // build a field manager object
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
          = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    
    // Choose model sublist for this element block
    std::string closure_model_name = "";
    if (ic_block_closure_models.isSublist(blockId))
      closure_model_name = blockId;
    else if (ic_block_closure_models.isSublist("Default"))
      closure_model_name = "Default";
    else 
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to find initial condition for element block \"" << blockId 
                                                      << "\".  You must provide an initial condition for each element block or set a default!" 
                                                      << ic_block_closure_models);

    // Only the residual is used by initial conditions
    std::vector<bool> active_evaluation_types(Sacado::mpl::size<panzer::Traits::EvalTypes>::value,false);
    int residual_index = Sacado::mpl::find<panzer::Traits::EvalTypes,panzer::Traits::Residual>::value;
    active_evaluation_types[residual_index] = true;
    pb->setActiveEvaluationTypes(active_evaluation_types);

    // build and register all closure models
    pb->buildAndRegisterClosureModelEvaluators(*fm,cm_factory,closure_models,user_data);
     
    // use the physics block to register evaluators
    pb->buildAndRegisterInitialConditionEvaluators(*fm, cm_factory, closure_model_name, ic_block_closure_models, lo_factory, user_data);

    pb->activateAllEvaluationTypes();

    // build the setup data using passed in information
    Traits::SD setupData;
    const WorksetDescriptor wd = blockDescriptor(blockId);
    setupData.worksets_ = wkstContainer.getWorksets(wd);
    setupData.orientations_ = wkstContainer.getOrientations();

    int i=0;
    for (auto eval_type=fm->begin(); eval_type != fm->end(); ++eval_type,++i) {
      if (active_evaluation_types[i])
        eval_type->postRegistrationSetup(setupData,*fm,false,false,nullptr);
    }

    phx_ic_field_managers[blockId] = fm;

    if (write_graphviz_file)
      fm->writeGraphvizFile(graphviz_file_prefix+"_IC_"+blockId);
  }
}

void 
evaluateInitialCondition(WorksetContainer & wkstContainer,
                         const std::map< std::string,Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers,
                         Teuchos::RCP<panzer::LinearObjContainer> loc,
                         const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
                         const double time_stamp)
{
  typedef LinearObjContainer LOC;
  panzer::Traits::PED ped;

  // allocate a ghosted container for the initial condition
  Teuchos::RCP<LOC> ghostedloc = lo_factory.buildGhostedLinearObjContainer();
  lo_factory.initializeGhostedContainer(LOC::X,*ghostedloc);

  // allocate a counter to keep track of where this processor set initial conditions
  Teuchos::RCP<LOC> localCounter = lo_factory.buildPrimitiveGhostedLinearObjContainer();
  Teuchos::RCP<LOC> globalCounter = lo_factory.buildPrimitiveLinearObjContainer();
  Teuchos::RCP<LOC> summedGhostedCounter = lo_factory.buildPrimitiveGhostedLinearObjContainer();

  lo_factory.initializeGhostedContainer(LOC::F,*localCounter); // store counter in F
  localCounter->initialize();

  ped.gedc->addDataObject("Residual Scatter Container",ghostedloc);
  ped.gedc->addDataObject("Dirichlet Counter",localCounter);
  ped.first_sensitivities_name = "";

  for(std::map< std::string,Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >::const_iterator itr=phx_ic_field_managers.begin();
      itr!=phx_ic_field_managers.end();++itr) {
    std::string blockId = itr->first;
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = itr->second;

    fm->preEvaluate<panzer::Traits::Residual>(ped);

    // Loop over worksets in this element block
    const WorksetDescriptor wd = blockDescriptor(blockId);
    std::vector<panzer::Workset>& w = *wkstContainer.getWorksets(wd);
    for (std::size_t i = 0; i < w.size(); ++i) {
      panzer::Workset& workset = w[i];
      workset.time = time_stamp;
      
      fm->evaluateFields<panzer::Traits::Residual>(workset);
    }
  }

  lo_factory.initializeGhostedContainer(LOC::F,*summedGhostedCounter); // store counter in F
  summedGhostedCounter->initialize();

  // do communication to build summed ghosted counter for dirichlet conditions
  {
    lo_factory.initializeContainer(LOC::F,*globalCounter); // store counter in F
    globalCounter->initialize();
    lo_factory.ghostToGlobalContainer(*localCounter,*globalCounter,LOC::F);
        // Here we do the reduction across all processors so that the number of times
        // a dirichlet condition is applied is summed into the global counter

   lo_factory.globalToGhostContainer(*globalCounter,*summedGhostedCounter,LOC::F);
        // finally we move the summed global vector into a local ghosted vector
        // so that the dirichlet conditions can be applied to both the ghosted
        // right hand side and the ghosted matrix
  }

  panzer::GlobalEvaluationDataContainer gedc;
  gedc.addDataObject("Residual Scatter Container",ghostedloc);

  // adjust ghosted system for boundary conditions
  for(GlobalEvaluationDataContainer::iterator itr=gedc.begin();itr!=gedc.end();itr++) {
    Teuchos::RCP<LOC> loc2 = Teuchos::rcp_dynamic_cast<LOC>(itr->second);
    if(loc2!=Teuchos::null) {
      bool zeroVectorRows = false;
      bool adjustX = true;
      lo_factory.adjustForDirichletConditions(*localCounter,*summedGhostedCounter,*loc2, zeroVectorRows, adjustX);
    }
  }

  // gather the ghosted solution back into the global container, which performs the sum
  lo_factory.ghostToGlobalContainer(*ghostedloc,*loc,LOC::X);
}

// This is an anonymous namespace that hides a few helper classes for the control
// initial condition construction. In particual an IC Equation set is defined that
// is useful for building initial condition vectors.
namespace {

template <typename EvalT>
class EquationSet_IC : public EquationSet_DefaultImpl<EvalT> {
public:    

   /** In the constructor you set all the fields provided by this
     * equation set. 
     */
   EquationSet_IC(const Teuchos::RCP<Teuchos::ParameterList>& params,
                  const int& default_integration_order,
                  const CellData& cell_data,
                  const Teuchos::RCP<GlobalData>& global_data,
                  const bool build_transient_support);
    
   /** The specific evaluators are registered with the field manager argument.
     */
   void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<Traits>& /* fm */,
                                              const FieldLibrary& /* field_library */,
                                              const Teuchos::ParameterList& /* user_data */) const {}

};

// ***********************************************************************
template <typename EvalT>
EquationSet_IC<EvalT>::
EquationSet_IC(const Teuchos::RCP<Teuchos::ParameterList>& params,
               const int& default_integration_order,
               const CellData& cell_data,
               const Teuchos::RCP<GlobalData>& global_data,
               const bool build_transient_support) :
  EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support )
{
  // ********************
  // Validate and parse parameter list
  // ********************
  Teuchos::ParameterList valid_parameters_sublist;
  valid_parameters_sublist.set("Basis Type","HGrad","Type of Basis to use");
  valid_parameters_sublist.set("Basis Order",1,"Order of the basis");

  for(auto itr=params->begin();itr!=params->end();++itr) {
     
    const std::string field = params->name(itr);
    const Teuchos::ParameterEntry & entry = params->entry(itr);

    // only process lists
    if(!entry.isList()) continue;

    Teuchos::ParameterList & basisPL = entry.getValue((Teuchos::ParameterList *) 0);
    basisPL.validateParametersAndSetDefaults(valid_parameters_sublist);

    std::string basis_type = basisPL.get<std::string>("Basis Type");
    int basis_order = basisPL.get<int>("Basis Order");

    this->addDOF(field,basis_type,basis_order,default_integration_order);
  }
  
  this->addClosureModel("");

  this->setupDOFs();
}

PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(EquationSet_IC, EquationSet_IC)

// A user written factory that creates each equation set.  The key member here
// is buildEquationSet
class IC_EquationSetFactory : public EquationSetFactory {
public:

   Teuchos::RCP<EquationSet_TemplateManager<Traits> >
   buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		    const int& default_integration_order,
                    const CellData& cell_data,
		    const Teuchos::RCP<GlobalData>& global_data,
                    const bool build_transient_support) const
   {
      Teuchos::RCP<EquationSet_TemplateManager<Traits> > eq_set= 
         Teuchos::rcp(new EquationSet_TemplateManager<Traits>);
         
      bool found = false; // this is used by PANZER_BUILD_EQSET_OBJECTS
         
      // macro checks if(ies.name=="Poisson") then an EquationSet_Energy object is constructed
      PANZER_BUILD_EQSET_OBJECTS("IC", EquationSet_IC)
         
      // make sure your equation set has been found
      if(!found) {
	std::string msg = "Error - the \"Equation Set\" called \"" + params->get<std::string>("Type") +
                           "\" is not a valid equation set identifier. Please supply the correct factory.\n";
         TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
         
      return eq_set;
   }
    
};

} // end anonymous namespace

void 
setupControlInitialCondition(const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_ids_to_cell_topo,
                             const std::map<std::string,std::vector<ICFieldDescriptor> > & block_ids_to_fields,
                             WorksetContainer & wkstContainer,
                             const LinearObjFactory<Traits> & lof,
                             const ClosureModelFactory_TemplateManager<Traits> & cm_factory,
                             const Teuchos::ParameterList & ic_closure_models,
                             const Teuchos::ParameterList & user_data,
                             int workset_size,
                             double t0,
                             const Teuchos::RCP<Thyra::VectorBase<double> > & vec)
{
  std::vector<Teuchos::RCP<PhysicsBlock> > physics_blocks;
  buildICPhysicsBlocks(block_ids_to_cell_topo,block_ids_to_fields,workset_size,physics_blocks);

  std::map<std::string, Teuchos::RCP< PHX::FieldManager<Traits> > > phx_ic_field_managers;
  setupInitialConditionFieldManagers(wkstContainer,
                                               physics_blocks,
                                               cm_factory,
                                               ic_closure_models,
                                               lof,
                                               user_data,
                                               false,
                                               "initial_condition_control_test",
                                               phx_ic_field_managers);

  
  Teuchos::RCP<LinearObjContainer> loc  = lof.buildLinearObjContainer();
  Teuchos::rcp_dynamic_cast<ThyraObjContainer<double> >(loc)->set_x_th(vec);

  evaluateInitialCondition(wkstContainer, phx_ic_field_managers, loc, lof, t0);
}


void
buildICPhysicsBlocks(const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_ids_to_cell_topo,
                     const std::map<std::string,std::vector<ICFieldDescriptor> > & block_ids_to_fields,
                     int workset_size,
                     std::vector<Teuchos::RCP<PhysicsBlock> > & physics_blocks)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::map<std::string,std::string> block_ids_to_physics_ids;

  RCP<Teuchos::ParameterList> ipb = rcp(new Teuchos::ParameterList);

  for(auto itr=block_ids_to_cell_topo.begin();itr!=block_ids_to_cell_topo.end();++itr) {
    std::string eblock                 = itr->first;
    RCP<const shards::CellTopology> ct = itr->second;

    // get the field descriptor vector, check to make sure block is represented
    auto fds_itr = block_ids_to_fields.find(eblock);
    TEUCHOS_ASSERT(fds_itr!=block_ids_to_fields.end());

    const std::vector<ICFieldDescriptor> & fd_vec = fds_itr->second;

    std::string physics_id = "ic_"+eblock;
    block_ids_to_physics_ids[eblock] = physics_id;

    // start building a physics block named "ic_"+eblock, with one anonymous list
    Teuchos::ParameterList & physics_block = ipb->sublist(physics_id).sublist("");
    physics_block.set("Type","IC");  // point IC type

    for(std::size_t i=0;i<fd_vec.size();i++) {
      const ICFieldDescriptor & fd = fd_vec[i];

      // number the lists, these should be anonymous
      physics_block.sublist(fd.fieldName).set("Basis Type", fd.basisName);
      physics_block.sublist(fd.fieldName).set("Basis Order",fd.basisOrder);
    }
  }

  RCP<EquationSetFactory> eqset_factory = Teuchos::rcp(new IC_EquationSetFactory);

  RCP<GlobalData> gd = createGlobalData();
  buildPhysicsBlocks(block_ids_to_physics_ids,
                     block_ids_to_cell_topo,
                     ipb,1,workset_size,eqset_factory,gd,false,physics_blocks);
}

} // end namespace panzer
