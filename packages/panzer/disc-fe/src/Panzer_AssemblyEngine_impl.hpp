// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ASSEMBLY_ENGINE_IMPL_HPP
#define PANZER_ASSEMBLY_ENGINE_IMPL_HPP

#include "Phalanx_FieldManager.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"
#include <sstream>

//===========================================================================
//===========================================================================
template <typename EvalT>
panzer::AssemblyEngine<EvalT>::
AssemblyEngine(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof)
  : m_field_manager_builder(fmb), m_lin_obj_factory(lof), countersInitialized_(false)
{ 

}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluate(const panzer::AssemblyEngineInArgs& in, const EvaluationFlags flags)
{
  typedef LinearObjContainer LOC;

  // make sure this container gets a dirichlet adjustment
  in.ghostedContainer_->setRequiresDirichletAdjustment(true);

  GlobalEvaluationDataContainer gedc;

  if ( flags.getValue() & EvaluationFlags::Initialize ) {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::evaluate_gather("+PHX::print<EvalT>()+")", eval_gather);

    in.fillGlobalEvaluationDataContainer(gedc);
    gedc.initialize(); // make sure all ghosted data is ready to go
    gedc.globalToGhost(LOC::X | LOC::DxDt);

    // Push solution, x and dxdt into ghosted domain
    m_lin_obj_factory->globalToGhostContainer(*in.container_,*in.ghostedContainer_,LOC::X | LOC::DxDt);
    m_lin_obj_factory->beginFill(*in.ghostedContainer_);
  }

  // *********************
  // Volumetric fill
  // *********************
  if ( flags.getValue() & EvaluationFlags::VolumetricFill) {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::evaluate_volume("+PHX::print<EvalT>()+")", eval_vol);
    this->evaluateVolume(in);
  }

  // *********************
  // BC fill
  // *********************
  // NOTE: We have to split neumann and dirichlet bcs since dirichlet
  // bcs overwrite equations where neumann sum into equations.  Make
  // sure all neumann are done before dirichlet.

  if ( flags.getValue() & EvaluationFlags::BoundaryFill) {
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::evaluate_neumannbcs("+PHX::print<EvalT>()+")",eval_neumannbcs);
      this->evaluateNeumannBCs(in);
    }

    {
      PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::evaluate_interfacebcs("+PHX::print<EvalT>()+")",eval_interfacebcs);
      this->evaluateInterfaceBCs(in);
    }

    // Dirchlet conditions require a global matrix
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::evaluate_dirichletbcs("+PHX::print<EvalT>()+")",eval_dirichletbcs);
      this->evaluateDirichletBCs(in);
    }
  }

  if ( flags.getValue() & EvaluationFlags::Scatter) {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::evaluate_scatter("+PHX::print<EvalT>()+")",eval_scatter);
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::lof->ghostToGlobalContainer("+PHX::print<EvalT>()+")",lof_gtgc);
      m_lin_obj_factory->ghostToGlobalContainer(*in.ghostedContainer_,*in.container_,LOC::F | LOC::Mat);
    }
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("panzer::AssemblyEngine::gedc.ghostToGlobal("+PHX::print<EvalT>()+")",gedc_gtg);
      m_lin_obj_factory->beginFill(*in.container_);
      gedc.ghostToGlobal(LOC::F | LOC::Mat);
      m_lin_obj_factory->endFill(*in.container_);
    }
    m_lin_obj_factory->endFill(*in.ghostedContainer_);
  }

  return;
}

//===========================================================================
//===========================================================================
template <typename EvalT>
Teuchos::RCP<panzer::LinearObjContainer> panzer::AssemblyEngine<EvalT>::
evaluateOnlyDirichletBCs(const panzer::AssemblyEngineInArgs& in)
{
  typedef LinearObjContainer LOC;

  // make sure this container gets a dirichlet adjustment
  in.ghostedContainer_->setRequiresDirichletAdjustment(true);

  GlobalEvaluationDataContainer gedc;
  in.fillGlobalEvaluationDataContainer(gedc);
  gedc.initialize(); // make sure all ghosted data is ready to go
  gedc.globalToGhost(LOC::X | LOC::DxDt);

  // Push solution, x and dxdt into ghosted domain
  m_lin_obj_factory->globalToGhostContainer(*in.container_,*in.ghostedContainer_,LOC::X | LOC::DxDt);
  m_lin_obj_factory->beginFill(*in.ghostedContainer_);

  // Dirchlet conditions require a global matrix
  Teuchos::RCP<LOC> counter = this->evaluateDirichletBCs(in);

  m_lin_obj_factory->ghostToGlobalContainer(*in.ghostedContainer_,*in.container_,LOC::F | LOC::Mat);

  m_lin_obj_factory->beginFill(*in.container_);
  gedc.ghostToGlobal(LOC::F | LOC::Mat);
  m_lin_obj_factory->endFill(*in.container_);

  m_lin_obj_factory->endFill(*in.ghostedContainer_);

  return counter;
}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluateVolume(const panzer::AssemblyEngineInArgs& in)
{
  const std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > &
    volume_field_managers = m_field_manager_builder->getVolumeFieldManagers();
  const std::vector<WorksetDescriptor> & wkstDesc = m_field_manager_builder->getVolumeWorksetDescriptors();

  Teuchos::RCP<panzer::WorksetContainer> wkstContainer = m_field_manager_builder->getWorksetContainer();

  panzer::Traits::PED ped;
  ped.gedc->addDataObject("Solution Gather Container",in.ghostedContainer_);
  ped.gedc->addDataObject("Residual Scatter Container",in.ghostedContainer_);
  ped.first_sensitivities_name  = in.first_sensitivities_name;
  ped.second_sensitivities_name = in.second_sensitivities_name;
  in.fillGlobalEvaluationDataContainer(*(ped.gedc));

  // Loop over volume field managers
  for (std::size_t block = 0; block < volume_field_managers.size(); ++block) {
    const WorksetDescriptor & wd = wkstDesc[block];
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = volume_field_managers[block];
    std::vector<panzer::Workset>& w = *wkstContainer->getWorksets(wd);

    fm->template preEvaluate<EvalT>(ped);

    // Loop over worksets in this element block
    for (std::size_t i = 0; i < w.size(); ++i) {
      panzer::Workset& workset = w[i];

      workset.alpha = in.alpha;
      workset.beta = in.beta;
      workset.time = in.time;
      workset.step_size = in.step_size;
      workset.stage_number = in.stage_number;
      workset.gather_seeds = in.gather_seeds;
      workset.evaluate_transient_terms = in.evaluate_transient_terms;


      fm->template evaluateFields<EvalT>(workset);
    }

    // double s = 0.;
    // double p = 0.;
    // fm->template analyzeGraph<EvalT>(s,p);
    // std::cout << "Analyze Graph: " << PHX::print<EvalT>() << ",b=" << block << ", s=" << s << ", p=" << p << std::endl;

    fm->template postEvaluate<EvalT>(NULL);
  }
}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluateNeumannBCs(const panzer::AssemblyEngineInArgs& in)
{
  this->evaluateBCs(panzer::BCT_Neumann, in);
}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluateInterfaceBCs(const panzer::AssemblyEngineInArgs& in)
{
  this->evaluateBCs(panzer::BCT_Interface, in);
}

//===========================================================================
//===========================================================================
template <typename EvalT>
Teuchos::RCP<panzer::LinearObjContainer> panzer::AssemblyEngine<EvalT>::
evaluateDirichletBCs(const panzer::AssemblyEngineInArgs& in)
{
  typedef LinearObjContainer LOC;

  if(!countersInitialized_) {
    localCounter_ = m_lin_obj_factory->buildPrimitiveGhostedLinearObjContainer();
    globalCounter_ = m_lin_obj_factory->buildPrimitiveLinearObjContainer();
    summedGhostedCounter_ = m_lin_obj_factory->buildPrimitiveGhostedLinearObjContainer();
    countersInitialized_ = true;
 
    m_lin_obj_factory->initializeGhostedContainer(LinearObjContainer::F,*localCounter_); // store counter in F
    m_lin_obj_factory->initializeContainer(       LinearObjContainer::F,*globalCounter_); // store counter in X
    m_lin_obj_factory->initializeGhostedContainer(LinearObjContainer::F,*summedGhostedCounter_); // store counter in X
  }

  {
    localCounter_->initialize();
    summedGhostedCounter_->initialize();
    globalCounter_->initialize();
  }

  // apply dirichlet conditions, make sure to keep track of the local counter
  this->evaluateBCs(panzer::BCT_Dirichlet, in,localCounter_);

  // do communication to build summed ghosted counter for dirichlet conditions
  {
     m_lin_obj_factory->ghostToGlobalContainer(*localCounter_,*globalCounter_,LOC::F);
        // Here we do the reduction across all processors so that the number of times
        // a dirichlet condition is applied is summed into the global counter

     m_lin_obj_factory->globalToGhostContainer(*globalCounter_,*summedGhostedCounter_,LOC::F);
        // finally we move the summed global vector into a local ghosted vector
        // so that the dirichlet conditions can be applied to both the ghosted
        // right hand side and the ghosted matrix
  }

  panzer::GlobalEvaluationDataContainer gedc;
  gedc.addDataObject("Residual Scatter Container",in.ghostedContainer_);
  in.fillGlobalEvaluationDataContainer(gedc);

  // adjust ghosted system for boundary conditions
  for(GlobalEvaluationDataContainer::iterator itr=gedc.begin();itr!=gedc.end();itr++) {

    if(itr->second->requiresDirichletAdjustment()) {
      Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LinearObjContainer>(itr->second);
      if(loc!=Teuchos::null) {
        m_lin_obj_factory->adjustForDirichletConditions(*localCounter_,*summedGhostedCounter_,*loc);
      }
      else {
        // it was not a linear object container, so if you want an adjustment it better be a GED_BCAdjustment object
        Teuchos::RCP<GlobalEvaluationData_BCAdjustment> bc_adjust = Teuchos::rcp_dynamic_cast<GlobalEvaluationData_BCAdjustment>(itr->second,true);
        bc_adjust->adjustForDirichletConditions(*localCounter_,*summedGhostedCounter_);
      }
    }
  }

  return globalCounter_;
}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluateBCs(const panzer::BCType bc_type,
            const panzer::AssemblyEngineInArgs& in,
            const Teuchos::RCP<LinearObjContainer> preEval_loc)
{
  Teuchos::RCP<panzer::WorksetContainer> wkstContainer = m_field_manager_builder->getWorksetContainer();

  panzer::Traits::PED ped;
  ped.gedc->addDataObject("Dirichlet Counter",preEval_loc);
  ped.gedc->addDataObject("Solution Gather Container",in.ghostedContainer_);
  ped.gedc->addDataObject("Residual Scatter Container",in.ghostedContainer_);
  ped.first_sensitivities_name  = in.first_sensitivities_name;
  ped.second_sensitivities_name = in.second_sensitivities_name;
  in.fillGlobalEvaluationDataContainer(*(ped.gedc));

  // this helps work around issues when constructing a mass
  // matrix using an evaluation of only the transient terms.
  // In particular, the terms associated with the dirichlet
  // conditions.
  double betaValue = in.beta; // default to the passed in beta
  if(bc_type==panzer::BCT_Dirichlet && in.apply_dirichlet_beta) {
    betaValue = in.dirichlet_beta;
  }

  {
    const std::map<panzer::BC, 
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC>& bc_field_managers = 
      m_field_manager_builder->getBCFieldManagers();
  
    // Must do all neumann before all dirichlet so we need a double loop
    // here over all bcs
    typedef typename std::map<panzer::BC, 
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC>::const_iterator bcfm_it_type;

    // loop over bcs
    for (bcfm_it_type bcfm_it = bc_field_managers.begin(); 
         bcfm_it != bc_field_managers.end(); ++bcfm_it) {
      
      const panzer::BC& bc = bcfm_it->first;
      const std::map<unsigned,PHX::FieldManager<panzer::Traits> > bc_fm = 
        bcfm_it->second;
   
      panzer::WorksetDescriptor desc = panzer::bcDescriptor(bc);
      Teuchos::RCP<const std::map<unsigned,panzer::Workset> > bc_wkst_ptr = wkstContainer->getSideWorksets(desc);
      TEUCHOS_TEST_FOR_EXCEPTION(bc_wkst_ptr == Teuchos::null, std::logic_error,
                         "Failed to find corresponding bc workset!");
      const std::map<unsigned,panzer::Workset>& bc_wkst = *bc_wkst_ptr;

      // Only process bcs of the appropriate type (neumann or dirichlet)
      if (bc.bcType() == bc_type) {
        std::ostringstream timerName;
        timerName << "panzer::AssemblyEngine::evaluateBCs: " << bc.identifier();
#ifdef PANZER_TEUCHOS_TIME_MONITOR
        auto timer1 = Teuchos::TimeMonitor::getNewTimer(timerName.str());
        Teuchos::TimeMonitor tm1(*timer1);
#endif

        // Loop over local faces
        for (std::map<unsigned,PHX::FieldManager<panzer::Traits> >::const_iterator side = bc_fm.begin(); side != bc_fm.end(); ++side) {
          std::ostringstream timerSideName;
          timerSideName << "panzer::AssemblyEngine::evaluateBCs: " << bc.identifier() << ", side=" << side->first;
#ifdef PANZER_TEUCHOS_TIME_MONITOR
        auto timer2 = Teuchos::TimeMonitor::getNewTimer(timerSideName.str());
        Teuchos::TimeMonitor tm2(*timer2);
#endif

          // extract field manager for this side  
          unsigned local_side_index = side->first;
          PHX::FieldManager<panzer::Traits>& local_side_fm = 
            const_cast<PHX::FieldManager<panzer::Traits>& >(side->second);
          
          // extract workset for this side: only one workset per face
          std::map<unsigned,panzer::Workset>::const_iterator wkst_it = 
            bc_wkst.find(local_side_index);
          
          TEUCHOS_TEST_FOR_EXCEPTION(wkst_it == bc_wkst.end(), std::logic_error,
                             "Failed to find corresponding bc workset side!");
          
          panzer::Workset& workset = 
            const_cast<panzer::Workset&>(wkst_it->second); 

          // run prevaluate
          local_side_fm.template preEvaluate<EvalT>(ped);

          // build and evaluate fields for the workset: only one workset per face
          workset.alpha = in.alpha;
          workset.beta = betaValue;
          workset.time = in.time;
          workset.gather_seeds = in.gather_seeds;
          workset.evaluate_transient_terms = in.evaluate_transient_terms;
          
          local_side_fm.template evaluateFields<EvalT>(workset);

          // run postevaluate for consistency
          local_side_fm.template postEvaluate<EvalT>(NULL);
          
        }
      }
    } 
  }

}

#endif
