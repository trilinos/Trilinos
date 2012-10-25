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

#ifndef PANZER_ASSEMBLY_ENGINE_IMPL_HPP
#define PANZER_ASSEMBLY_ENGINE_IMPL_HPP

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"

//===========================================================================
//===========================================================================
template <typename EvalT>
panzer::AssemblyEngine<EvalT>::
AssemblyEngine(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof)
  : m_field_manager_builder(fmb), m_lin_obj_factory(lof)
{ 

}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluate(const panzer::AssemblyEngineInArgs& in)
{
  typedef LinearObjContainer LOC;

  GlobalEvaluationDataContainer gedc;
  in.fillGlobalEvaluationDataContainer(gedc);
  gedc.initialize(); // make sure all ghosted data is ready to go
  gedc.globalToGhost(LOC::X | LOC::DxDt);

  // Push solution, x and dxdt into ghosted domain
  m_lin_obj_factory->globalToGhostContainer(*in.container_,*in.ghostedContainer_,LOC::X | LOC::DxDt);

  // *********************
  // Volumetric fill
  // *********************
  this->evaluateVolume(in);

  // *********************
  // BC fill
  // *********************
  // NOTE: We have to split neumann and dirichlet bcs since dirichlet
  // bcs overwrite equations where neumann sum into equations.  Make
  // sure all neumann are done before dirichlet.

  this->evaluateNeumannBCs(in);

  // Dirchlet conditions require a global matrix
  this->evaluateDirichletBCs(in);

  m_lin_obj_factory->ghostToGlobalContainer(*in.ghostedContainer_,*in.container_,LOC::F | LOC::Mat);

  gedc.ghostToGlobal(LOC::F | LOC::Mat);

  return;
}

//===========================================================================
//===========================================================================
template <typename EvalT>
void panzer::AssemblyEngine<EvalT>::
evaluateVolume(const panzer::AssemblyEngineInArgs& in)
{
  const std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > &
    volume_field_managers = m_field_manager_builder->getVolumeFieldManagers();
  const std::vector<std::string> & elmt_blk_names = m_field_manager_builder->getElementBlockNames();

  Teuchos::RCP<panzer::WorksetContainer> wkstContainer = m_field_manager_builder->getWorksetContainer();

  GlobalEvaluationDataContainer gedc;
  gedc.addDataObject("Solution Gather Container",in.ghostedContainer_);
  gedc.addDataObject("Residual Scatter Container",in.ghostedContainer_);
  in.fillGlobalEvaluationDataContainer(gedc);

  // Loop over volume field managers
  for (std::size_t block = 0; block < volume_field_managers.size(); ++block) {
    const std::string & blockId = elmt_blk_names[block];
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = volume_field_managers[block];
    std::vector<panzer::Workset>& w = *wkstContainer->getVolumeWorksets(blockId);

    fm->template preEvaluate<EvalT>(gedc);

    // Loop over worksets in this element block
    for (std::size_t i = 0; i < w.size(); ++i) {
      panzer::Workset& workset = w[i];

      workset.ghostedLinContainer = in.ghostedContainer_;
      workset.linContainer = in.container_;
      workset.alpha = in.alpha;
      workset.beta = in.beta;
      workset.time = in.time;
      workset.evaluate_transient_terms = in.evaluate_transient_terms;

      fm->template evaluateFields<EvalT>(workset);
    }

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
evaluateDirichletBCs(const panzer::AssemblyEngineInArgs& in)
{
  typedef LinearObjContainer LOC;

  // allocate a counter to keep track of where this processor set dirichlet boundary conditions
  Teuchos::RCP<LinearObjContainer> localCounter = m_lin_obj_factory->buildPrimitiveGhostedLinearObjContainer();
  m_lin_obj_factory->initializeGhostedContainer(LinearObjContainer::X,*localCounter); // store counter in X
  localCounter->initialize();
     // this has only an X vector. The evaluate BCs will add a one to each row
     // that has been set as a dirichlet condition on this processor

  // apply dirichlet conditions, make sure to keep track of the local counter
  this->evaluateBCs(panzer::BCT_Dirichlet, in,localCounter);

  Teuchos::RCP<LinearObjContainer> summedGhostedCounter = m_lin_obj_factory->buildPrimitiveGhostedLinearObjContainer();
  m_lin_obj_factory->initializeGhostedContainer(LinearObjContainer::X,*summedGhostedCounter); // store counter in X
  summedGhostedCounter->initialize();

  // do communication to build summed ghosted counter for dirichlet conditions
  {
     Teuchos::RCP<LinearObjContainer> globalCounter = m_lin_obj_factory->buildPrimitiveLinearObjContainer();
     m_lin_obj_factory->initializeContainer(LinearObjContainer::X,*globalCounter); // store counter in X
     globalCounter->initialize();
     m_lin_obj_factory->ghostToGlobalContainer(*localCounter,*globalCounter,LOC::X);
        // Here we do the reduction across all processors so that the number of times
        // a dirichlet condition is applied is summed into the global counter

     m_lin_obj_factory->globalToGhostContainer(*globalCounter,*summedGhostedCounter,LOC::X);
        // finally we move the summed global vector into a local ghosted vector
        // so that the dirichlet conditions can be applied to both the ghosted
        // right hand side and the ghosted matrix
  }

  // adjust ghosted system for boundary conditions
  m_lin_obj_factory->adjustForDirichletConditions(*localCounter,*summedGhostedCounter,*in.ghostedContainer_);
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

  panzer::GlobalEvaluationDataContainer gedc;

  gedc.addDataObject("Dirichlet Counter",preEval_loc);
  gedc.addDataObject("Solution Gather Container",in.ghostedContainer_);
  gedc.addDataObject("Residual Scatter Container",in.ghostedContainer_);
  in.fillGlobalEvaluationDataContainer(gedc);

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
   
      Teuchos::RCP<const std::map<unsigned,panzer::Workset> > bc_wkst_ptr = wkstContainer->getSideWorksets(bc);
      TEUCHOS_TEST_FOR_EXCEPTION(bc_wkst_ptr == Teuchos::null, std::logic_error,
			 "Failed to find corresponding bc workset!");
      const std::map<unsigned,panzer::Workset>& bc_wkst = *bc_wkst_ptr;

      // Only process bcs of the appropriate type (neumann or dirichlet)
      if (bc.bcType() == bc_type) {

	// Loop over local faces
	for (std::map<unsigned,PHX::FieldManager<panzer::Traits> >::const_iterator side = bc_fm.begin(); side != bc_fm.end(); ++side) {

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
          local_side_fm.template preEvaluate<EvalT>(gedc);

          // build and evaluate fields for the workset: only one workset per face
          workset.ghostedLinContainer = in.ghostedContainer_;
          workset.linContainer = in.container_;
	  workset.alpha = in.alpha;
	  workset.beta = in.beta;
	  workset.time = in.time;
	  
	  local_side_fm.template evaluateFields<EvalT>(workset);

          // run postevaluate for consistency
	  local_side_fm.template postEvaluate<EvalT>(NULL);
	  
	}
      }
    } 
  }

}

#endif
