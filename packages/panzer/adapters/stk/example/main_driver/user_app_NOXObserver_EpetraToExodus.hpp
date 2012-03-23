// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef USER_APP_NOX_OBSERVER_EPETRA_TO_EXODUS_HPP
#define USER_APP_NOX_OBSERVER_EPETRA_TO_EXODUS_HPP

#include "NOX_Abstract_PrePostOperator.H"
#include "Teuchos_RCP.hpp"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseAggregator_SolutionWriter.hpp"

#include "NOX_Epetra_Vector.H"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"

namespace user_app {
  
  class NOXObserver_EpetraToExodus : public NOX::Abstract::PrePostOperator {
    
  public:
    
    NOXObserver_EpetraToExodus(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			       const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
			       const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof,
                               const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & response_library) :
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof),
      m_response_library(response_library)
    { 
      // register solution writer response aggregator with this library
      // this is an "Action" only response aggregator
      panzer::ResponseAggregator_Manager<panzer::Traits> & aggMngr = m_response_library->getAggregatorManager();
      panzer_stk::ResponseAggregator_SolutionWriter_Builder builder(mesh);
      builder.setLinearObjFactory(aggMngr.getLinearObjFactory());
      builder.setGlobalIndexer(aggMngr.getGlobalIndexer());
      aggMngr.defineAggregatorTypeFromBuilder("Solution Writer",builder);

      // require a particular "Solution Writer" response
      panzer::ResponseId rid("Main Field Output","Solution Writer");
      std::list<std::string> eTypes;
      eTypes.push_back("Residual");
      #ifdef HAVE_STOKHOS
         eTypes.push_back("SGResidual");
      #endif

      std::list<std::string> eBlocks;
      {
         // get all element blocks and add them to the list
         std::vector<std::string> eBlockNames;
         mesh->getElementBlockNames(eBlockNames);
         for(std::size_t i=0;i<eBlockNames.size();i++)
            eBlocks.push_back(eBlockNames[i]);
      }

      // reserve response guranteeing that we can evaluate it (assuming things are done correctly elsewhere)
      response_library->reserveLabeledBlockAggregatedVolumeResponse("Main Field Output",rid,eBlocks,eTypes);
    }
      
    void runPreIterate(const NOX::Solver::Generic& solver)
    {

    }
    
    void runPostIterate(const NOX::Solver::Generic& solver)
    {

    }
    
    void runPreSolve(const NOX::Solver::Generic& solver)
    {

    }
    
    void runPostSolve(const NOX::Solver::Generic& solver)
    {
      const NOX::Abstract::Vector& x = solver.getSolutionGroup().getX();
      const NOX::Thyra::Vector* n_th_x = dynamic_cast<const NOX::Thyra::Vector*>(&x);
      TEUCHOS_TEST_FOR_EXCEPTION(n_th_x == NULL, std::runtime_error, "Failed to dynamic_cast to NOX::Thyra::Vector!")
      const ::Thyra::VectorBase<double>& th_x = n_th_x->getThyraVector(); 

      Teuchos::RCP<const Epetra_Vector> ep_x = Thyra::get_Epetra_Vector(*(m_lof->getMap()), Teuchos::rcp(&th_x, false));
      Teuchos::MpiComm<int> comm(Teuchos::opaqueWrapper(dynamic_cast<const Epetra_MpiComm &>(ep_x->Comm()).Comm()));

      // initialize the assembly container
      panzer::AssemblyEngineInArgs ae_inargs;
      ae_inargs.container_ = m_lof->buildLinearObjContainer();
      ae_inargs.ghostedContainer_ = m_lof->buildGhostedLinearObjContainer();
      ae_inargs.alpha = 0.0;
      ae_inargs.beta = 1.0;
      ae_inargs.evaluate_transient_terms = false;

      // initialize the ghosted container
      m_lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

      const Teuchos::RCP<panzer::EpetraLinearObjContainer> epGlobalContainer
         = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs.container_,true);
      epGlobalContainer->x = Teuchos::rcp_const_cast<Epetra_Vector>(ep_x);

      // do import
      m_lof->globalToGhostContainer(*ae_inargs.container_,*ae_inargs.ghostedContainer_,panzer::LinearObjContainer::X);

      // fill STK mesh objects
      m_response_library->evaluateVolumeFieldManagers<panzer::Traits::Residual>(ae_inargs,comm);
      
      // write to disk
      m_mesh->writeToExodus(0.0);
    }
    
  protected:

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > m_dof_manager;
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > m_lof;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;

  };
}

#endif
