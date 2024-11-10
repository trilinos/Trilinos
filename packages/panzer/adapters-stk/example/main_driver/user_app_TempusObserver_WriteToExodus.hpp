// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_TEMPUS_OBSERVER_WRITE_TO_EXODUS_HPP
#define USER_APP_TEMPUS_OBSERVER_WRITE_TO_EXODUS_HPP

#include "Tempus_Integrator.hpp"
#include "Tempus_IntegratorObserver.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_STK_Utilities.hpp"

namespace user_app {

  class TempusObserver_WriteToExodus : public Tempus::IntegratorObserver<double> {

  public:

    TempusObserver_WriteToExodus(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                                 const Teuchos::RCP<const panzer::GlobalIndexer>& dof_manager,
                                 const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof,
                                 const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & response_library) :
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof),
      m_response_library(response_library)
    {
      // get all element blocks and add them to the list
      std::vector<std::string> eBlocks;
      mesh->getElementBlockNames(eBlocks);

      panzer_stk::RespFactorySolnWriter_Builder builder;
      builder.mesh = mesh;
      m_response_library->addResponse("Main Field Output",eBlocks,builder);
    }

    void observeStartIntegrator(const Tempus::Integrator<double>& integrator) override {this->writeToExodus(integrator);};
    void observeStartTimeStep(const Tempus::Integrator<double>& ) override {};
    void observeNextTimeStep(const Tempus::Integrator<double>& ) override {};
    void observeBeforeTakeStep(const Tempus::Integrator<double>& ) override {};
    void observeAfterTakeStep(const Tempus::Integrator<double>& integrator) override {this->writeToExodus(integrator);}
    void observeAfterCheckTimeStep(const Tempus::Integrator<double>& ) override {};
    void observeEndTimeStep(const Tempus::Integrator<double>& ) override {};
    void observeEndIntegrator(const Tempus::Integrator<double>& ) override {};

    void writeToExodus(const Tempus::Integrator<double>& integrator)
    {

      Teuchos::RCP<const Thyra::VectorBase<double> > solution = integrator.getSolutionHistory()->getStateTimeIndexN()->getX();

      // initialize the assembly container
      panzer::AssemblyEngineInArgs ae_inargs;
      ae_inargs.container_ = m_lof->buildLinearObjContainer();
      ae_inargs.ghostedContainer_ = m_lof->buildGhostedLinearObjContainer();
      ae_inargs.alpha = 0.0;
      ae_inargs.beta = 1.0;
      ae_inargs.evaluate_transient_terms = false;

      // initialize the ghosted container
      m_lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

      {
         // initialize the x vector
         const Teuchos::RCP<panzer::ThyraObjContainer<double> > thyraContainer
            = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.container_,true);
         thyraContainer->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(solution));
      }

      m_response_library->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
      m_response_library->evaluate<panzer::Traits::Residual>(ae_inargs);
      m_mesh->writeToExodus(integrator.getSolutionHistory()->getCurrentTime());
    }

  protected:

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<const panzer::GlobalIndexer> m_dof_manager;
    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > m_lof;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;
  };

}

#endif
