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

#ifndef PANZER_STK_RYTHMOS_OBSERVER_HPP
#define PANZER_STK_RYTHMOS_OBSERVER_HPP

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeRange.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "Panzer_STK_Utilities.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"

namespace user_app {

  class RythmosObserver_Epetra : 
    public Rythmos::IntegrationObserverBase<double> {

  public:
    
    RythmosObserver_Epetra(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			   const RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
			   const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof) :
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof)
    { }
    
    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    cloneIntegrationObserver() const
    {
      return Teuchos::rcp(new RythmosObserver_Epetra(m_mesh, m_dof_manager, m_lof));
    }

    void resetIntegrationObserver(const Rythmos::TimeRange<double> &integrationTimeDomain)
    { }

    void observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
				  const Rythmos::StepControlInfo<double> &stepCtrlInfo,
				  const int timeStepIter)
    { 
      cout << "*************************ROGER in Time************************************"  << endl;
      cout << "time = " << stepper.getStepStatus().time << endl;
      cout << *(stepper.getStepStatus().solution);
      Teuchos::RCP<const Thyra::VectorBase<double> > solution = stepper.getStepStatus().solution;
      
      // Next few lines are inefficient, but we can revisit later
      Teuchos::RCP<const Epetra_Vector> ep_solution = Thyra::get_Epetra_Vector(*(m_lof->getMap()), solution);
      Epetra_Vector ghosted_solution(*(m_lof->getGhostedMap()));
      Teuchos::RCP<Epetra_Import> importer = m_lof->getGhostedImport();
      ghosted_solution.PutScalar(0.0);
      ghosted_solution.Import(*ep_solution,*importer,Insert);

      panzer_stk::write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(m_dof_manager),*m_mesh,
			  ghosted_solution);
      
      m_mesh->writeToExodus(stepper.getStepStatus().time);
    }
    
  protected:

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > m_dof_manager;
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > m_lof;

  };

}

#endif
