#ifndef PANZER_STK_RYTHMOS_OBSERVER_HPP
#define PANZER_STK_RYTHMOS_OBSERVER_HPP

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeRange.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "Panzer_STK_Utilities.hpp"

namespace panzer_stk {

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

      write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(m_dof_manager),*m_mesh,
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
