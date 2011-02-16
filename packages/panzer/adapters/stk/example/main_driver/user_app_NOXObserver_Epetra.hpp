#ifndef USER_APP_NOX_OBSERVER_EPETRA_HPP
#define USER_APP_NOX_OBSERVER_EPETRA_HPP

#include "NOX_Abstract_PrePostOperator.H"
#include "Teuchos_RCP.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "NOX_Epetra_Vector.H"
#include "Epetra_Vector.h"

namespace user_app {
  
  class NOXObserver_Epetra : public NOX::Abstract::PrePostOperator {
    
  public:
    
    NOXObserver_Epetra(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
		       const RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
		       const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof) :
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof)
    { }
      
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
      TEST_FOR_EXCEPTION(n_th_x == NULL, std::runtime_error, "Failed to dynamic_cast to NOX::Thyra::Vector!")
      const ::Thyra::VectorBase<double>& th_x = n_th_x->getThyraVector(); 

      Teuchos::RCP<const Epetra_Vector> ep_x = Thyra::get_Epetra_Vector(*(m_lof->getMap()), Teuchos::rcp(&th_x, false));

      Epetra_Vector ghosted_solution(*(m_lof->getGhostedMap()));
      Teuchos::RCP<Epetra_Import> importer = m_lof->getGhostedImport();
      ghosted_solution.PutScalar(0.0);
      ghosted_solution.Import(*ep_x,*importer,Insert);

      write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(m_dof_manager),*m_mesh,
			  ghosted_solution);
      
      m_mesh->writeToExodus(0.0);
    }
    
  protected:

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > m_dof_manager;
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > m_lof;

  };
}

#endif
