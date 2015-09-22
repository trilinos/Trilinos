
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"

#include "Teuchos_UnitTestHarness.hpp"

namespace Thyra {

  /* This test is to make sure that stratimikos wrappers to belos
     correctly catch belos exceptions thrown and report back a failed
     solve.  If belos detects a nan in the linear solve, it will throw
     an exception.  This was causing the termination of the
     simulation.  We have codes where cutting the time step allows
     recovery so it makes sense that stratimikos should catch the
     thrown exception and report that the solve failed.  Then the time
     integrator can retake a time step.  This test makes sure the
     stratimikos interface catches the exception and reports a failed
     solve.
  */
  TEUCHOS_UNIT_TEST(belos, nan_handling)
  {
    
    Epetra_SerialComm comm;
    Teuchos::RCP<Epetra_CrsMatrix> epetra_A;
    std::string matrixFile = "FourByFour.mtx";
    EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

    Teuchos::RCP<const LinearOpBase<double> > A = epetraLinearOp(epetra_A);

    Teuchos::RCP<LinearOpWithSolveFactoryBase<double> >
      lowsFactory;
    {
      Teuchos::RCP<BelosLinearOpWithSolveFactory<double> >
        belosLowsFactory = Teuchos::rcp(new BelosLinearOpWithSolveFactory<double>());
      lowsFactory = belosLowsFactory;
    }

    Teuchos::ParameterList belosLOWSFPL;
    {
      belosLOWSFPL.set("Solver Type","Block GMRES");

      Teuchos::ParameterList& belosLOWSFPL_solver =
	belosLOWSFPL.sublist("Solver Types");
      
      Teuchos::ParameterList& belosLOWSFPL_gmres =
	belosLOWSFPL_solver.sublist("Block GMRES");
      
      belosLOWSFPL_gmres.set("Maximum Iterations",int(4));
      belosLOWSFPL_gmres.set("Convergence Tolerance",double(1.0e-4));
      belosLOWSFPL_gmres.set("Maximum Restarts",int(0));
      belosLOWSFPL_gmres.set("Block Size",int(1));
      belosLOWSFPL_gmres.set("Num Blocks",int(4));
      belosLOWSFPL_gmres.set("Output Frequency",int(1));
      belosLOWSFPL_gmres.set("Show Maximum Residual Norm Only",bool(false));
      
      lowsFactory->setParameterList(Teuchos::rcp(&belosLOWSFPL,false));
    }

    Teuchos::RCP<LinearOpWithSolveBase<double> > nsA = lowsFactory->createOp();
    Thyra::initializeOp<double>(*lowsFactory,  A, nsA.ptr());

    Teuchos::RCP<Thyra::VectorBase<double> > x = Thyra::createMember(A->domain());
    Teuchos::RCP<Thyra::VectorBase<double> > f = Thyra::createMember(A->range());

    Thyra::put_scalar(0.0,x.ptr());
    Thyra::put_scalar(1.0,f.ptr());
    
    // Insert a nan
    Thyra::set_ele(0,Teuchos::ScalarTraits<double>::nan(),x.ptr());

    Thyra::SolveStatus<double> status = Thyra::solve<double>(*nsA,Thyra::NOTRANS,*f,x.ptr());

    TEST_EQUALITY(status.solveStatus, Thyra::SOLVE_STATUS_UNCONVERGED);
  }
  
}
