#include "simpleStratimikosSolve.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"  // Contains create_MultiVector
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"  // Contains LinearOpWithSolve ?

/*
  simpleStratimikosSolve performs x = A \ b; 
  and returns 0 on success
 */

int simpleStratimikosSolve( 
   Epetra_CrsMatrix const& epetra_A,  // non-persisting, non-changeable
   Epetra_MultiVector const& epetra_B,  // non-persisting, non-changeable
   Epetra_MultiVector *epetra_X,  // non-persisting, changeable
   Teuchos::ParameterList *paramList  // non-persisting, changeable
   )
  {

    using Teuchos::RCP;
    using Teuchos::rcp;


    //
    // C) The "Glue" code that takes Epetra objects and wraps them as Thyra
    // objects
    //
    // This next set of code wraps the Epetra objects that define the linear
    // system to be solved as Thyra objects so that they can be passed to the
    // linear solver.
    //
 
    // Create RCPs that will be used to hold the Thyra wrappers

    typedef RCP<const Thyra::LinearOpBase<double> > LinearOpPtr;
    typedef RCP<Thyra::MutliVectorBase<double> > MultiVectorPtr;

    LinearOpPtr A = Thyra::epetraLinearOp( epetra_A );
    VectorPtr   X = Thyra::create_Vector( rcp(epetra_X,false), A->domain() );
    VectorPtr   B = Thyra::create_Vector( rcp(&epetra_B,false), A->range() );
 
    // Note that above Thyra is only interacted with in the most trival of
    // ways.  For most users, Thyra will only be seen as a thin wrapper that
    // they need know little about in order to wrap their objects in order to
    // pass them to Thyra-enabled solvers.

     
    //
    // D) Thyra-specific code for solving the linear system
    //
    // Note that this code has no mention of any concrete implementation and
    // therefore can be used in any use case.
    //
 
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    // Set the parameter list
    linearSolverBuilder.setParameterList( rcp(paramList,false) ) ; 
    // Create a linear solver factory given information read from the
    // parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");
    // Setup the verbosity level
    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);
    // Create a linear solver based on the forward operator A
    RCP<Thyra::LinearOpWithSolveBase<double> >
      lows = Thyra::linearOpWithSolve(*lowsFactory,A);
      //      lows = Thyra::linearOpWithSolve(*lowsFactory,rcp(&A,false));
    // Solve the linear system (note: the initial guess in 'X' is critical)
    Thyra::SolveStatus<double>
      status = Thyra::solve(*lows,Thyra::NOTRANS,*B,&*X);
    // Write the linear solver parameters after they were read
    linearSolverBuilder.writeParamsFile(*lowsFactory);


#if 0
    std::cout << __FILE__ << "::" << __LINE__ << 
      " paramlist  = " << *(linearSolverBuilder.getParameterList( ))   << std::endl ; 
#endif

    return (status.solveStatus!=Thyra::SOLVE_STATUS_CONVERGED);
 
  }
