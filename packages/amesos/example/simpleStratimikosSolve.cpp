#include "simpleStratimikosSolve.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"  // Contains create_MultiVector
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"  // Contains LinearOpWithSolve ?

/*
  simpleStratimikosSolve performs x = A \ b; 
  and returns 0 on success
 */

int simpleStratimikosSolve( 
   Epetra_CrsMatrix                          const& epetra_A,  // non-persisting, non-changeable
   Epetra_MultiVector                        const& epetra_B,  // non-persisting, non-changeable
   Epetra_MultiVector                             * epetra_X,  // non-persisting, changeable
   Teuchos::ParameterList                         * paramList  // non-persisting, changeable
   )
  {

    using namespace Teuchos;
    //
    // C) The "Glue" code that takes Epetra objects and wraps them as Thyra
    // objects
    //
    // This next set of code wraps the Epetra objects that define the linear
    // system to be solved as Thyra objects so that they can be passed to the
    // linear solver.
    //
 
    // Create RCPs that will be used to hold the Thyra wrappers

    RefCountPtr<Epetra_Operator >  E_OP_A = rcp( const_cast<Epetra_Operator *>(dynamic_cast<const Epetra_Operator *>( &epetra_A )), false );

    RefCountPtr<const Thyra::LinearOpBase<double> >      A;
    RefCountPtr<Thyra::MultiVectorBase<double> >         X;
    RefCountPtr<const Thyra::MultiVectorBase<double> >   B;
 
    // Create the Thyra wrappers
    if(1) {
      // Create an RCP directly to the EpetraLinearOp so that we can access the
      // right range and domains spaces to use to create the wrappers for the
      // vector objects.
      RefCountPtr<const Thyra::EpetraLinearOp>
	_A = rcp(new Thyra::EpetraLinearOp(E_OP_A));
      // Create almost any old vector space for the domain space for the
      // multi-vectors (this helps in reuse).
      RefCountPtr<const Thyra::ScalarProdVectorSpaceBase<double> >
        mvDomain = rcp(new Thyra::DefaultSpmdVectorSpace<double>(epetra_B.NumVectors()));
      // Create Thyra wrappers for the multi-vector objects that will automatically
      // update the Epetra objects.
      B = Thyra::create_MultiVector(rcp(&epetra_B,false),_A->spmdRange(),mvDomain);
      X = Thyra::create_MultiVector(rcp(epetra_X,false),_A->spmdDomain(),mvDomain);
      // Set the RCP to the base linear operator interface
      A = _A;
    }
 
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
 
    Thyra::DefaultRealLinearSolverBuilder linearSolverBuilder;
    // Set the parameter list
    const Teuchos::RefCountPtr<Teuchos::ParameterList> RCP_PL = rcp( new Teuchos::ParameterList( *paramList ) );

    linearSolverBuilder.setParameterList( RCP_PL ) ; 
    // Create a linear solver factory given information read from the
    // parameter list.
    RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");
    // Setup the verbosity level
    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);
    // Create a linear solver based on the forward operator A
    RefCountPtr<Thyra::LinearOpWithSolveBase<double> >
      lows = Thyra::linearOpWithSolve(*lowsFactory,A);
      //      lows = Thyra::linearOpWithSolve(*lowsFactory,rcp(&A,false));
    // Solve the linear system (note: the initial guess in 'X' is critical)
    Thyra::SolveStatus<double>
      status = Thyra::solve(*lows,Thyra::NOTRANS,*B,&*X);
    // Write the linear solver parameters after they were read
    linearSolverBuilder.writeParamsFile(*lowsFactory);

#if 0
    cout << __FILE__ << "::" << __LINE__ << 
      " paramlist  = " << *(linearSolverBuilder.getParameterList( ))   << endl ; 
#endif

    *paramList = *(linearSolverBuilder.getParameterList( ));

    return (status.solveStatus!=Thyra::SOLVE_STATUS_CONVERGED);
 
  }
