// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This driver reads a problem from a Harwell-Boeing (HB) file into an 
// Epetra_CrsMatrix.  This matrix is then converted into a Thyra linear operator
// through the Thyra-Epetra adapters.
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero. 
//
//

// Thyra includes
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

// Epetra includes
#include "Epetra_SerialComm.h"
#include "EpetraExt_readEpetraLinearSystem.h"

// Ifpack includes
#ifdef HAVE_BELOS_IFPACK
  #include "Thyra_IfpackPreconditionerFactory.hpp"
#endif

// Teuchos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[])
{
  //
  // Get a default output stream from the Teuchos::VerboseObjectBase
  //
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //
  // Set the parameters for the Belos LOWS Factory and create a parameter list.
  //
  int             blockSize              = 2;
  int             maxIterations          = 400;
  int             maxRestarts            = 25;
  int             gmresKrylovLength      = 25;
  int             outputFrequency        = 1;
  bool            outputMaxResOnly       = true;
  double          maxResid               = 1e-6;
  
  Teuchos::RCP<Teuchos::ParameterList>
    belosLOWSFPL = Teuchos::rcp( new Teuchos::ParameterList() );
  
  belosLOWSFPL->set("Solver Type","Block GMRES");
  
  Teuchos::ParameterList& belosLOWSFPL_solver = 
    belosLOWSFPL->sublist("Solver Types");

  Teuchos::ParameterList& belosLOWSFPL_gmres = 
    belosLOWSFPL_solver.sublist("Block GMRES");

  belosLOWSFPL_gmres.set("Maximum Iterations",int(maxIterations));
  belosLOWSFPL_gmres.set("Convergence Tolerance",double(maxResid));
  belosLOWSFPL_gmres.set("Maximum Restarts",int(maxRestarts));
  belosLOWSFPL_gmres.set("Block Size",int(blockSize));
  belosLOWSFPL_gmres.set("Num Blocks",int(gmresKrylovLength));
  belosLOWSFPL_gmres.set("Output Frequency",int(outputFrequency));
  belosLOWSFPL_gmres.set("Show Maximum Residual Norm Only",bool(outputMaxResOnly));

#ifdef HAVE_BELOS_IFPACK
  //
  // Set the parameters for the Ifpack Preconditioner Factory and create parameter list
  //
  Teuchos::ParameterList &ifpackPFSL = belosLOWSFPL->sublist("IfpackPreconditionerFactory");
  
  ifpackPFSL.set("Overlap",int(2));
  ifpackPFSL.set("Prec Type","ILUT");
#endif

  // Whether the linear solver succeeded.
  // (this will be set during the residual check at the end)
  bool success = true;

  // Number of random right-hand sides we will be solving for.
  int numRhs = 5;

  // Name of input matrix file
  std::string matrixFile = "orsirr1.hb";

  // Read in the matrix file (can be *.mtx, *.hb, etc.)
  Epetra_SerialComm comm;
  Teuchos::RCP<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );
  
  // Create a Thyra linear operator (A) using the Epetra_CrsMatrix (epetra_A).
  Teuchos::RCP<const Thyra::LinearOpBase<double> > 
    A = Thyra::epetraLinearOp(epetra_A);

  // Get the domain space for the Thyra linear operator 
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domain = A->domain();

  // Create the Belos LOWS factory.
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    belosLOWSFactory = Teuchos::rcp(new Thyra::BelosLinearOpWithSolveFactory<double>());

#ifdef HAVE_BELOS_IFPACK

  // Set the preconditioner factory for the LOWS factory.
  belosLOWSFactory->setPreconditionerFactory(
					     Teuchos::rcp(new Thyra::IfpackPreconditionerFactory())
					     ,"IfpackPreconditionerFactory"
					     );
#endif

  // Set the parameter list to specify the behavior of the factory.
  belosLOWSFactory->setParameterList( belosLOWSFPL );

  // Set the output stream and the verbosity level (prints to std::cout by defualt)
  belosLOWSFactory->setVerbLevel(Teuchos::VERB_LOW);

  // Create a BelosLinearOpWithSolve object from the Belos LOWS factory.
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> >
    nsA = belosLOWSFactory->createOp();

  // Initialize the BelosLinearOpWithSolve object with the Thyra linear operator.
  Thyra::initializeOp<double>( *belosLOWSFactory, A, &*nsA );

  // Create a right-hand side with numRhs vectors in it and randomize it.
  Teuchos::RCP< Thyra::MultiVectorBase<double> > 
    b = Thyra::createMembers(domain, numRhs);
  Thyra::seed_randomize<double>(0);
  Thyra::randomize(-1.0, 1.0, &*b);

  // Create an initial std::vector with numRhs vectors in it and initialize it to zero.
  Teuchos::RCP< Thyra::MultiVectorBase<double> >
    x = Thyra::createMembers(domain, numRhs);
  Thyra::assign(&*x, 0.0);

  // Perform solve using the linear operator to get the approximate solution of Ax=b,
  // where b is the right-hand side and x is the left-hand side.
  Thyra::SolveStatus<double> solveStatus;
  solveStatus = Thyra::solve( *nsA, Thyra::NONCONJ_ELE, *b, &*x );

  // Print out status of solve.
  *out << "\nBelos LOWS Status: "<< solveStatus << std::endl;

  //
  // Compute residual and double check convergence.
  //
  std::vector<double> norm_b(numRhs), norm_res(numRhs);
  Teuchos::RCP< Thyra::MultiVectorBase<double> >
    y = Thyra::createMembers(domain, numRhs);

  // Compute the column norms of the right-hand side b.
  Thyra::norms_2( *b, &norm_b[0] );

  // Compute y=A*x, where x is the solution from the linear solver.
  A->apply( Thyra::NONCONJ_ELE, *x, &*y );
  
  // Compute A*x-b = y-b
  Thyra::update( -1.0, *b, &*y );

  // Compute the column norms of A*x-b.
  Thyra::norms_2( *y, &norm_res[0] );

  // Print out the final relative residual norms.
  double rel_res = 0.0;
  *out << "Final relative residual norms" << std::endl;  
  for (int i=0; i<numRhs; ++i) {
    rel_res = norm_res[i]/norm_b[i];
    if (rel_res > maxResid)
      success = false;
    *out << "RHS " << i+1 << " : " 
         << std::setw(16) << std::right << rel_res << std::endl;
  }

  return ( success ? 0 : 1 );
}
