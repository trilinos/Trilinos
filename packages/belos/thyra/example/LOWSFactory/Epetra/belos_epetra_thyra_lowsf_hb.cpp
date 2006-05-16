
// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file into an 
// Epetra_CrsMatrix.  This matrix is then converted into a Thyra linear operator
// through the Thyra-Epetra adapters.
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero. 
//
//

// Thyra includes
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
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
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //
  // Set the parameters for the Belos LOWS Factory and create a parameter list.
  //
  int             blockSize              = 2;
  int             maxIterations          = 400;
  int             maxRestarts            = 25;
  int             gmresKrylovLength      = 25;
  int             outputFrequency        = 10;
  bool            outputMaxResOnly       = true;
  double          maxResid               = 1e-6;
  
  Teuchos::RefCountPtr<Teuchos::ParameterList>
    belosLOWSFPL = Teuchos::rcp( new Teuchos::ParameterList() );
  
  belosLOWSFPL->set("Solver Type","GMRES");
  belosLOWSFPL->set("Max Iters",int(maxIterations));
  belosLOWSFPL->set("Default Rel Res Norm",double(maxResid));
  belosLOWSFPL->set("Max Restarts",int(maxRestarts));
  belosLOWSFPL->set("Block Size",int(blockSize));
  belosLOWSFPL->sublist("GMRES").set("Max Number of Krylov Vectors",int(gmresKrylovLength*blockSize));
  belosLOWSFPL->sublist("GMRES").set("Variant","Standard");
  Teuchos::ParameterList &outputterSL = belosLOWSFPL->sublist("Outputter");
  outputterSL.set("Output Frequency",int(outputFrequency));
  outputterSL.set("Output Max Res Only",bool(outputMaxResOnly));

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
  Teuchos::RefCountPtr<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );
  
  // Create a Thyra linear operator (A) using the Epetra_CrsMatrix (epetra_A).
  Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > 
    A = Teuchos::rcp(new Thyra::EpetraLinearOp(epetra_A));

  // Get the domain space for the Thyra linear operator 
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> > domain = A->domain();

  // Create the Belos LOWS factory.
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
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

  // Create a BelosLinearOpWithSolve object from the Belos LOWS factory.
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<double> >
    nsA = belosLOWSFactory->createOp();

  // Initialize the BelosLinearOpWithSolve object with the Thyra linear operator.
  belosLOWSFactory->initializeOp( A, &*nsA );

  // Create a right-hand side with numRhs vectors in it and randomize it.
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<double> > 
    b = Thyra::createMembers(domain, numRhs);
  Thyra::seed_randomize<double>(0);
  Thyra::randomize(-1.0, 1.0, &*b);

  // Create an initial vector with numRhs vectors in it and initialize it to zero.
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<double> >
    x = Thyra::createMembers(domain, numRhs);
  Thyra::assign(&*x, 0.0);

  // Perform solve using the linear operator to get the approximate solution of Ax=b,
  // where b is the right-hand side and x is the left-hand side.
  Thyra::SolveStatus<double> solveStatus;
  solveStatus = Thyra::solve( *nsA, Thyra::NONCONJ_ELE, *b, &*x );

  // Print out status of solve.
  *out << "\nBelos LOWS Status: "<< solveStatus << endl;

  //
  // Compute residual and double check convergence.
  //
  std::vector<double> norm_b(numRhs), norm_res(numRhs);
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<double> >
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
  *out << "Final relative residual norms" << endl;  
  for (int i=0; i<numRhs; ++i) {
    rel_res = norm_res[i]/norm_b[i];
    if (rel_res > maxResid)
      success = false;
    *out << "RHS " << i+1 << " : " 
	 << std::setw(16) << std::right << rel_res << endl;
  }

  return ( success ? 0 : 1 );
}
