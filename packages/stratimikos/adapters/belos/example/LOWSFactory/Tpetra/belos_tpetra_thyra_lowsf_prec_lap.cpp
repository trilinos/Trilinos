
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
// This driver creates a 2D Laplacian operator as a Tpetra::CisMatrix
// and a preconditioner.  These matrices are then converted into Thyra 
// linear operators through the Thyra-Tpetra adapters.
//
// The right-hand sides are all random vectors.
// The initial guesses are all set to zero. 
//
//

// Thyra includes
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

// Tpetra includes
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"

// Teuchos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
#  include "Tpetra_MpiPlatform.hpp"
#else
#  include "Tpetra_SerialPlatform.hpp"
#endif

int main(int argc, char* argv[])
{
  //
  // Get a default output stream from the Teuchos::VerboseObjectBase
  //
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  //
  // Create the scalar-type typedefs
  //

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;  // Scalar-type typedef
#elif HAVE_COMPLEX_H
  typedef std::complex<double> ST;     // Scalar-type typedef
#else
  typedef double ST;                // Scalar-type typedef
#endif
  typedef Teuchos::ScalarTraits<ST>::magnitudeType MT;  // Magnitude-type typdef
  typedef int OT;                   // Ordinal-type typdef

  ST one = Teuchos::ScalarTraits<ST>::one();
  ST zero = Teuchos::ScalarTraits<ST>::zero();
  
  //
  // Create 2D Laplacian using Tpetra::CisMatrix
  //

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
  const Tpetra::MpiPlatform<OT,OT>  ordinalPlatform(mpiComm);
  const Tpetra::MpiPlatform<OT,ST>   scalarPlatform(mpiComm);
#else
  const Tpetra::SerialPlatform<OT,OT>  ordinalPlatform;
  const Tpetra::SerialPlatform<OT,ST>   scalarPlatform;
#endif

  // Declare global dimension of the linear operator
  OT globalDim = 500;

  // Create the element space and std::vector space
  const Tpetra::ElementSpace<OT> elementSpace(globalDim,0,ordinalPlatform);
  const Tpetra::VectorSpace<OT,ST> vectorSpace(elementSpace,scalarPlatform);
  
  // Allocate the Tpetra::CisMatrix object.
  RCP<Tpetra::CisMatrix<OT,ST> >
    tpetra_A = rcp(new Tpetra::CisMatrix<OT,ST>(vectorSpace));

  // Allocate the Tpetra::CisMatrix preconditioning object.
  RCP<Tpetra::CisMatrix<OT,ST> >
    tpetra_Prec = rcp(new Tpetra::CisMatrix<OT,ST>(vectorSpace));

  // Get the indexes of the rows on this processor
  const int numMyElements = vectorSpace.getNumMyEntries();
  const std::vector<int> &myGlobalElements = vectorSpace.elementSpace().getMyGlobalElements();

  // Fill the local matrix entries one row at a time.
  const ST offDiag = -1.0, diag = 2.0;
  int numEntries; ST values[3]; int indexes[3];
  ST prec_value[1]; int prec_index[1];
  prec_value[0] = one/diag;   // Diagonal preconditioning
  for( int k = 0; k < numMyElements; ++k ) {
    const int rowIndex = myGlobalElements[k];
    prec_index[0] = rowIndex;
    if( rowIndex == 0 ) {                     // First row
      numEntries = 2;
      values[0]  = diag;             values[1]  = offDiag;
      indexes[0] = 0;                indexes[1] = 1; 
    }
    else if( rowIndex == globalDim - 1 ) {    // Last row
      numEntries = 2;
      values[0]  = offDiag;         values[1]  = diag;
      indexes[0] = globalDim-2;     indexes[1] = globalDim-1; 
    }
    else {                                    // Middle rows
      numEntries = 3;
      values[0]  = offDiag;         values[1]  = diag;          values[2]  = offDiag;
      indexes[0] = rowIndex-1;      indexes[1] = rowIndex;      indexes[2] = rowIndex+1; 
    }
    tpetra_A->submitEntries(Tpetra::Insert,rowIndex,numEntries,values,indexes);
    tpetra_Prec->submitEntries(Tpetra::Insert,rowIndex,1,prec_value,prec_index);
  }

  // Finish the construction of the Tpetra::CisMatrix
  tpetra_A->fillComplete();
  tpetra_Prec->fillComplete();

  // Create a Thyra linear operator (A) using the Tpetra::CisMatrix (tpetra_A).
  RCP<Thyra::LinearOpBase<ST> >
    A = Teuchos::rcp( new Thyra::TpetraLinearOp<OT,ST>(
                      Teuchos::rcp_implicit_cast<Tpetra::Operator<OT,ST> >(tpetra_A) ) );

  // Create a Thyra linear operator (Prec) using the Tpetra::CisMatrix (tpetra_Prec)
  RCP<Thyra::LinearOpBase<ST> >
    Prec = Teuchos::rcp( new Thyra::TpetraLinearOp<OT,ST>(
                    Teuchos::rcp_implicit_cast<Tpetra::Operator<OT,ST> >(tpetra_Prec) ) );

  // Create a Thyra default preconditioner (DefPrec) using the Thyra linear operator (Prec)
  RCP<const Thyra::DefaultPreconditioner<ST> >
    DefPrec = Thyra::unspecifiedPrec<ST>(Prec);

  //
  // Set the parameters for the Belos LOWS Factory and create a parameter list.
  // NOTE:  All the code below only uses Thyra and is independent of Tpetra.
  //
  int             blockSize              = 1;  // can only be 1 right now.
  int             maxIterations          = globalDim;
  int             maxRestarts            = 0;
  int             gmresKrylovLength      = globalDim;
  int             outputFrequency        = 100;
  bool            outputMaxResOnly       = true;
  MT              maxResid               = 1e-3;

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
 
  // Whether the linear solver succeeded.
  // (this will be set during the residual check at the end)
  bool success = true;

  // Number of random right-hand sides we will be solving for.
  int numRhs = 5;

  // Get the domain space for the Thyra linear operator 
  Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > domain = A->domain();

  // Create the Belos LOWS factory.
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ST> >
    belosLOWSFactory = Teuchos::rcp(new Thyra::BelosLinearOpWithSolveFactory<ST>());

  // Set the parameter list to specify the behavior of the factory.
  belosLOWSFactory->setParameterList( belosLOWSFPL );

  // Set the output stream and the verbosity level (prints to std::cout by defualt)
  // NOTE:  Set to VERB_NONE for no output from the solver.
  belosLOWSFactory->setVerbLevel(Teuchos::VERB_LOW);

  // Create a BelosLinearOpWithSolve object from the Belos LOWS factory.
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<ST> >
    nsA = belosLOWSFactory->createOp();

  // Initialize the BelosLinearOpWithSolve object with the Thyra linear operator (A)
  // and preconditioner (DefPrec).
  Thyra::initializeOp<ST>( *belosLOWSFactory, A, &*nsA );

  // Create a right-hand side with numRhs vectors in it and randomize it.
  Teuchos::RCP< Thyra::MultiVectorBase<ST> > 
    b = Thyra::createMembers(domain, numRhs);
  Thyra::seed_randomize<ST>(0);
  Thyra::randomize(-one, one, &*b);

  // Create an initial std::vector with numRhs vectors in it and initialize it to zero.
  Teuchos::RCP< Thyra::MultiVectorBase<ST> >
    x = Thyra::createMembers(domain, numRhs);
  Thyra::assign(&*x, zero);

  // Perform solve using the linear operator to get the approximate solution of Ax=b,
  // where b is the right-hand side and x is the left-hand side.
  Thyra::SolveStatus<ST> solveStatus;
  solveStatus = Thyra::solve( *nsA, Thyra::NONCONJ_ELE, *b, &*x );

  // Print out status of solve.
  *out << "\nBelos LOWS Status: "<< solveStatus << std::endl;

  //
  // Compute residual and ST check convergence.
  //
  std::vector<MT> norm_b(numRhs), norm_res(numRhs);
  Teuchos::RCP< Thyra::MultiVectorBase<ST> >
    y = Thyra::createMembers(domain, numRhs);

  // Compute the column norms of the right-hand side b.
  Thyra::norms_2( *b, &norm_b[0] );

  // Compute y=A*x, where x is the solution from the linear solver.
  A->apply( Thyra::NONCONJ_ELE, *x, &*y );
  
  // Compute A*x-b = y-b
  Thyra::update( -one, *b, &*y );

  // Compute the column norms of A*x-b.
  Thyra::norms_2( *y, &norm_res[0] );

  // Print out the final relative residual norms.
  MT rel_res = 0.0;
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
