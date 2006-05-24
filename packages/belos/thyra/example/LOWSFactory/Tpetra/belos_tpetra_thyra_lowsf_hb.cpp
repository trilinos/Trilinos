
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
// Tpetra::CisMatrix.  This matrix is then converted into a Thyra linear operator
// through the Thyra-Tpetra adapters.
//
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero. 
//
//

// Thyra includes
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_TpetraLinearOp.hpp"

// Tpetra includes
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"

// I/O for Harwell-Boeing files
#ifdef HAVE_BELOS_TRIUTILS
#include "iohb.h"
#endif

// Teuchos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
#  include "Tpetra_MpiPlatform.hpp"
#else
#  include "Tpetra_SerialPlatform.hpp"
#endif

// My Tpetra::Operator include
#include "MyOperator.hpp"

int main(int argc, char* argv[])
{
  //
  // Get a default output stream from the Teuchos::VerboseObjectBase
  //
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;  // Scalar-type typedef
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;     // Scalar-type typedef
#else
  typedef double ST;                // Scalar-type typedef
#endif
  
  typedef Teuchos::ScalarTraits<ST>::magnitudeType MT;  // Magnitude-type typedef
  typedef int OT;                   // Ordinal-type typedef
  ST one = Teuchos::ScalarTraits<ST>::one(); 
  ST zero = Teuchos::ScalarTraits<ST>::zero(); 
  
#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
  const Tpetra::MpiPlatform<OT,OT>  ordinalPlatform(mpiComm);
  const Tpetra::MpiPlatform<OT,ST>   scalarPlatform(mpiComm);
#else
  const Tpetra::SerialPlatform<OT,OT>  ordinalPlatform;
  const Tpetra::SerialPlatform<OT,ST>   scalarPlatform;
#endif
  
  //
  // Get the data from the HB file
  //
  
  // Name of input matrix file
  std::string matrixFile = "mhd1280b.cua";
  
  int info=0;
  int dim,dim2,nnz;
  MT *dvals;
  int *colptr,*rowind;
  ST *cvals;
  nnz = -1;
  info = readHB_newmat_double(matrixFile.c_str(),&dim,&dim2,&nnz,
                              &colptr,&rowind,&dvals);

  if (info == 0 || nnz < 0) {
    *out << "Error reading '" << matrixFile << "'" << endl;
  }
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Convert interleaved doubles to complex values
  cvals = new ST[nnz];
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
  }
  
  // Declare global dimension of the linear operator
  OT globalDim = dim;
  
  // Create the element space and vector space
  const Tpetra::ElementSpace<OT> elementSpace(globalDim,0,ordinalPlatform);
  const Tpetra::VectorSpace<OT,ST> vectorSpace(elementSpace,scalarPlatform);
  
  // Create my implementation of a Tpetra::Operator
  RefCountPtr<Tpetra::Operator<OT,ST> >
    tpetra_A = rcp( new MyOperator<OT,ST>(vectorSpace,dim,colptr,nnz,rowind,cvals) );

  // Create a Thyra linear operator (A) using the Tpetra::CisMatrix (tpetra_A).
  RefCountPtr<Thyra::LinearOpBase<ST> >
    A = Teuchos::rcp( new Thyra::TpetraLinearOp<OT,ST>(tpetra_A) );

  //
  // Set the parameters for the Belos LOWS Factory and create a parameter list.
  //
  int             blockSize              = 1;
  int             maxIterations          = globalDim;
  int             maxRestarts            = 15;
  int             gmresKrylovLength      = 50;
  int             outputFrequency        = 100;
  bool            outputMaxResOnly       = true;
  MT              maxResid               = 1e-5;

  Teuchos::RefCountPtr<Teuchos::ParameterList>
    belosLOWSFPL = Teuchos::rcp( new Teuchos::ParameterList() );
  
  belosLOWSFPL->set("Solver Type","GMRES");
  belosLOWSFPL->set("Max Iters",int(maxIterations));
  belosLOWSFPL->set("Default Rel Res Norm",MT(maxResid));
  belosLOWSFPL->set("Max Restarts",int(maxRestarts));
  belosLOWSFPL->set("Block Size",int(blockSize));
  belosLOWSFPL->sublist("GMRES").set("Max Number of Krylov Vectors",int(gmresKrylovLength*blockSize));
  belosLOWSFPL->sublist("GMRES").set("Variant","Standard");
  Teuchos::ParameterList &outputterSL = belosLOWSFPL->sublist("Outputter");
  outputterSL.set("Output Frequency",int(outputFrequency));
  outputterSL.set("Output Max Res Only",bool(outputMaxResOnly));

  // Whether the linear solver succeeded.
  // (this will be set during the residual check at the end)
  bool success = true;

  // Number of random right-hand sides we will be solving for.
  int numRhs = 1;

  // Get the domain space for the Thyra linear operator 
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<ST> > domain = A->domain();

  // Create the Belos LOWS factory.
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<ST> >
    belosLOWSFactory = Teuchos::rcp(new Thyra::BelosLinearOpWithSolveFactory<ST>());

  // Set the parameter list to specify the behavior of the factory.
  belosLOWSFactory->setParameterList( belosLOWSFPL );

  // Set the output stream and the verbosity level (prints to std::cout by defualt)
  // NOTE:  Set to VERB_NONE for no output from the solver.
  belosLOWSFactory->setVerbLevel(Teuchos::VERB_LOW);

  // Create a BelosLinearOpWithSolve object from the Belos LOWS factory.
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<ST> >
    nsA = belosLOWSFactory->createOp();

  // Initialize the BelosLinearOpWithSolve object with the Thyra linear operator.
  belosLOWSFactory->initializeOp( A, &*nsA );

  // Create a right-hand side with numRhs vectors in it.
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<ST> > 
    b = Thyra::createMembers(domain, numRhs);

  // Create an initial vector with numRhs vectors in it and initialize it to one.
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<ST> >
    x = Thyra::createMembers(domain, numRhs);
  Thyra::assign(&*x, one);

  // Initialize the right-hand side so that the solution is a vector of ones.
  A->apply( Thyra::NONCONJ_ELE, *x, &*b );
  Thyra::assign(&*x, zero);

  // Perform solve using the linear operator to get the approximate solution of Ax=b,
  // where b is the right-hand side and x is the left-hand side.
  Thyra::SolveStatus<ST> solveStatus;
  solveStatus = Thyra::solve( *nsA, Thyra::NONCONJ_ELE, *b, &*x );

  // Print out status of solve.
  *out << "\nBelos LOWS Status: "<< solveStatus << endl;

  //
  // Compute residual and ST check convergence.
  //
  std::vector<MT> norm_b(numRhs), norm_res(numRhs);
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<ST> >
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
