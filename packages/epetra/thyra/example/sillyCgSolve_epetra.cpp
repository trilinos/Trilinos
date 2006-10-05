// @HEADER
// ***********************************************************************
// 
//                Epetra: Linear Algebra Services Package 
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

#include "createTridiagEpetraLinearOp.hpp"
#include "sillyCgSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

//
// Main driver program for epetra implementation of CG.
//
int main(int argc, char *argv[])
{
  using Teuchos::CommandLineProcessor;
  using Teuchos::RefCountPtr;
 
  bool success = true;
  bool verbose = true;
  bool result;
  
  //
  // (A) Setup and get basic MPI info
  //

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int procRank = Teuchos::GlobalMPISession::getRank();
  const int numProc  = Teuchos::GlobalMPISession::getNProc();
#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
#endif

  //
  // (B) Setup the output stream (do output only on root process!)
  //
  
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {

    //
    // (C) Read in commandline options
    //

    int    globalDim                  = 500;
    double diagScale                  = 1.001;
    bool   useWithNonEpetraVectors    = false;
    double tolerance                  = 1e-4;
    int    maxNumIters                = 300;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "global-dim", &globalDim, "Global dimension of the linear system." );
    clp.setOption( "diag-scale", &diagScale, "Scaling of the diagonal to improve conditioning." );
    clp.setOption( "use-with-non-epetra-vectors", "use-with-epetra-vectors", &useWithNonEpetraVectors
                   , "Use non-epetra vectors with Epetra operator or not." );
    clp.setOption( "tol", &tolerance, "Relative tolerance for linear system solve." );
    clp.setOption( "max-num-iters", &maxNumIters, "Maximum of CG iterations." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPTION( globalDim < 2, std::logic_error, "Error, globalDim=" << globalDim << " < 2 is not allowed!" );

    if(verbose) *out << "\n***\n*** Running CG example using Epetra implementation\n***\n" << std::scientific;

    //
    // (D) Setup a simple linear system with tridagonal operator:
    //
    //       [  a*2   -1                ]
    //       [ -1    a*2  -1            ]
    //  A =  [         .   .    .       ]
    //       [            -1  a*2    -1 ]
    //       [                 -1   a*2 ]
    //
    // (D.1) Create the tridagonal Epetra matrix operator
    if(verbose) *out << "\n(1) Constructing tridagonal Epetra matrix A_epetra of global dimension = " << globalDim << " ...\n";
    RefCountPtr<Epetra_Operator>
      A_epetra = createTridiagEpetraLinearOp(
        globalDim
#ifdef HAVE_MPI
        ,mpiComm
#endif
        ,diagScale,verbose,*out
        );
    // (D.2) Wrap the Epetra_Opertor in a Thyra::EpetraLinearOp object
    RefCountPtr<Thyra::LinearOpBase<double> >
      A = rcp(new Thyra::EpetraLinearOp(A_epetra));
    // (D.3) Create RHS vector b and set to a random value
    RefCountPtr<const Thyra::VectorSpaceBase<double> >
      b_space = A->range();
    // (D.4) Create the RHS vector b and initialize it to a random vector
    RefCountPtr<Thyra::VectorBase<double> > b = createMember(b_space);
    Thyra::seed_randomize<double>(0);
    Thyra::randomize( -1.0, +1.0, &*b );
    // (D.5) Create LHS vector x and set to zero
    RefCountPtr<Thyra::VectorBase<double> > x = createMember(A->domain());
    Thyra::assign( &*x, 0.0 );
    //
    // (E) Solve the linear system with the silly CG solver
    //
    result = sillyCgSolve(*A,*b,maxNumIters,tolerance,&*x,verbose?&*out:0);
    if(!result) success = false;
    //
    // (F) Check that the linear system was solved to the specified tolerance
    //
    RefCountPtr<Thyra::VectorBase<double> > r = createMember(A->range());                     
    Thyra::assign(&*r,*b);                                        // r = b
    Thyra::apply(*A,Thyra::NOTRANS,*x,&*r,-1.0,+1.0);             // r = -A*x + r
    const double r_nrm = Thyra::norm(*r), b_nrm =Thyra:: norm(*b);
    const double rel_err = r_nrm/b_nrm, relaxTol = 2.0*tolerance;
    result = rel_err <= relaxTol;
    if(!result) success = false;
    if(verbose)
      *out
        << "\n||b-A*x||/||b|| = "<<r_nrm<<"/"<<b_nrm<<" = "<<rel_err<<(result?" <= ":" > ")
        <<"2.0*tolerance = "<<relaxTol<<": "<<(result?"passed":"failed")<<std::endl;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return success ? 0 : 1;

} // end main()
