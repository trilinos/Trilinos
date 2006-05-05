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

#include "sillyPowerMethod.hpp"
#include "createTridiagEpetraLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

//
// Increase the first diagonal element of your tridiagonal matrix.
//
void scaleFirstDiagElement( const double diagScale, Thyra::LinearOpBase<double> *A )
{
  using Teuchos::RefCountPtr;
  TEST_FOR_EXCEPT(A==NULL);
  // (A) Get at the underlying Epetra_Operator object that the EpetraLinearOp
  // object directly maintains.
  const RefCountPtr<Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator(*A);
  // (B) Perform a dynamic cast to Epetra_CrsMatrix.
  // Note, the dyn_cast<>() template function will throw std::bad_cast
  // with a nice error message if the cast fails! 
  Epetra_CrsMatrix &crsMatrix = Teuchos::dyn_cast<Epetra_CrsMatrix>(*epetra_op);
  // (C) Query if the first row of the matrix is on this processor and if
  // it is get it and scale it.
  if(crsMatrix.MyGlobalRow(0)) {
    const int numRowNz = crsMatrix.NumGlobalEntries(0);
    TEST_FOR_EXCEPT( numRowNz != 2 );
    int returndNumRowNz; double rowValues[2]; int rowIndexes[2];
    crsMatrix.ExtractGlobalRowCopy( 0, numRowNz, returndNumRowNz, &rowValues[0], &rowIndexes[0] );
    TEST_FOR_EXCEPT( returndNumRowNz != 2 );
    for( int k = 0; k < numRowNz; ++k) if (rowIndexes[k] == 0) rowValues[k] *= diagScale;
    crsMatrix.ReplaceGlobalValues( 0, numRowNz, rowValues, rowIndexes );
  }
} // end scaleFirstDiagElement()

//
// Main driver program for epetra implementation of the power method.
//
int main(int argc, char *argv[])
{
  using Teuchos::CommandLineProcessor;
  using Teuchos::RefCountPtr;
 
  bool success = true;
  bool verbose = true;
  bool result;
  int procRank = 0;
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif
  
  //
  // (A) Get basic MPI info
  //
  
#ifdef HAVE_MPI
  int numProc;
  MPI_Comm mpiComm = MPI_COMM_WORLD;
  MPI_Comm_size( mpiComm, &numProc );
  MPI_Comm_rank( mpiComm, &procRank );
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
    bool   dumpAll                    = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "global-dim", &globalDim, "Global dimension of the linear system." );
    clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPTION( globalDim < 2, std::logic_error, "Error, globalDim=" << globalDim << " < 2 is not allowed!" );

    if(verbose) *out << "\n***\n*** Running power method example using Epetra implementation\n***\n" << std::scientific;

    //
    // (D) Setup the operator and run the power method!
    //

    //
    // (1) Setup the initial tridiagonal operator
    //
    //       [  2  -1             ]
    //       [ -1   2  -1         ]
    //  A =  [      .   .   .     ]
    //       [          -1  2  -1 ]
    //       [             -1   2 ]
    //
    if(verbose) *out << "\n(1) Constructing tridagonal Epetra matrix A of global dimension = " << globalDim << " ...\n";
    RefCountPtr<Thyra::LinearOpBase<double> >
      A = createTridiagEpetraLinearOp(
        globalDim
#ifdef HAVE_MPI
        ,mpiComm
#endif
        ,1.0,verbose,*out
        );
    if( verbose && dumpAll ) *out << "\nA =\n" << *A; // This works even in parallel!
    
    //
    // (2) Run the power method ANA
    //
    if(verbose) *out << "\n(2) Running the power method on matrix A ...\n";
    double  lambda      = 0.0;
    double  tolerance   = 1e-3;
    int     maxNumIters = 10*globalDim;
    result = sillyPowerMethod<double>(*A,maxNumIters,tolerance,&lambda,(verbose?&*out:NULL));
    if(!result) success = false;
    if(verbose) *out << "\n  Estimate of dominate eigenvalue lambda = " << lambda << std::endl;
    
    //
    // (3) Increase dominance of first eigenvalue
    //
    if(verbose) *out << "\n(3) Scale the diagonal of A by a factor of 10 ...\n";
    scaleFirstDiagElement( 10.0, &*A );
    
    //
    // (4) Run the power method ANA again
    //
    if(verbose) *out << "\n(4) Running the power method again on matrix A ...\n";
    result = sillyPowerMethod<double>(*A,maxNumIters,tolerance,&lambda,(verbose?&*out:NULL));
    if(!result) success = false;
    if(verbose) *out << "\n  Estimate of dominate eigenvalue lambda = " << lambda << std::endl;
    
  }
  catch( const std::exception &excpt ) {
    std::cerr << "p="<<procRank<<": *** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    std::cerr << "p="<<procRank<<": *** Caught an unknown exception\n";
    success = false;
  }
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

#ifdef HAVE_MPI
   MPI_Finalize();
#endif

  return success ? 0 : 1;

} // end main()
