// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "sillyPowerMethod.hpp"
#include "createTridiagEpetraLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

//
// Increase the first diagonal element of your tridiagonal matrix.
//
void scaleFirstDiagElement( const double diagScale, Thyra::LinearOpBase<double> *A )
{
  using Teuchos::RCP;
  TEST_FOR_EXCEPT(A==NULL);
  // (A) Get at the underlying Epetra_Operator object that the EpetraLinearOp
  // object directly maintains.
  const RCP<Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator(*A);
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

  using Teuchos::outArg;
  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;
 
  bool success = true;
  bool result;
  
  //
  // (A) Setup and get basic MPI info
  //
  
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
#endif
  
  //
  // (B) Setup the output stream (do output only on root process!)
  //

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {

    //
    // (C) Read in commandline options
    //


    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int globalDim = 500;
    clp.setOption( "global-dim", &globalDim, "Global dimension of the linear system." );

    bool dumpAll = false;
    clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if (parse_return != CommandLineProcessor::PARSE_SUCCESSFUL)
      return parse_return;

    TEST_FOR_EXCEPTION( globalDim < 2, std::logic_error, "Error, globalDim=" << globalDim << " < 2 is not allowed!" );

    *out << "\n***\n*** Running power method example using Epetra implementation\n***\n" << std::scientific;

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
    *out << "\n(1) Constructing tridagonal Epetra matrix A of global dimension = " << globalDim << " ...\n";
    RCP<Epetra_Operator>
      A_epetra = createTridiagEpetraLinearOp(
        globalDim,
#ifdef HAVE_MPI
        mpiComm,
#endif
        1.0, true, *out
        );
    // Wrap in an Thyra::EpetraLinearOp object
    RCP<Thyra::LinearOpBase<double> >
      A = Thyra::nonconstEpetraLinearOp(A_epetra);
    //
    if (dumpAll) *out << "\nA =\n" << *A; // This works even in parallel!
    
    //
    // (2) Run the power method ANA
    //
    *out << "\n(2) Running the power method on matrix A ...\n";
    double  lambda      = 0.0;
    double  tolerance   = 1e-3;
    int     maxNumIters = 10*globalDim;
    result = sillyPowerMethod<double>(*A, maxNumIters, tolerance, outArg(lambda), *out);
    if(!result) success = false;
    *out << "\n  Estimate of dominate eigenvalue lambda = " << lambda << std::endl;
    
    //
    // (3) Increase dominance of first eigenvalue
    //
    *out << "\n(3) Scale the diagonal of A by a factor of 10 ...\n";
    scaleFirstDiagElement( 10.0, &*A );
    
    //
    // (4) Run the power method ANA again
    //
    *out << "\n(4) Running the power method again on matrix A ...\n";
    result = sillyPowerMethod<double>(*A, maxNumIters, tolerance, outArg(lambda), *out);
    if(!result) success = false;
    *out << "\n  Estimate of dominate eigenvalue lambda = " << lambda << std::endl;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)
  
  if (success)
    *out << "\nCongratulations! All of the tests checked out!\n";
  else
    *out << "\nOh no! At least one of the tests failed!\n";

  return success ? 0 : 1;

} // end main()
