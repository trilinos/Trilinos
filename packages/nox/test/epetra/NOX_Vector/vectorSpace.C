//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
// NOX headers
#include "NOX.H"  // Required headers
#include "NOX_Epetra.H" // Epetra Interface headers
#include "NOX_TestError.H" // Test Suite headers
#include "NOX_TestCompare.H" // Test Suite headers
#include "Teuchos_RCP.hpp"

// Trilinos headers
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

int main(int argc, char *argv[]) {

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
 
  int globalLength = 100; // This should suffice

  bool verbose = false;

  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();

  // Set up the printing utilities
  Teuchos::RCP<Teuchos::ParameterList> noxParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& noxParams = *(noxParamsPtr.get());
  // Only print output if the "-v" flag is set on the command line
  Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 5);
  printParams.set("Output Processor", 0);
  if( verbose )
    printParams.set("Output Information", 
		NOX::Utils::OuterIteration + 
		NOX::Utils::OuterIterationStatusTest + 
		NOX::Utils::InnerIteration +
		NOX::Utils::Parameters + 
		NOX::Utils::Details + 
		NOX::Utils::Warning +
		NOX::Utils::TestDetails);
  else
    printParams.set("Output Information", NOX::Utils::Error);
  
  NOX::Utils printing(printParams);

  // Return value
  int status = 0;

  // *** Start Testing Here!!! ***

  // Create a map describing data distribution
  Teuchos::RCP<Epetra_Map> standardMap = 
    Teuchos::rcp(new Epetra_Map(globalLength, 0, Comm));
  
  // Create the L2 vector space
  Teuchos::RCP<NOX::Epetra::VectorSpace> vectorSpace1 =
    Teuchos::rcp(new NOX::Epetra::VectorSpaceL2);

  // Create the scaled L2 vector space
  Teuchos::RCP<Epetra_Vector> scaleVec =
    Teuchos::rcp(new Epetra_Vector(*standardMap,true));
  Teuchos::RCP<NOX::Epetra::Scaling> scaling = 
    Teuchos::rcp(new NOX::Epetra::Scaling);
  scaleVec->PutScalar(2.0);
  scaling->addUserScaling(NOX::Epetra::Scaling::Left, scaleVec);
  Teuchos::RCP<NOX::Epetra::VectorSpace> vectorSpace2 =
    Teuchos::rcp(new NOX::Epetra::VectorSpaceScaledL2(scaling));
  
  // Create a base Epetra_Vector 
  Teuchos::RCP<Epetra_Vector> epVec = 
    Teuchos::rcp(new Epetra_Vector(*standardMap, true));

  // Create NOX::Epetra::Vectors
  Teuchos::RCP<NOX::Epetra::Vector> vecL2 = 
    Teuchos::rcp(new NOX::Epetra::Vector(epVec, 
					 NOX::Epetra::Vector::CreateView, 
					 NOX::DeepCopy, 
					 vectorSpace1));
  Teuchos::RCP<NOX::Epetra::Vector> vecScaledL2 = 
    Teuchos::rcp(new NOX::Epetra::Vector(*epVec, NOX::DeepCopy, vectorSpace2));
  
  // Test our norms
  NOX::Abstract::Vector::NormType oneNorm = NOX::Abstract::Vector::OneNorm;
  NOX::Abstract::Vector::NormType twoNorm = NOX::Abstract::Vector::TwoNorm;
  NOX::Abstract::Vector::NormType infNorm = NOX::Abstract::Vector::MaxNorm;

  double expectedInfNorm = 1.0;
  double expectedOneNorm = static_cast<double>(globalLength);
  double expectedTwoNorm = sqrt(static_cast<double>(globalLength));

  double expectedScaledInfNorm = 0.5;
  double expectedScaledOneNorm = static_cast<double>(globalLength) / 2.0;
  double expectedScaledTwoNorm = 0.0;
  for (int i = 0; i < globalLength; i++)
    expectedScaledTwoNorm += 0.5 * 0.5;
  expectedScaledTwoNorm = sqrt(expectedScaledTwoNorm);

  // Create a TestCompare class
  NOX::TestCompare tester(printing.out(), printing);
  double tolerance = 1.e-12;
  NOX::TestCompare::CompareType aComp = NOX::TestCompare::Absolute;

  vecL2->init(1.0);
  vecScaledL2->init(1.0);

  status += tester.testValue(vecL2->norm(infNorm), expectedInfNorm,
			     tolerance, "Max-Norm Test", aComp);
  status += tester.testValue(vecL2->norm(oneNorm), expectedOneNorm,
			     tolerance, "One-Norm Test", aComp);
  status += tester.testValue(vecL2->norm(twoNorm), expectedTwoNorm,
			     tolerance, "Two-Norm Test", aComp);
  
  status += tester.testValue(vecScaledL2->norm(infNorm), expectedScaledInfNorm,
			     tolerance, "Scaled Max-Norm Test", aComp);
  status += tester.testValue(vecScaledL2->norm(oneNorm), expectedScaledOneNorm,
			     tolerance, "Scaled One-Norm Test", aComp);
  status += tester.testValue(vecScaledL2->norm(twoNorm), expectedScaledTwoNorm,
			     tolerance, "Scaled Two-Norm Test", aComp);

  // Test the inner products for L2
  Teuchos::RCP<NOX::Abstract::Vector> tmpVector = 
    vecL2->clone();
  vecL2->init(2.0);
  tmpVector->init(1.0);
  double expectedInnerProduct = static_cast<double>(globalLength) * 2.0;

  status += tester.testValue(vecL2->innerProduct(*tmpVector), 
			     expectedInnerProduct,
			     tolerance, 
			     "Inner Product Test", 
			     aComp);

  // Test the inner products for Scaled L2
  tmpVector = vecScaledL2->clone();
  vecScaledL2->init(2.0);
  tmpVector->init(1.0);
  expectedInnerProduct = static_cast<double>(globalLength) * 0.5;

  status += tester.testValue(vecScaledL2->innerProduct(*tmpVector), 
			     expectedInnerProduct,
			     tolerance, 
			     "Scaled Inner Product Test", 
			     aComp);

  if (status == 0)
    printing.out() << "Test passed!" << std::endl;
  else 
    printing.out() << "Test failed!" << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // return 0 for a successful test
  return status;
}

/*
  end of file test.C
*/
