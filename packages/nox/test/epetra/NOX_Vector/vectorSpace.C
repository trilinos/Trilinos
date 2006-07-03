//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
// NOX headers
#include "NOX.H"  // Required headers
#include "NOX_Epetra.H" // Epetra Interface headers
#include "NOX_TestError.H" // Test Suite headers
#include "NOX_TestCompare.H" // Test Suite headers
#include "Teuchos_RefCountPtr.hpp"

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
  Teuchos::RefCountPtr<Teuchos::ParameterList> noxParamsPtr =
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
  Teuchos::RefCountPtr<Epetra_Map> standardMap = 
    Teuchos::rcp(new Epetra_Map(globalLength, 0, Comm));
  
  // Create the L2 vector space
  Teuchos::RefCountPtr<NOX::Epetra::VectorSpace> vectorSpace1 =
    Teuchos::rcp(new NOX::Epetra::VectorSpaceL2);

  // Create the scaled L2 vector space
  Teuchos::RefCountPtr<Epetra_Vector> scaleVec =
    Teuchos::rcp(new Epetra_Vector(*standardMap,true));
  Teuchos::RefCountPtr<NOX::Epetra::Scaling> scaling = 
    Teuchos::rcp(new NOX::Epetra::Scaling);
  scaleVec->PutScalar(2.0);
  scaling->addUserScaling(NOX::Epetra::Scaling::Left, scaleVec);
  Teuchos::RefCountPtr<NOX::Epetra::VectorSpace> vectorSpace2 =
    Teuchos::rcp(new NOX::Epetra::VectorSpaceScaledL2(scaling));
  
  // Create a base Epetra_Vector 
  Teuchos::RefCountPtr<Epetra_Vector> epVec = 
    Teuchos::rcp(new Epetra_Vector(*standardMap, true));

  // Create NOX::Epetra::Vectors
  Teuchos::RefCountPtr<NOX::Epetra::Vector> vecL2 = 
    Teuchos::rcp(new NOX::Epetra::Vector(epVec, 
					 NOX::Epetra::Vector::CreateView, 
					 NOX::DeepCopy, 
					 vectorSpace1));
  Teuchos::RefCountPtr<NOX::Epetra::Vector> vecScaledL2 = 
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
  Teuchos::RefCountPtr<NOX::Abstract::Vector> tmpVector = 
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
    printing.out() << "Test passed!" << endl;
  else 
    printing.out() << "Test failed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // return 0 for a successful test
  return status;
}

/*
  end of file test.C
*/
