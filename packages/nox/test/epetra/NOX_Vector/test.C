//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

  // Identify the test problem
  if (printing.isPrintType(NOX::Utils::TestDetails))
    printing.out() << "Starting epetra/NOX_Vector/NOX_Vector.exe" << endl;

  // Create a TestCompare class
  NOX::TestCompare tester( printing.out(), printing);
  double tolerance = 1.e-12;
  NOX::TestCompare::CompareType aComp = NOX::TestCompare::Absolute;

  // Identify processor information
#ifdef HAVE_MPI
  printing.out() << "Parallel Run" << endl;
  printing.out() << "Number of processors = " << Comm.NumProc() << endl;
  printing.out() << "Print Process = " << MyPID << endl;
  Comm.Barrier();
  if (printing.isPrintType(NOX::Utils::TestDetails))
    printing.out() << "Process " << MyPID << " is alive!" << endl;
  Comm.Barrier();
#else
  printing.out() << "Serial Run" << endl;
#endif

  // Create a map describing data distribution
  Epetra_Map * standardMap = new Epetra_Map(globalLength, 0, Comm);
  
  // Return value
  int status = 0;

  // *** Start Testing Here!!! ***

  // First create the Epetra_Vector needed to construct our NOX vector
  Epetra_Vector * epetraVec = new Epetra_Vector(*standardMap, true);

  NOX::Epetra::Vector * noxVec1 = new NOX::Epetra::Vector(*epetraVec, NOX::DeepCopy);
  delete epetraVec; epetraVec = 0;

  NOX::Epetra::Vector * noxVec2 = new NOX::Epetra::Vector(*noxVec1);
  noxVec2->init(1.0);

  // Test our norms
  NOX::Abstract::Vector::NormType 
    oneNorm = NOX::Abstract::Vector::OneNorm,
    twoNorm = NOX::Abstract::Vector::TwoNorm,
    infNorm = NOX::Abstract::Vector::MaxNorm;

  double expectedOneNorm = (double) globalLength,
         expectedTwoNorm = sqrt( (double) globalLength),
         expectedInfNorm = 1.0;

  status += tester.testValue( noxVec2->norm(oneNorm), expectedOneNorm,
                              tolerance, "One-Norm Test", aComp);
  status += tester.testValue( noxVec2->norm(twoNorm), expectedTwoNorm,
                              tolerance, "Two-Norm Test", aComp);
  status += tester.testValue( noxVec2->norm(infNorm), expectedInfNorm,
                              tolerance, "Max-Norm Test", aComp);


  // Test random, reciprocal and dot methods
  noxVec1->random();
  // Threshold values since we want to do a reciprocal
  int myLength = standardMap->NumMyElements();
  for( int i = 0; i < myLength; ++i )
    if( fabs(noxVec1->getEpetraVector()[i]) < 1.e-8 ) noxVec1->getEpetraVector()[i] = 1.e-8;

  noxVec2->reciprocal(*noxVec1);
  double product = noxVec1->innerProduct(*noxVec2);

  status += tester.testValue( product, expectedOneNorm,
                              tolerance, "Random, Reciprocal and Dot Test", aComp);


  // Test abs and weighted-norm methods
  /*  ----------------------------
      NOT SUPPORTED AT THIS TIME
      ----------------------------
  noxVec2->abs(*noxVec2);
  double wNorm = noxVec1->norm(*noxVec2);
  status += tester.testValue( wNorm, noxVec1->norm(oneNorm),
                              tolerance, "Abs and Weighted-Norm Test", aComp);
  */

  // Test operator= , abs, update and scale methods
  (*noxVec2) = (*noxVec1);
  noxVec2->abs(*noxVec2);
  double sumAll = noxVec1->norm(oneNorm);
  noxVec2->update( 1.0, *noxVec1, 1.0 );
  noxVec2->scale(0.5);
  double sumPositive = noxVec2->norm(oneNorm);
  (*noxVec2) = (*noxVec1);
  noxVec2->abs(*noxVec2);
  noxVec2->update( 1.0, *noxVec1, -1.0 );
  noxVec2->scale(0.5);
  double sumNegative = noxVec2->norm(oneNorm);

  status += tester.testValue( (sumPositive + sumNegative), sumAll,
                              tolerance, "Abs, Operator= , Update and Scale Test", aComp);



  if (status == 0)
    printing.out() << "Test passed!" << endl;
  else 
    printing.out() << "Test failed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  delete noxVec2;
  delete noxVec1;
  delete standardMap;

  // return 0 for a successful test
  return status;
}

/*
  end of file test.C
*/
