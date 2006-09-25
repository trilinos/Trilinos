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
#include "NOX.H"  
#include "NOX_Epetra.H"
#include "NOX_TestCompare.H" // Test Suite headers
#include "NOX_Epetra_DebugTools.H"

// Trilinos headers
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"


int main(int argc, char *argv[]) 
{
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
 
  int * testInt = new int[100];
  delete [] testInt;

  bool verbose = false;

  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Set up theolver options parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> noxParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList & noxParams = *(noxParamsPtr.get());

  // Set up the printing utilities
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
    printParams.set("Output Information", NOX::Utils::Error +
		NOX::Utils::TestDetails);

  Teuchos::RefCountPtr<NOX::Utils> printing = Teuchos::rcp( new NOX::Utils(printParams) );

  // Identify the test problem
  if (printing->isPrintType(NOX::Utils::TestDetails))
    printing->out() << "Starting epetra/NOX_Operators/NOX_BroydenOp.exe" << endl;

  // Identify processor information
#ifdef HAVE_MPI
  if (printing->isPrintType(NOX::Utils::TestDetails)) 
  {
    printing->out() << "Parallel Run" << endl;
    printing->out() << "Number of processors = " << NumProc << endl;
    printing->out() << "Print Process = " << MyPID << endl;
  }
  Comm.Barrier();
  if (printing->isPrintType(NOX::Utils::TestDetails))
    printing->out() << "Process " << MyPID << " is alive!" << endl;
  Comm.Barrier();
#else
  if (printing->isPrintType(NOX::Utils::TestDetails))
    printing->out() << "Serial Run" << endl;
#endif

  int status = 0;

  // Create a TestCompare class
  NOX::Epetra::TestCompare tester( printing->out(), *printing);
  double abstol = 1.e-4;
  double reltol = 1.e-4 ;

  // Test NOX::Epetra::BroydenOperator
  int numGlobalElems = 3 * NumProc;
  Epetra_Map      broydenRowMap   ( numGlobalElems, 0, Comm );
  Epetra_Vector   broydenWorkVec  ( broydenRowMap );
  Epetra_CrsGraph broydenWorkGraph( Copy, broydenRowMap, 0 );
  vector<int> globalIndices(3);
  for( int lcol = 0; lcol < 3; ++lcol )
    globalIndices[lcol] = 3 * MyPID + lcol;

  vector<int> myGlobalIndices(2);

  // Row 1 structure
  myGlobalIndices[0] = globalIndices[0];
  myGlobalIndices[1] = globalIndices[2];
  broydenWorkGraph.InsertGlobalIndices( globalIndices[0], 2, &myGlobalIndices[0] );
  // Row 2 structure
  myGlobalIndices[0] = globalIndices[0];
  myGlobalIndices[1] = globalIndices[1];
  broydenWorkGraph.InsertGlobalIndices( globalIndices[1], 2, &myGlobalIndices[0] );
  // Row 3 structure
  myGlobalIndices[0] = globalIndices[1];
  myGlobalIndices[1] = globalIndices[2];
  broydenWorkGraph.InsertGlobalIndices( globalIndices[2], 2, &myGlobalIndices[0] );

  broydenWorkGraph.FillComplete();

  Teuchos::RefCountPtr<Epetra_CrsMatrix> broydenWorkMatrix = 
    Teuchos::rcp( new Epetra_CrsMatrix( Copy, broydenWorkGraph ) );

  // Create an identity matrix
  broydenWorkVec.PutScalar(1.0);
  broydenWorkMatrix->ReplaceDiagonalValues(broydenWorkVec);

  NOX::Epetra::BroydenOperator broydenOp( noxParams, printing, broydenWorkVec, broydenWorkMatrix, true );

  broydenWorkVec[0] =  1.0;
  broydenWorkVec[1] = -1.0;
  broydenWorkVec[2] =  2.0;
  broydenOp.setStepVector( broydenWorkVec );
  
  broydenWorkVec[0] =  2.0;
  broydenWorkVec[1] =  1.0;
  broydenWorkVec[2] =  3.0;
  broydenOp.setYieldVector( broydenWorkVec );

  broydenOp.computeSparseBroydenUpdate();

  // Create the gold matrix for comparison
  Teuchos::RefCountPtr<Epetra_CrsMatrix> goldMatrix = Teuchos::rcp( new Epetra_CrsMatrix( Copy, broydenWorkGraph ) );

  int      numCols ;
  double * values  ;

  // Row 1 answers
  goldMatrix->ExtractMyRowView( 0, numCols, values );
  values[0] =  6.0 ;
  values[1] =  2.0 ;
  // Row 2 answers
  goldMatrix->ExtractMyRowView( 1, numCols, values );
  values[0] =  5.0 ;
  values[1] =  0.0 ;
  // Row 3 structure
  goldMatrix->ExtractMyRowView( 2, numCols, values );
  values[0] = -1.0 ;
  values[1] =  7.0 ;

  goldMatrix->Scale(0.2);

  status += tester.testCrsMatrices( broydenOp.getBroydenMatrix(), *goldMatrix, reltol, abstol,
                              "Broyden Sparse Operator Update Test" );


  // Now try a dense Broyden Update
  Epetra_CrsGraph broydenWorkGraph2( Copy, broydenRowMap, 0 );

  myGlobalIndices.resize(3);

  // All Rowsstructure
  myGlobalIndices[0] = globalIndices[0];
  myGlobalIndices[1] = globalIndices[1];
  myGlobalIndices[2] = globalIndices[2];
  broydenWorkGraph2.InsertGlobalIndices( globalIndices[0], 3, &myGlobalIndices[0] );
  broydenWorkGraph2.InsertGlobalIndices( globalIndices[1], 3, &myGlobalIndices[0] );
  broydenWorkGraph2.InsertGlobalIndices( globalIndices[2], 3, &myGlobalIndices[0] );

  broydenWorkGraph2.FillComplete();

  Teuchos::RefCountPtr<Epetra_CrsMatrix> broydenWorkMatrix2 = Teuchos::rcp( new Epetra_CrsMatrix( Copy, broydenWorkGraph2 ) );

  // Create an identity matrix
  broydenWorkVec.PutScalar(1.0);
  broydenWorkMatrix2->ReplaceDiagonalValues(broydenWorkVec);

  NOX::Epetra::BroydenOperator broydenOp2( noxParams, printing, broydenWorkVec, broydenWorkMatrix2, true );

  broydenWorkVec[0] =  1.0;
  broydenWorkVec[1] = -1.0;
  broydenWorkVec[2] =  2.0;
  broydenOp2.setStepVector( broydenWorkVec );
  
  broydenWorkVec[0] =  2.0;
  broydenWorkVec[1] =  1.0;
  broydenWorkVec[2] =  3.0;
  broydenOp2.setYieldVector( broydenWorkVec );

  broydenOp2.computeSparseBroydenUpdate();

  // Create the gold matrix for comparison
  Teuchos::RefCountPtr<Epetra_CrsMatrix> goldMatrix2 = Teuchos::rcp( new Epetra_CrsMatrix( Copy, broydenWorkGraph2 ) );

  // Row 1 answers
  goldMatrix2->ExtractMyRowView( 0, numCols, values );
  values[0] =  7.0 ;
  values[1] = -1.0 ;
  values[2] =  2.0 ;
  // Row 2 answers
  goldMatrix2->ExtractMyRowView( 1, numCols, values );
  values[0] =  2.0 ;
  values[1] =  4.0 ;
  values[2] =  4.0 ;
  // Row 3 structure
  goldMatrix2->ExtractMyRowView( 2, numCols, values );
  values[0] =  1.0 ;
  values[1] = -1.0 ;
  values[2] =  8.0 ;

  double scaleF = 1.0 / 6.0;
  goldMatrix2->Scale( scaleF );

  status += tester.testCrsMatrices( broydenOp2.getBroydenMatrix(), *goldMatrix2, reltol, abstol,
                              "Broyden Sparse Operator Update Test (Dense)" );

  // Now test the ability to remove active entries in the Broyden update
  Epetra_CrsGraph inactiveGraph( Copy, broydenRowMap, 0 );

  // Row 1 structure
  inactiveGraph.InsertGlobalIndices( globalIndices[0], 1, &myGlobalIndices[1] );
  // Row 2 structure
  inactiveGraph.InsertGlobalIndices( globalIndices[1], 1, &myGlobalIndices[2] );
  // Row 3 structure
  inactiveGraph.InsertGlobalIndices( globalIndices[2], 1, &myGlobalIndices[0] );

  inactiveGraph.FillComplete();

  // Inactivate entries in dense matrix to arrive again at the original sparse structure
  broydenOp2.removeEntriesFromBroydenUpdate( inactiveGraph );

#ifdef HAVE_NOX_DEBUG
  if( verbose )
    broydenOp2.outputActiveEntries();
#endif

  // Reset to the identity matrix
  broydenOp2.resetBroydenMatrix( *broydenWorkMatrix2 );

  // Step and Yield vectors are already set
  broydenOp2.computeSparseBroydenUpdate();


  status += tester.testCrsMatrices( broydenOp2.getBroydenMatrix(), *goldMatrix, reltol, abstol,
                              "Broyden Sparse Operator Update Test (Entry Removal)", false );


  // Summarize test results  
  if( status == 0 )
    printing->out() << "Test passed!" << endl;
  else 
    printing->out() << "Test failed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return status;

}

/*
  end of file test.C
*/
