// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
//  This test is for the BlockDavidson solver
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBlockDavidson.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziStatusTestMaxIters.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};

using namespace Anasazi;

  
typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;
typedef double ScalarType;
typedef MultiVecTraits<ScalarType,MV>    MVT;
typedef OperatorTraits<ScalarType,MV,OP> OPT;


int main(int argc, char *argv[]) 
{
  
#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  Epetra_SerialComm Comm;

#endif

  int MyPID = Comm.MyPID();

  bool testFailed = false;
  bool verbose = 0;
  std::string which("LM");
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
    else {
      which = argv[1];
    }
  }
  if (argc>2) {
    if (argv[2][0]=='-' && argv[2][1]=='v') {
      verbose = true;
    }
    else {
      which = argv[2];
    }
  }

  if (verbose && MyPID == 0) {
    cout << Anasazi_Version() << endl << endl;
  }

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  int nev = 4;
  //
  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create the initial vectors
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), nev) );
  ivec->Random();
  //
  // Create eigenproblem
  Teuchos::RefCountPtr<BasicEigenproblem<ScalarType, MV, OP> > problem =
    Teuchos::rcp( new BasicEigenproblem<ScalarType, MV, OP>(K, M, ivec) );
  //  problem->SetPrec( Teuchos::rcp( const_cast<Epetra_Operator *>(opStiffness->getPreconditioner()), false ) );
  //
  // Inform the eigenproblem that the operator A is symmetric
  problem->setHermitian(true);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are finishing passing it information
  if (problem->setProblem() != true) {  
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned false." << endl;
  }

  // create the output manager
  Teuchos::RefCountPtr< OutputManager<ScalarType> > printer = Teuchos::rcp( new BasicOutputManager<ScalarType>() );
  printer->setVerbosity( Anasazi::Debug );
  // create the orthogonalization manager
  Teuchos::RefCountPtr< MatOrthoManager<ScalarType,MV,OP> > ortho = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(M) );
  // create the sort manager
  Teuchos::RefCountPtr< SortManager<ScalarType,MV,OP> > sorter = Teuchos::rcp( new BasicSort<ScalarType,MV,OP>(which) );
  // create the status tester
  Teuchos::RefCountPtr< StatusTest<ScalarType,MV,OP> > tester = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>(1) );
  // create the parameter list
  Teuchos::ParameterList pls;
  pls.set<int>("Block Size",5);
  pls.set<int>("Num Blocks",6);
  // create the solver
  Teuchos::RefCountPtr< BlockDavidson<ScalarType,MV,OP> > solver = Teuchos::rcp( new BlockDavidson<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,pls) );

  try 
  {
    BlockDavidsonState<ScalarType,MV> state;
    // solver should be uninitialized, with block size and full ortho as specified by the ParameterList
    state = solver->getState();
    TEST_FOR_EXCEPTION(solver->isInitialized() != false,get_out,"Solver should be un-initialized after instantiation.");  
    TEST_FOR_EXCEPTION(solver->getBlockSize() != pls.get<int>("Block Size"),get_out,"Solver block size does not match ParameterList.");  
    TEST_FOR_EXCEPTION(solver->getSize()/solver->getBlockSize() != pls.get<int>("Num Blocks"),get_out,"Solver num blocks does not match ParameterList.");  
    TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero.")
    TEST_FOR_EXCEPTION(&solver->getProblem() != problem.get(),get_out,"getProblem() did not return the submitted problem.");
    TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.V) != solver->getSize(),get_out,"blockSize() does not match allocated size for V");
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.X) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for X");
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.R) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for R");

    // initialize solver and perform checks: initial iterate from eigenproblem, residual computed correctly, etc
    solver->initialize();
    state = solver->getState();
    TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");  
    TEST_FOR_EXCEPTION(solver->getBlockSize() != pls.get<int>("Block Size"),get_out,"Solver block size does not match ParameterList.");  
    TEST_FOR_EXCEPTION(solver->getSize()/solver->getBlockSize() != pls.get<int>("Num Blocks"),get_out,"Solver num blocks does not match ParameterList.");  
    TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero.")
    TEST_FOR_EXCEPTION(&solver->getProblem() != problem.get(),get_out,"getProblem() did not return the submitted problem.");
    TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.V) != solver->getSize(),get_out,"blockSize() does not match allocated size for V");
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.X) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for X");
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.R) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for R");

  }
  catch (get_out go) {
    printer->stream(Errors) << "Test failed: " << go.what() << endl;
    testFailed = true;
  }
  catch (std::exception e) {
    printer->stream(Errors) << "Caught unexpected exception: " << e.what() << endl;
    testFailed = true;
  }

  
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}	
