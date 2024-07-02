// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   quick_solve.cpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Thu Jul 14 16:24:46 MDT 2011

   \brief  Intended to quickly check a solver interface

   This example solves a simple sparse system of linear equations
   using a given Amesos2 solver interface.
*/

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Map.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <EpetraExt_CrsMatrixIn.h>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef HAVE_MPI
  const Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  const Epetra_SerialComm comm;
#endif
  size_t myRank = comm.MyPID();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  *fos << myRank << " : " << Amesos2::version() << std::endl << std::endl;

  bool printMatrix   = false;
  bool printSolution = false;
  bool printTiming   = false;
  bool verbose       = false;
  std::string solver_name = "SuperLU";
  std::string filedir = "../test/matrices/";
  std::string filename = "arc130.mtx";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filedir",&filedir,"Directory where matrix-market files are located");
  cmdp.setOption("filename",&filename,"Filename for Matrix-Market test matrix.");
  cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("print-solution","no-print-solution",&printSolution,"Print solution vector after solve.");
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");
  cmdp.setOption("solver", &solver_name, "Which TPL solver library to use.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Before we do anything, check that the solver is enabled
  if( !Amesos2::query(solver_name) ){
    std::cerr << solver_name << " not enabled.  Exiting..." << std::endl;
    return EXIT_SUCCESS;	// Otherwise CTest will pick it up as
				// failure, which it isn't really
  }

  const size_t numVectors = 1;

  std::string mat_pathname = filedir + filename;
  
  MAT* A;
  int ret = EpetraExt::MatrixMarketFileToCrsMatrix(mat_pathname.c_str(), comm, A, false, false);
  if( ret == -1 ){
    *fos << "error reading matrix file from disk, aborting..." << std::endl;
    return EXIT_FAILURE;
  }

  if( printMatrix ){
    A->Print(*(fos->getOStream()));
  }

  // get the maps
  const Epetra_Map dmnmap = A->DomainMap();
  const Epetra_Map rngmap = A->RangeMap();

  // Create random X
  RCP<MV> X = rcp( new MV(dmnmap,numVectors) );
  X->Random();

  RCP<MV> B = rcp(new MV(rngmap,numVectors));
  B->Random();

  // Constructor from Factory
  RCP<Amesos2::Solver<MAT,MV> > solver;
  try{
    solver = Amesos2::create<MAT,MV>(solver_name, rcp(A), X, B);
  } catch (const std::invalid_argument& e){
    *fos << e.what() << std::endl;
    return 0;
  }

  #ifdef HAVE_AMESOS2_SHYLU_NODEBASKER
  if( Amesos2::query("shylubasker") ) {
    Teuchos::ParameterList amesos2_params("Amesos2");
      amesos2_params.sublist(solver_name).set("num_threads", 1, "Number of threads");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  }
  #endif

  solver->solve();

  if( printSolution ){
    // Print the solution
    X->Print(*(fos->getOStream()));
  }
  
  if( printTiming ){
    // Print some timing statistics
    solver->printTiming(*fos);
  }

  // We are done.
  return 0;
}
