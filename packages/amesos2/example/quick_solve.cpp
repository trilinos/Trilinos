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
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <MatrixMarket_Tpetra.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  typedef Tpetra::CrsMatrix<>::scalar_type Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  size_t myRank = comm->getRank();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  *fos << Amesos2::version() << std::endl << std::endl;

  bool printMatrix   = false;
  bool printSolution = false;
  bool printTiming   = false;
  bool printResidual = false;
  bool printLUStats  = false;
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
  cmdp.setOption("print-residual","no-print-residual",&printResidual,"Print solution residual");
  cmdp.setOption("print-lu-stats","no-print-lu-stats",&printLUStats,"Print nnz in L and U factors");
  cmdp.setOption("solver", &solver_name, "Which TPL solver library to use.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Before we do anything, check that the solver is enabled
  if( !Amesos2::query(solver_name) ){
    std::cerr << solver_name << " not enabled.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  const size_t numVectors = 1;

  // create a Map
  global_size_t nrows = 6;
  RCP<Tpetra::Map<LO,GO> > map
    = rcp( new Tpetra::Map<LO,GO>(nrows,0,comm) );

  std::string mat_pathname = filedir + filename;
  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(mat_pathname,comm);

  if( printMatrix ){
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if( verbose && myRank==0 ){
    *fos << std::endl << A->description() << std::endl << std::endl;
  }

  // get the maps
  RCP<const Tpetra::Map<LO,GO> > dmnmap = A->getDomainMap();
  RCP<const Tpetra::Map<LO,GO> > rngmap = A->getRangeMap();

  // Create random X
  RCP<MV> Xhat = rcp( new MV(dmnmap,numVectors) );
  RCP<MV> X = rcp( new MV(dmnmap,numVectors) );
  X->setObjectLabel("X");
  Xhat->setObjectLabel("Xhat");
  X->randomize();

  RCP<MV> B = rcp(new MV(rngmap,numVectors));
  A->apply(*X, *B);

  // Constructor from Factory
  RCP<Amesos2::Solver<MAT,MV> > solver;
  try{
    solver = Amesos2::create<MAT,MV>(solver_name, A, Xhat, B);
  } catch (const std::invalid_argument & e){
    *fos << e.what() << std::endl;
    return 0;
  }

  #ifdef HAVE_AMESOS2_SHYLU_NODEBASKER
  if( Amesos2::query("ShyLUBasker") ) {
    Teuchos::ParameterList amesos2_params("Amesos2");
      amesos2_params.sublist(solver_name).set("num_threads", 1, "Number of threads");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  }
  #endif

  solver->numericFactorization();

  if( printLUStats && myRank == 0 ){
    Amesos2::Status solver_status = solver->getStatus();
    *fos << "L+U nnz = " << solver_status.getNnzLU() << std::endl;
  }

  solver->solve();

  if( printSolution ){
    // Print the solution
    Xhat->describe(*fos,Teuchos::VERB_EXTREME);
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  if( printTiming ){
    // Print some timing statistics
    solver->printTiming(*fos);
  }

  if( printResidual ){
    Teuchos::Array<Magnitude> xhatnorms(numVectors);
    Xhat->update(-1.0, *X, 1.0);
    Xhat->norm2(xhatnorms());
    if( myRank == 0 ){
      *fos << "Norm2 of Ax - b = " << xhatnorms << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
