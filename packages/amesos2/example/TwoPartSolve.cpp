// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
 * \file   TwoPartSolve.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Wed Jun 22 11:56:12 2011
 *
 * \brief  An example of using Amesos2 to factor a matrix before the X
 *         and B vectors are known, then solving once they are available.
 */

#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Import.hpp>

#include <MatrixMarket_Tpetra.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Tpetra::Import;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;


  //
  // Get the default communicator
  //
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int myRank  = comm->getRank();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  bool printMatrix   = false;
  bool printSolution = false;
  bool printTiming   = false;
  bool verbose       = false;
  std::string filename("arc130.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Matrix-Market test matrix.");
  cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("print-solution","no-print-solution",&printSolution,"Print solution vector after solve.");
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Say hello
  if( myRank == 0 ) *fos << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename,comm);
  if( printMatrix ){
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if( verbose && myRank==0 ){
    *fos << std::endl << A->description() << std::endl << std::endl;
  }

  // We have our matrix, create our solver and factor
  RCP<Amesos2::Solver<MAT,MV> > solver;
  try{
    solver = Amesos2::create<MAT,MV>("Superlu", A);
  } catch(std::invalid_argument e){
    // This solver is not supported/enabled.  This is not really a
    // "failure", so we exit with success.
    return EXIT_SUCCESS;
  }
  solver->symbolicFactorization();
  solver->numericFactorization();

  // Now create X and B vectors

  // get the matrix maps
  RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

  // Create random X
  RCP<MV> X = rcp( new MV(dmnmap,numVectors) );
  RCP<MV> Xhat = rcp(new MV(dmnmap, numVectors));
  X->randomize();

  // Initialize Xhat[0..n] = 10
  Xhat->putScalar(10);

  /* Create B
   *
   * Use RHS:
   *
   *  [[10]
   *   [10]
   *   [10]
   *   [10]
   *   [10]
   *   [10]]
   */
  RCP<MV> B = rcp(new MV(rngmap,numVectors));
  B->putScalar(10);

  try{
    solver->setX(X);
    solver->setB(B);
    solver->solve();

    if( printSolution ){
      // Print the solution
      X->describe(*fos,Teuchos::VERB_EXTREME);
    }

    // Create a new B vector from a give Xhat, and solve with it.
    RCP<MV> B_new = rcp(new MV(rngmap,numVectors));
    A->apply(*Xhat, *B_new);

    if( verbose ){
      if( myRank == 0) *fos << "New RHS vector:" << std::endl;
      B_new->describe(*fos,Teuchos::VERB_EXTREME);
    }

    solver->setB(B_new);
    solver->solve();

    if( printSolution ){
      // Print the solution
      X->describe(*fos,Teuchos::VERB_EXTREME);
    }

    if( printTiming ){
      // Print some timing statistics
      solver->printTiming(*fos);
    }

    if( verbose ){
      Array<Magnitude> xhatnorms(numVectors);
      Xhat->update(-1.0, *X, 1.0);
      Xhat->norm2(xhatnorms());
      *fos << "Norm2 of Ax - b = " << xhatnorms << std::endl;
    }

  } catch (std::invalid_argument e){
    *fos << "The solver does not support the matrix shape" << std::endl;
  }

  // We are done.
  return 0;
}
