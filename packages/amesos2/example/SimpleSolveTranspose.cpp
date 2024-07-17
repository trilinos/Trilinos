// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \brief  Simple example of Amesos2 usage.

   This example solves a simple sparse system of linear equations using the
   Amesos2 interface to the ShyLUBasker solver.
*/

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  typedef Tpetra::CrsMatrix<>::scalar_type Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;


  // Before we do anything, check that ShyLUBasker is enabled
  if( !Amesos2::query("ShyLUBasker") ){
    std::cerr << "ShyLUBasker not enabled.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Tpetra::getDefaultComm();

  size_t myRank = comm->getRank();

  std::ostream &out = std::cout;
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  out << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  // create a Map
  const global_size_t nrows = 3;
  RCP<Tpetra::Map<LO,GO> > map
    = rcp( new Tpetra::Map<LO,GO>(nrows,0,comm) );

  RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row

  /*
   * We will solve a tiny system with a known solution, for which we will be using
   * the following matrix:
   *
   * [ [ 1,  0,  0 ]
   *   [ 1,  2,  0 ]
   *   [ 1,  0,  3 ] ]
   *
   */
  // Construct matrix
  if( myRank == 0 ){
    A->insertGlobalValues(0,tuple<GO>(0),tuple<Scalar>(1));
    A->insertGlobalValues(1,tuple<GO>(0,1),tuple<Scalar>(1,2));
    A->insertGlobalValues(2,tuple<GO>(0,2),tuple<Scalar>(1,3));
  }
  A->fillComplete();
  A->describe(*fos,Teuchos::VERB_EXTREME);

  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->randomize();

  /* Create B
   *
   * Use RHS (A*ones):
   *
   *  [[ 1 ]
   *   [ 3 ]
   *   [ 4 ]]
   */

  RCP<MV> B = rcp(new MV(map,numVectors));
  int data[nrows] = {1,3,4};
  for( global_size_t i = 0; i < nrows; ++i ){
    if( B->getMap()->isNodeGlobalElement(i) ){
      B->replaceGlobalValue(i,0,data[i]);
    }
  }
  B->describe(*fos,Teuchos::VERB_EXTREME);


  /* Create BT
   *
   * Use RHS (A'*ones):
   *
   *  [[ 3 ]
   *   [ 2 ]
   *   [ 3 ]]
   */

  RCP<MV> BT = rcp(new MV(map,numVectors));
  int dataT[nrows] = {3,2,3};
  for( global_size_t i = 0; i < nrows; ++i ){
    if( BT->getMap()->isNodeGlobalElement(i) ){
      BT->replaceGlobalValue(i,0,dataT[i]);
    }
  }

  {
    // Solve A \ B
    // Create solver interface to Superlu with Amesos2 factory method
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("ShyLUBasker", A, X, B);
    *fos << "\nAmesos2: Solver description:  " << solver->description() << std::endl;

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("ShyLUBasker").set("num_threads", 1, "Num threads == 1 by default");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();


    /* Print the solution
     *
     * Should be:
     *
     *  [[1]
     *   [1]
     *   [1]]
     */

    *fos << "\nA^-1 * B Solution :" << std::endl;
    X->describe(*fos,Teuchos::VERB_EXTREME);
    *fos << std::endl;
  }

  {
    // Solve A' \ BT
    // Create solver interface to Superlu with Amesos2 factory method
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("ShyLUBasker", A, X, BT);

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("ShyLUBasker").set("num_threads", 1, "Num threads == 1 by default");
    amesos2_params.sublist("ShyLUBasker").set("transpose", true, "Transpose solve");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();

    *fos << "\nAT^-1 * BT Solution :" << std::endl;
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  {
    // Solve A \ B then Solve A' \ BT - same solver instance
    // Create solver interface to Superlu with Amesos2 factory method
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("ShyLUBasker", A, X, B);

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("ShyLUBasker").set("num_threads", 1, "Num threads == 1 by default");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();

    *fos << "\nA^-1 * B Solution Try2 :" << std::endl;
    X->describe(*fos,Teuchos::VERB_EXTREME);


    solver->setB(BT);
    amesos2_params.set("Transpose", true, "transpose solve");
    //amesos2_params.sublist("ShyLUBasker").set("transpose", true, "Transpose solve");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
//    solver->symbolicFactorization();
//    solver->numericFactorization();
    solver->solve();

    *fos << "\nAT^-1 * BT Solution Try2 :" << std::endl;
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  // We are done.
  return 0;
}
