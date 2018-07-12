// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//
// @HEADER

/**
   \file   SimpleSolveNonContigMap.cpp
   \author Eric Bavier <etbavie@sandia.gov>
   \author Nathan Ellingwood <ndellin@sandia.gov>
   \date   Wed Apr 26 10:08:39 2017

   \brief  Simple example of Amesos2 usage with a map with non-contiguous GIDs.

   This example solves a simple sparse system of linear equations using the
   Amesos2 interface

   The example is hard-coded for 2 MPI processes
*/

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Teuchos_CommandLineProcessor.hpp>


#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

#ifdef HAVE_AMESOS2_SHYLUBASKER
  std::string solver_name = "ShyLUBasker";
#else
  std::string solver_name = "Basker";
#endif
  Teuchos::CommandLineProcessor cmdp(false,true);
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


  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  size_t myRank = comm->getRank();

  std::ostream &out = std::cout;

  out << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  const global_size_t numProcs = comm->getSize();

  global_size_t nrows = 6;

  if( numProcs != 2 ){
    std::cerr << "Example should be run with number of processes \
      equal to 2 (hard-coded example).  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  const GO numGlobalEntries = nrows;
  const LO numLocalEntries = nrows / numProcs;

  // Create non-contiguous Map
  // This example: np 2 leads to GIDS: proc0 - 0,2,4  proc 1 - 1,3,5
  Teuchos::Array<GO> elementList(numLocalEntries);
  for ( LO k = 0; k < numLocalEntries; ++k ) {
        elementList[k] = myRank + k*numProcs + 4*myRank;
  }

  typedef Tpetra::Map<LO,GO>  map_type;
  RCP< const map_type > map = rcp( new map_type(numGlobalEntries, elementList, 0, comm) );
  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () > 1 && map->isContiguous (),
    std::logic_error,
    "The non-contiguous map claims to be contiguous.");


  //RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row
  RCP<MAT> A = rcp( new MAT(map,0) );

  /*
   * We will solve a system with a known solution, for which we will be using
   * the following matrix:
   *
   *  GID  0   2   4   5   7   9
   * [ 0 [ 7,  0, -3,  0, -1,  0 ]
   *   2 [ 2,  8,  0,  0,  0,  0 ]
   *   4 [ 0,  0,  1,  0,  0,  0 ]
   *   5 [-3,  0,  0,  5,  0,  0 ]
   *   7 [ 0, -1,  0,  0,  4,  0 ]
   *   9 [ 0,  0,  0, -2,  0,  6 ] ]
   *
   */

  // Construct matrix
    if( myRank == 0 ){
      A->insertGlobalValues(0,tuple<GO>(0,4,7),tuple<Scalar>(7,-3,-1));
      A->insertGlobalValues(2,tuple<GO>(0,2),tuple<Scalar>(2,8));
      A->insertGlobalValues(4,tuple<GO>(4),tuple<Scalar>(1));
      A->insertGlobalValues(5,tuple<GO>(0,5),tuple<Scalar>(-3,5));
      A->insertGlobalValues(7,tuple<GO>(2,7),tuple<Scalar>(-1,4));
      A->insertGlobalValues(9,tuple<GO>(5,9),tuple<Scalar>(-2,6));
    }

  A->fillComplete();

  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () > 1 && A->getMap()->isContiguous (),
    std::logic_error,
    "The non-contiguous map of A claims to be contiguous.");


  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->randomize();

  /* Create B, use same GIDs
   *
   * Use RHS:
   *
   *  [[-7]
   *   [18]
   *   [ 3]
   *   [17]
   *   [18]
   *   [28]]
   */
  RCP<MV> B = rcp(new MV(map,numVectors));
  int rowids[6] = {0,2,4,5,7,9};
  int data[6] = {-7,18,3,17,18,28};
  for( int i = 0; i < 6; ++i ){
    if( B->getMap()->isNodeGlobalElement(rowids[i]) ){
      B->replaceGlobalValue(rowids[i],0,data[i]);
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () > 1 && X->getMap()->isContiguous () && B->getMap()->isContiguous (),
    std::logic_error,
    "The non-contiguous maps of X or B claims to be contiguous.");


  // Create solver interface with Amesos2 factory method
  RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>(solver_name, A, X, B);

  if( Amesos2::query(solver_name) ) {

    // Create a Teuchos::ParameterList to hold solver parameters
    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist(solver_name).set("IsContiguous", false, "Are GIDs Contiguous");

//    Teuchos::ParameterList & sublist_params = amesos2_params.sublist(solver_name);
//     sublist_params.set("IsContiguous", false, "Are GIDs Contiguous");
#ifdef HAVE_AMESOS2_SHYLUBASKER
    if ( solver_name == "ShyLUBasker" ) {
      amesos2_params.sublist(solver_name).set("num_threads", 1, "Number of threads");
//      sublist_params.set("num_threads", 1, "Number of threads");
    }
#endif

    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  }

  solver->symbolicFactorization().numericFactorization().solve();


  /* Print the solution
   *
   * Should be:
   *
   *  [[1]
   *   [2]
   *   [3]
   *   [4]
   *   [5]
   *   [6]]
   */
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  *fos << "\nSolution :" << std::endl;
  X->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;

  // We are done.
  return 0;
}
