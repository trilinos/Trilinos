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
   \file   SimpleSolve_WithParameters.cpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Sat Jul 17 14:07:36 2010

   \brief  Small example of using a Teuchos::ParameterList to specify Amesos2
           solver parameters.
*/

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

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
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));


  size_t myRank = comm->getRank();

  if( myRank == 0 ) *fos << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  // create a Map
  global_size_t nrows = 6;
  RCP<Tpetra::Map<LO,GO> > map = rcp( new Tpetra::Map<LO,GO>(nrows,0,comm) );
  RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row

  /*
   * We will solve a system with a known solution, for which we will be using
   * the following matrix:
   *
   * [ [ 7,  2,  0, -3,  0,  0 ]
   *   [ 0,  8,  0,  0, -1,  0 ]
   *   [ -3, 0,  1,  0,  0,  0 ]
   *   [ 0,  0,  0,  5,  0, -2 ]
   *   [ -1, 0,  0,  0,  4,  0 ]
   *   [ 0,  0,  0,  0,  0,  6 ] ]
   *
   * And we will solve with A^T
   */
  // Construct matrix
  if( myRank == 0 ){
    A->insertGlobalValues(0,tuple<GO>(0,1,3),tuple<Scalar>(7,2,-3));
    A->insertGlobalValues(1,tuple<GO>(1,4),tuple<Scalar>(8,-1));
    A->insertGlobalValues(2,tuple<GO>(0,2),tuple<Scalar>(-3,1));
    A->insertGlobalValues(3,tuple<GO>(3,5),tuple<Scalar>(5,-2));
    A->insertGlobalValues(4,tuple<GO>(0,4),tuple<Scalar>(-1,4));
    A->insertGlobalValues(5,tuple<GO>(5),tuple<Scalar>(6));
  }
  A->fillComplete();

  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->randomize();

  /* Create B
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
  int data[6] = {-7,18,3,17,18,28};
  for( int i = 0; i < 6; ++i ){
    if( B->getMap()->isNodeGlobalElement(i) ){
      B->replaceGlobalValue(i,0,data[i]);
    }
  }

  // Check first whether SuperLU is supported
  if( Amesos2::query("Superlu") ){

    // Constructor from Factory
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Superlu", A, X, B);

    // Create a Teuchos::ParameterList to hold solver parameters
    Teuchos::ParameterList amesos2_params("Amesos2");
    Teuchos::ParameterList superlu_params = amesos2_params.sublist("SuperLU");
    superlu_params.set("Trans","TRANS","Whether to solve with A^T");
    superlu_params.set("Equil",false,"Whether to equilibrate the system before solve");
    superlu_params.set("ColPerm","NATURAL","Use 'natural' ordering of columns");

    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
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
    X->describe(*fos,Teuchos::VERB_EXTREME);
  } else {
    *fos << "SuperLU solver not enable.  Exiting..." << std::endl;
  }

  // We are done.
  return EXIT_SUCCESS;
}
