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

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
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
  std::string solver_name = "SuperLU";
  std::string filename("arc130.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Matrix-Market test matrix.");
  cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("print-solution","no-print-solution",&printSolution,"Print solution vector after solve.");
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");
  cmdp.setOption("solver", &solver_name, "Which TPL solver library to use.");
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

  // get the maps
  RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

  // Create random X
  RCP<MV> X = rcp( new MV(dmnmap,numVectors) );
  X->randomize();

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

  // Constructor from Factory
  RCP<Amesos2::Solver<MAT,MV> > solver;
  try{
    solver = Amesos2::create<MAT,MV>(solver_name, A, X, B);
  } catch (std::invalid_argument e){
    *fos << e.what() << std::endl;
    return 0;
  }

  solver->symbolicFactorization().numericFactorization().solve();

  if( printSolution ){
    // Print the solution
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  // change one of the matrix values and re-solve.
  //
  // Replace the lowest column index and lowest row index entry with "20"
  A->resumeFill();
  A->replaceGlobalValues(Teuchos::as<GO>(A->getRowMap()->getMinGlobalIndex()),
                         tuple<GO>(A->getColMap()->getMinGlobalIndex()),
                         tuple<Scalar>(20));
  A->fillComplete();

  solver->numericFactorization().solve();

  if( printSolution ){
    // Print the solution
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  // change the RHS vector and re-solve.
  B->randomize();
  if( verbose ){
    if( myRank == 0) *fos << "New RHS vector:" << std::endl;
    B->describe(*fos,Teuchos::VERB_EXTREME);
  }

  solver->solve();

  if( printSolution ){
    // Print the solution
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  if( printTiming ){
    // Print some timing statistics
    solver->printTiming(*fos);
  }

  // We are done.
  return 0;
}
