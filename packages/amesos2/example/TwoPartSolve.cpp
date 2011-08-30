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
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Import.hpp>

#include <MatrixMarket_Tpetra.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef int Ordinal;

  typedef double Scalar;
  typedef int LO;
  typedef int GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

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
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();
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

  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename,comm,node);
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
  RCP<const Map<LO,GO,Node> > dmnmap = A->getDomainMap();		
  RCP<const Map<LO,GO,Node> > rngmap = A->getRangeMap();

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
