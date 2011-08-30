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

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <MatrixMarket_Tpetra.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef int LO;
  typedef int GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();
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
    return EXIT_SUCCESS;	// Otherwise CTest will pick it up as
				// failure, which it isn't really
  }

  const size_t numVectors = 1;

  // create a Map
  global_size_t nrows = 6;
  RCP<Tpetra::Map<LO,GO> > map
    = rcp( new Tpetra::Map<LO,GO>(nrows,0,comm) );

  std::string mat_pathname = filedir + filename;
  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(mat_pathname,comm,node);

  if( printMatrix ){
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if( verbose && myRank==0 ){
    *fos << std::endl << A->description() << std::endl << std::endl;
  }

  // get the maps
  RCP<const Tpetra::Map<LO,GO,Node> > dmnmap = A->getDomainMap();      	
  RCP<const Tpetra::Map<LO,GO,Node> > rngmap = A->getRangeMap();

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
  } catch (std::invalid_argument e){
    *fos << e.what() << std::endl;
    return 0;
  }

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
