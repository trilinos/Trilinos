// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_Import.hpp>

#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>


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
  using Teuchos::RCP;
  using Teuchos::rcp;


  // 
  // Get the default communicator
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();
  int myRank  = comm->getRank();

  Teuchos::oblackholestream blackhole;
  std::ostream &out = ( myRank == 0 ? std::cout : blackhole );
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  bool printMatrix   = false;
  bool printSolution = false;
  bool printTiming   = false;
  bool verbose = (myRank==0);
  std::string filename("bcsstk14.hb");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("print_matrix","no_print_matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("print_solution","no_print_solution",&printSolution,"Print solution vector after solve.");
  cmdp.setOption("print_timing","no_print_timing",&printTiming,"Print solver timing statistics");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Say hello
  out << Amesos::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  RCP<MAT> A;
  Tpetra::Utils::readHBMatrix(filename,comm,node,A);
  if( printMatrix ){
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if( verbose ){
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  // create a Map
  global_size_t nrows = A->getGlobalNumRows();
  RCP<Tpetra::Map<LO,GO,Node> > map = rcp( new Tpetra::Map<LO,GO,Node>(nrows,0,comm) );
  RCP<Map<LO,GO,Node> > root_map
    = rcp( new Map<LO,GO,Node>(nrows,myRank == 0 ? nrows : 0,0,comm) );
  RCP<MV> Xhat = rcp( new MV(root_map,numVectors) );
  RCP<Import<LO,GO,Node> > importer = rcp( new Import<LO,GO,Node>(map,root_map) );

  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
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
  RCP<MV> B = rcp(new MV(map,numVectors));
  B->putScalar(10);
  

  // Constructor from Factory
  RCP<Amesos::SolverBase> solver;
  try {
    solver = Amesos::Factory<MAT,MV>::create("Superlu",A,X,B);

    solver->symbolicFactorization().numericFactorization().solve();

    if( printSolution ){
      // Print the solution
      Xhat->doImport(*X,*importer,Tpetra::REPLACE);
      if( myRank == 0 ){
        *fos << "Solution :" << std::endl;
        Xhat->describe(*fos,Teuchos::VERB_EXTREME);
        *fos << std::endl;
      }
    }

    if( printTiming ){
      // Print some timing statistics
      solver->printTiming(*fos);
    }
  } catch ( std::invalid_argument e ){
    *fos << "Solver does not support this matrix shape" << std::endl;
  }
  
  // We are done.
  return 0;
}
