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

#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

namespace { // (anonymous)

using map_type = Tpetra::Map<>;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;
using Scalar = double;
using MAT = Tpetra::CrsMatrix<Scalar>;
using MV = Tpetra::MultiVector<Scalar>;
using reader_type = Tpetra::MatrixMarket::Reader<MAT>;

bool
readEntryFromFile (GO& gblRowInd, GO& gblColInd, Scalar& val, const std::string& s)
{
  if (s.size () == 0 || s.find ("%") != std::string::npos) {
    return false; // empty line or comment line
  }
  std::istringstream in (s);

  if (! in) {
    return false;
  }
  in >> gblRowInd;
  if (! in) {
    return false;
  }
  in >> gblColInd;
  if (! in) {
    return false;
  }
  in >> val;
  return true;
}

// NOTE: This is not a general routine, despite occasional efforts at
// being so.  It expects a particular input file and only works with a
// particular test, namely the test for Issue #3647.
Teuchos::RCP<MAT>
readCrsMatrixFromFile (const std::string& matrixFilename,
                       const Teuchos::RCP<const map_type>& rowMap,
                       const Teuchos::RCP<const map_type>& domainMap,
                       const Teuchos::RCP<const map_type>& rangeMap)
{
  auto comm = rowMap->getComm ();
  const int myRank = comm->getRank ();

  std::ifstream inFile;
  int opened = 0;
  if (myRank == 0) {
    try {
      inFile.open (matrixFilename);
      if (inFile) {
        opened = 1;
      }
    }
    catch (...) {
      opened = 0;
    }
  }
  Teuchos::broadcast<int, int> (*comm, 0, Teuchos::outArg (opened));
  TEUCHOS_TEST_FOR_EXCEPTION
    (opened == 0, std::runtime_error, "readCrsMatrixFromFile: "
     "Failed to open file \"" << matrixFilename << "\" on Process 0.");

  using Teuchos::RCP;
  RCP<MAT> A (new MAT (rowMap, 0));

  if (myRank == 0) {
    std::string line;

    // Skip the first four lines.  This is a hack, specific to the input file in question.

    // %%MatrixMarket matrix coordinate real general
    // %% the_linear_solver
    // %% Global Tpetra matrix for system: the_linear_solver :: numRows numColumns numEntries = 440 440 9360
    // 440 440 9360
    std::getline (inFile, line);
    std::getline (inFile, line);
    std::getline (inFile, line);
    std::getline (inFile, line);

    while (inFile) {
      std::getline (inFile, line);
      GO gblRowInd {};
      GO gblColInd {};
      Scalar val {};
      const bool gotLine = readEntryFromFile (gblRowInd, gblColInd, val, line);
      if (gotLine) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (! rowMap->isNodeGlobalElement (gblRowInd), std::invalid_argument,
           gblRowInd << " not in row Map");
        A->insertGlobalValues (gblRowInd, LO (1), &val, &gblColInd);
      }
    }
  }

  A->fillComplete (domainMap, rangeMap);
  return A;
}

} // namespace (anonymous)

int main(int argc, char *argv[]) {
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;

  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  {
    auto comm = Tpetra::getDefaultComm();
    TEUCHOS_TEST_FOR_EXCEPTION
      (comm->getSize () != 1, std::runtime_error, "This test must run on "
       "exactly 1 MPI process.");
    const int myRank = comm->getRank();

    bool printTiming   = false;
    bool allprint      = false;
    bool verbose = (myRank==0);
    std::string matrixFilename; // ("arc130.mtx");
    std::string rhsFilename;
    std::string mapFilename;
    std::string solverName;

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("mapFilename", &mapFilename, "Name of Map file");
    cmdp.setOption("matrixFilename", &matrixFilename, "Name of sparse matrix "
                   "file in Matrix Market format, to load as the linear system "
                   "matrix A");
    cmdp.setOption("rhsFilename", &rhsFilename, "Name of dense (multi)vector "
                   "file in Matrix Market format, to load as the linear system "
                   "right-hand side(s) B");
    cmdp.setOption("solverName", &solverName, "Amesos2 solver name");
    cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");
    cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }

    Teuchos::oblackholestream blackhole;
    std::ostream& out = ( (allprint || (myRank == 0)) ? std::cout : blackhole );
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

    // Say hello
    out << myRank << " : " << Amesos2::version() << std::endl << std::endl;

    RCP<MAT> A;
    RCP<const map_type> rowMap;
    RCP<const map_type> domainMap;
    RCP<const map_type> rangeMap;
    if (mapFilename == "") {
      if (verbose) {
        *fos << "Read sparse matrix from file \"" << matrixFilename << "\"" << endl;
      }
      A = reader_type::readSparseFile (matrixFilename, comm);
      domainMap = A->getDomainMap();
      rangeMap = A->getRangeMap();
    }
    else {
      if (verbose) {
        *fos << "Read row Map from file \"" << mapFilename << "\"" << endl;
      }
      rowMap = reader_type::readMapFile (mapFilename, comm);
      // Test readMapFile while we have the chance.
      {
        std::ostringstream outMap;
        Tpetra::MatrixMarket::Writer<MAT>::writeMap (outMap, *rowMap);
        std::istringstream inMap (outMap.str ());
        auto rowMap2 = reader_type::readMap (inMap, comm);

        const bool same = rowMap->isSameAs (*rowMap2);
        TEUCHOS_TEST_FOR_EXCEPTION
          (! same, std::logic_error, "map -> readMapFile -> writeMap -> readMap "
           "does not result in the same Map.");
      }
      //rowMap->describe (*fos, Teuchos::VERB_EXTREME);

      domainMap = rowMap;
      rangeMap = rowMap;
      if (verbose) {
        *fos << "Read sparse matrix from file \"" << matrixFilename << "\"" << endl;
      }
      A = readCrsMatrixFromFile (matrixFilename, rowMap, domainMap, rangeMap);
    }

    RCP<MV> B;
    if (rhsFilename == "") {
      if (verbose) {
        *fos << "Generate right-hand side(s)" << endl;
      }
      const size_t numVectors = 1;
      B = rcp (new MV (rangeMap, numVectors));
      B->putScalar (10.0);
    }
    else {
      if (verbose) {
        *fos << "Read right-hand side(s) from file \"" << rhsFilename << "\"" << endl;
      }
      B = Tpetra::MatrixMarket::Reader<MAT>::readDenseFile (rhsFilename, comm, rangeMap);
    }

    // Create random X
    RCP<MV> X = rcp (new MV (domainMap, B->getNumVectors ()));
    X->randomize();

    // Constructor from Factory
    RCP<Amesos2::Solver<MAT,MV> > solver;

    if (solverName != "" && ! Amesos2::query (solverName)) {
      *fos << "Solver \"" << solverName << "\" not enabled.  Exiting..." << endl;
      return EXIT_SUCCESS;
    }

    if (verbose) {
      *fos << "About to create solver named \"" << solverName << "\"" << endl;
    }
    solver = Amesos2::create<MAT,MV>(solverName, A, X, B);
    TEUCHOS_TEST_FOR_EXCEPTION(solver.is_null (), std::runtime_error, "Uh oh, solver is null");

    if (verbose) {
      *fos << "About to solve linear system" << endl;
    }
    solver->symbolicFactorization().numericFactorization().solve();

    if( printTiming ){
      // Print some timing statistics
      solver->printTiming(*fos);
    }
    Teuchos::TimeMonitor::summarize();
  }
  return 0;
}
