// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Matrix-Market files
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  typedef double Scalar;
  //typedef std::complex<double> Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;
  typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar>> reader_type;

  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Tpetra::Import;
  using Teuchos::RCP;
  using Teuchos::rcp;


  //
  // Get the default communicator
  //
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int myRank = comm->getRank();

  Teuchos::oblackholestream blackhole;

  bool printTiming   = false;
  bool equil         = false;
  bool tinyPivot     = true;
  bool multiSolves   = false;
  bool solveIR       = false;
  int  mc64_job      = 1;
  int  numIRs        = 5;
  bool verboseIR     = true;
  bool verbose       = true;
  std::string solverName("SuperLUDist");
  std::string rowPerm("NOROWPERM");
  std::string filename("arc130.mtx");
  std::string rhsFilename("");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("equil","noEquil",&equil,"Use Equil.");
  cmdp.setOption("rowperm",&rowPerm,"RowPerm.");
  cmdp.setOption("job",&mc64_job,"Option for MC64.");
  cmdp.setOption("tinyPivot","noTinyPivot",&tinyPivot,"Replace tiny pivot.");
  cmdp.setOption("multiSolves","noMultiSolves",&multiSolves,"Do numerical factor and solve twice.");
  cmdp.setOption("solveIR","noSolveIR",&solveIR,"Solve with IR.");
  cmdp.setOption("solver",&solverName,"Solver name");
  cmdp.setOption("filename",&filename,"Filename for Matrix-Market test matrix.");
  cmdp.setOption("rhsFilename",&rhsFilename,"Filename for RHS.");
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  std::ostream& out = ( (verbose && myRank == 0) ? std::cout : blackhole );
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  // Say hello
  out << myRank << " : " << Amesos2::version() << std::endl << std::endl;

  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename, comm);
  out << std::endl << A->description() << std::endl << std::endl;

  // get the maps
  RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

  const size_t numVectors = 1;
  RCP<MV> Xhat = rcp(new MV(rngmap,numVectors));
  RCP<MV> X = rcp(new MV(rngmap,numVectors));
  Xhat->putScalar(1.0);

  // Create random X
  X->randomize();

  // Create B
  RCP<MV> B = rcp(new MV(rngmap,numVectors));
  if (rhsFilename != "") {
    B = reader_type::readDenseFile (rhsFilename, comm, rngmap, false, false);
  } else {
    A->apply (*Xhat, *B);
    //B->putScalar(1.0);
  }

  // Constructor from Factory
  RCP<Amesos2::Solver<MAT,MV> > solver;
  solver = Amesos2::create<MAT,MV>(solverName, A, X, B);
  Teuchos::ParameterList amesos2_params("Amesos2");
  if (solverName == "SuperLUDist") {
    auto superlu_params = Teuchos::sublist(Teuchos::rcpFromRef(amesos2_params), "SuperLU_DIST");
    superlu_params->set("Equil", equil);
    superlu_params->set("RowPerm", rowPerm);
    superlu_params->set("LargeDiag_MC64-Options", mc64_job);
    superlu_params->set("ReplaceTinyPivot", tinyPivot);
  }
  if (solveIR) {
    amesos2_params.set("Iterative refinement", true);
    amesos2_params.set("Number of iterative refinements", numIRs);
    amesos2_params.set("Verboes for iterative refinement", verboseIR);
  }
  solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  solver->symbolicFactorization().numericFactorization().solve();
  if (multiSolves) {
    solver->numericFactorization().solve();
  }

  if (verbose) {
    using mag_type = MV::mag_type;
    MV R (B->getMap (), B->getNumVectors ());
    Teuchos::Array<mag_type> B_norms (B->getNumVectors ());
    Teuchos::Array<mag_type> R_norms (R.getNumVectors ());
    B->norm2 (B_norms ());
    A->apply (*X, R);
    R.update (1.0, *B, -1.0);
    R.norm2 (R_norms ());
    out << "normR = " << R_norms[0] << " / " << B_norms[0] << " = " << R_norms[0]/B_norms[0] << std::endl;
  }

  if( printTiming ){
    // Print some timing statistics
    solver->printTiming(*fos);
  }
  Teuchos::TimeMonitor::summarize();

  // We are done.
  return 0;
}
