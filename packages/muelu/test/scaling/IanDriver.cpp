// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <algorithm>  // shuffle
#include <vector>     // vector
#include <random>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_StackedTimer.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
// MueLu
#include "MueLu.hpp"
#include "MueLu_TestHelpers.hpp"
#include <MatrixLoad.hpp>

#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraImport.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "KokkosSparse_spmv.hpp"
#include "Kokkos_Core.hpp"

Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
Teuchos::RCP<Teuchos::TimeMonitor> globalTimeMonitor;

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using std::endl;
  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int numProc                         = comm->getSize();

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> pOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out      = *pOut;
  out.setOutputToRootOnly(0);

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 5, ny = 5, nz = 50;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  bool binaryFormat = false;
  clp.setOption("binary", "ascii", &binaryFormat, "print timings to screen");
  std::string matrixFile;
  clp.setOption("matrixfile", &matrixFile, "matrix market file containing matrix");
  std::string rowMapFile;
  clp.setOption("rowmap", &rowMapFile, "map data file");
  std::string colMapFile;
  clp.setOption("colmap", &colMapFile, "colmap data file");
  std::string domainMapFile;
  clp.setOption("domainmap", &domainMapFile, "domainmap data file");
  std::string rangeMapFile;
  clp.setOption("rangemap", &rangeMapFile, "rangemap data file");

  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  int nrepeat = 1;
  clp.setOption("nrepeat", &nrepeat, "repeat the experiment N times");

  bool describeMatrix = true;
  clp.setOption("showmatrix", "noshowmatrix", &describeMatrix, "describe matrix");
  bool useStackedTimer = false;
  clp.setOption("stackedtimer", "nostackedtimer", &useStackedTimer, "use stacked timer");
  bool verboseModel = false;
  clp.setOption("verbosemodel", "noverbosemodel", &verboseModel, "use stacked verbose performance model");



  std::ostringstream galeriStream;
  std::string rhsFile, coordFile, coordMapFile, nullFile, materialFile;  // unused
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  RCP<RealValuedMultiVector> coordinates;
  RCP<MultiVector> nullspace, material, x, b;
  RCP<Matrix> A;
  RCP<const Map> map;

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  // Load the matrix off disk (or generate it via Galeri), assuming only one right hand side is loaded.
  MatrixLoad<SC, LO, GO, NO>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile,
                             domainMapFile, rangeMapFile, coordFile, coordMapFile, nullFile, materialFile,
                             map, A, coordinates, nullspace, material, x, b, 1,
                             galeriParameters, xpetraParameters, galeriStream);


  out << "========================================================" << endl
      << xpetraParameters
      // << matrixParameters
      << "========================================================" << endl
      << "Template Types:" << endl
      << "  Scalar:        " << Teuchos::demangleName(typeid(SC).name()) << endl
      << "  LocalOrdinal:  " << Teuchos::demangleName(typeid(LO).name()) << endl
      << "  GlobalOrdinal: " << Teuchos::demangleName(typeid(GO).name()) << endl
      << "  Node:          " << Teuchos::demangleName(typeid(NO).name()) << endl
      << "Sizes:" << endl
      << "  Scalar:        " << sizeof(SC) << endl
      << "  LocalOrdinal:  " << sizeof(LO) << endl
      << "  GlobalOrdinal: " << sizeof(GO) << endl
      << "========================================================" << endl
      << "Matrix:        " << Teuchos::demangleName(typeid(Matrix).name()) << endl
      << "Vector:        " << Teuchos::demangleName(typeid(MultiVector).name()) << endl
      << "Hierarchy:     " << Teuchos::demangleName(typeid(Hierarchy).name()) << endl
      << "========================================================" << endl
      << " MPI Ranks:    " << numProc << endl
      << "========================================================" << endl;

  // =========================================================================
  // Problem construction
  // =========================================================================

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
  RCP<const crs_matrix_type> At;
  vector_type xt;
  vector_type yt;

  std::cout << "just before Op2TpetraCrs" << std::endl;

  A->describe(*pOut, Teuchos::VERB_EXTREME);
  //At                         = Utilities::Op2TpetraCrs(A);

}  // main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
