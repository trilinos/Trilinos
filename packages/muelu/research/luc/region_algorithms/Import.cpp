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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#include <MueLu_Utilities.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_TentativePFactory.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Export.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::Array; using Teuchos::ArrayView;
  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::TimeMonitor;

  typedef Tpetra::Map<LO,GO,NO> map_type;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO> crs_matrix_type;

  typedef typename crs_matrix_type::scalar_type scalar_type;
  typedef typename crs_matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename crs_matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename crs_matrix_type::node_type node_type;
  typedef typename crs_matrix_type::mag_type magnitude_type;

  bool success = false;
  bool verbose = true;

  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    const local_ordinal_type myRank = comm->getRank();

    // Manage the way output stream works
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    Teuchos::FancyOStream& out = *fancy;
    // out.setOutputToRootOnly(0);
    // // Do not print anything
    // Teuchos::oblackholestream blackhole;
    // std::ostream& out = blackhole;

    // =========================================================================
    // Parameters initialization
    // =========================================================================

    //GO nx = 100, ny = 100, nz = 100;
    //Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    std::string xmlFileName = "driver.xml";      clp.setOption("xml",                   &xmlFileName,     "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    global_ordinal_type gNumCompFineGIDs = 10, lCompFineGIDOffset = 0;
    global_ordinal_type gNumRegFineGIDs  = 12, lRegFineGIDOffset  = 0;
    local_ordinal_type  lNumCompFineGIDs =  0, lNumRegFineGIDs = 0;
    local_ordinal_type  lNumCompCoarseGIDs =  0, lNumRegCoarseGIDs = 0;
    Array<global_ordinal_type> quasiGIDs, compGIDs;
    if(myRank == 0) {
      compGIDs.resize(3);
      quasiGIDs.resize(3);
      for(int i = 0; i < 3; ++i) {
	quasiGIDs[i] = i;
	compGIDs[i] = i;
      }
    } else if(myRank == 1) {
      quasiGIDs.resize(5);
      quasiGIDs[0] = 2;
      quasiGIDs[1] = 3;
      quasiGIDs[2] = 4;
      quasiGIDs[3] = 4;
      quasiGIDs[4] = 5;

      compGIDs.resize(3);
      for(int i = 0; i < 3; ++i) {
	compGIDs[i] = i + 3;
      }
    }
    out << "quasiGIDs: " << quasiGIDs << std::endl;
    out << "compGIDs:  " << compGIDs << std::endl;
    RCP<const map_type> compMap  = rcp(new map_type(6,  compGIDs, 0, comm));
    RCP<const map_type> quasiMap = rcp(new map_type(8, quasiGIDs, 0, comm));

    Tpetra::MultiVector<SC,LO,GO,NO> quasiVector(quasiMap, 1);
    Teuchos::ArrayRCP<SC> inputVector = quasiVector.getDataNonConst(0);
    if(myRank == 0) {
      for(int i = 0; i < 3; ++i) {
	inputVector[i] = i;
      }
    } else if(myRank == 1) {
      for(int i = 0; i < 3; ++i) {
	inputVector[i] = i + 3;
      }
    }
    Tpetra::MultiVector<SC,LO,GO,NO> compVector(compMap, 1);

    Tpetra::Export<LO,GO,NO> myExporter(quasiMap, compMap);
    compVector.doExport(quasiVector, myExporter, Tpetra::CombineMode::ADD);
    Teuchos::ArrayRCP<const SC> myData = compVector.getData(0);
    std::cout << "p=" << comm->getRank() << " | compVector: " << myData() << std::endl;

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
