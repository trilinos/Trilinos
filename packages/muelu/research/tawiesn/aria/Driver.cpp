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
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#include <MueLu_VisualizationHelpers.hpp>
#include <MueLu_FacadeClassFactory.hpp>

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosXpetraStatusTestGenResSubNorm.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ExportVTK : public MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    ExportVTK() {};

  public:

    void writeFile(std::ofstream& fout, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& coordinates, Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& sol)
    {
      using namespace std;
      typedef MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> VH;
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > nodeMap = coordinates->getMap();
      std::vector<LocalOrdinal> vertices;
      std::vector<LocalOrdinal> geomSize;
      LocalOrdinal numFineNodes = Teuchos::as<LocalOrdinal>(coordinates->getLocalLength());

      vertices.reserve(numFineNodes);
      geomSize.reserve(numFineNodes);
      for(LocalOrdinal i = 0; i < numFineNodes; i++)
      {
        vertices.push_back(i);
        geomSize.push_back(1);
      }

      Teuchos::ArrayRCP<const double> xCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(0));
      Teuchos::ArrayRCP<const double> yCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(1));
      Teuchos::ArrayRCP<const double> zCoords = Teuchos::null;
      if(coordinates->getNumVectors() == 3) {
        zCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(2));
      }
      Teuchos::ArrayRCP<const double> solData = Teuchos::arcp_reinterpret_cast<const double>(sol->getData(0));

      std::vector<int> uniqueFine = this->makeUnique(vertices);
      this->writeFileVTKOpening(fout, uniqueFine, geomSize);
      this->writeFileVTKNodes(fout, uniqueFine, nodeMap);

      std::string indent = "      ";
      fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << std::endl;
      fout << indent;
      int myRank = coordinates->getMap()->getComm()->getRank();
      for(int i = 0; i < int(uniqueFine.size()); i++)
      {
        fout << myRank << " ";
        if(i % 20 == 19)
          fout << std::endl << indent;
      }
      fout << std::endl;
      fout << "        </DataArray>" << std::endl;
      // solution vector
      fout << "        <DataArray type=\"Float64\" Name=\"Solution\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      fout << indent;
      for(int i = 0; i < int(uniqueFine.size()); i++)
      {
        size_t numVec = coordinates->getNumVectors();
        for(int d=0; d<numVec; d++)
          fout << solData[i*numVec+d] << " ";
        if(i % 3 == 0)
          fout << std::endl << indent;
      }
      fout << std::endl;
      fout << "        </DataArray>" << std::endl;
      fout << "      </PointData>" << std::endl; // that is annoying

      this->writeFileVTKCoordinates(fout, uniqueFine, xCoords, yCoords, zCoords, coordinates->getNumVectors());
      this->writeFileVTKCells(fout, uniqueFine, vertices, geomSize);
      this->writeFileVTKClosing(fout);
      fout.close();

    };

  };

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& fancyout = *fancy;
    fancyout.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);

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

    // =========================================================================
    // Read in problem
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time"))), tm;



    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main
